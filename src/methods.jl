############
# timestep #
############

function safe_minimum(f, iter)
    ret = typemax(f(first(iter)))
    @inbounds @simd for x in iter
        y = f(x)
        if !(isnan(y) || isinf(y))
            ret = min(ret, y)
        end
    end
    ret
end

function timestep(model::DruckerPrager, p, dx) # currently support only LinearElastic
    ρ = p.m / p.V
    v = norm(p.v)
    vc = matcalc(Val(:sound_speed), model.elastic.K, model.elastic.G, ρ)
    dx / (vc + min(v, vc)) # set limit of `v` as `vc` to prevent using too high velocity
end

# https://doi.org/10.1016/j.ijnonlinmec.2011.10.007
function timestep(model::NewtonianFluid, p, dx)
    ρ = p.m / p.V
    ν = model.μ / ρ # kinemtatic viscosity
    v = norm(p.v)
    vc = model.pressure_model.c # speed of sound
    min(dx/(vc+v), dx^2/ν)
end

####################
# advance timestep #
####################

function advancestep!(grid::Grid, pointstate::AbstractVector, rigidbodies, cache, INPUT, dt)
    g = INPUT.General.gravity
    materials = INPUT.Material
    matmodels = map(create_materialmodel, materials)

    update!(cache, grid, pointstate)

    # Point-to-grid transfer
    P2G!(grid, pointstate, cache, dt)
    masks = map(rigidbodies) do rigidbody
        α = INPUT.Advanced.contact_threshold_scale
        ξ = INPUT.Advanced.contact_penalty_parameter
        P2G_contact!(grid, pointstate, cache, dt, rigidbody, α, ξ)
    end

    # Boundary conditions
    boundary_contacts = create_boundary_contacts(INPUT.BoundaryCondition)
    for bd in eachboundary(grid)
        @inbounds grid.state.v[bd.I] = boundary_velocity(grid.state.v[bd.I], bd.n, boundary_contacts)
    end

    # Grid-to-point transfer
    G2P!(pointstate, grid, cache, matmodels, materials, dt) # `materials` are for densities
    for (i, (rigidbody, mask)) in enumerate(zip(rigidbodies, masks))
        input = INPUT.RigidBody[i]
        if !getoftype(input, :control, true) # rigid bodies don't move freely by default
            G2P_contact!(pointstate, grid, cache, rigidbody, mask)
        end
    end

    # Update rigid bodies
    for (i, (rigidbody, mask)) in enumerate(zip(rigidbodies, masks))
        input = INPUT.RigidBody[i]
        b = getoftype(input, :body_force, zero(Vec{2}))
        if getoftype(input, :control, true) # rigid bodies don't move freely by default
            GeometricObjects.update!(rigidbody, b, zero(Vec{3}), dt)
        else
            inds = findall(mask)
            Fc, Mc = GeometricObjects.compute_force_moment(rigidbody, view(pointstate.fc, inds), view(pointstate.x, inds))
            Fc += rigidbody.m * Vec(0,-g) + b
            GeometricObjects.update!(rigidbody, Fc, Mc, dt)
        end
    end
end

##########################
# point-to-grid transfer #
##########################

function P2G!(grid::Grid, pointstate::AbstractVector, cache::MPCache, dt::Real)
    default_point_to_grid!(grid, pointstate, cache, dt)
end

function P2G_contact!(grid::Grid, pointstate::AbstractVector, cache::MPCache, dt::Real, rigidbody::GeometricObject, α::Real, ξ::Real)
    mask = @. distance($Ref(rigidbody), pointstate.x, α * mean(pointstate.r)) !== nothing
    point_to_grid!((grid.state.d, grid.state.vᵣ, grid.state.μ, grid.state.m_contacted), cache, mask) do it, p, i
        @_inline_meta
        @_propagate_inbounds_meta
        N = it.N
        mₚ = pointstate.m[p]
        xₚ = pointstate.x[p]
        vₚ = pointstate.v[p]
        m = N * mₚ
        d₀ = α * mean(pointstate.r[p]) # threshold
        if length(pointstate.μ[p]) == 1
            μ = only(pointstate.μ[p])
            dₚ = distance(rigidbody, xₚ, d₀)
        else
            # friction is interpolated
            dₚ, μ = distance(rigidbody, xₚ, d₀, pointstate.μ[p])
        end
        d = d₀*normalize(dₚ) - dₚ
        vᵣ = vₚ - velocityat(rigidbody, xₚ)
        m*d, m*vᵣ, m*μ, m
    end
    mᵢ = @dot_lazy grid.state.m / (grid.state.m/rigidbody.m + 1)
    @dot_threads grid.state.d /= grid.state.m
    @dot_threads grid.state.vᵣ /= grid.state.m_contacted
    @dot_threads grid.state.μ /= grid.state.m_contacted
    @dot_threads grid.state.fc = compute_contact_force(grid.state.d, grid.state.vᵣ, grid.state.m, dt, grid.state.μ, ξ)
    @dot_threads grid.state.v += (grid.state.fc / grid.state.m) * dt
    mask
end

function compute_contact_force(d::Vec{dim, T}, vᵣ::Vec{dim, T}, m::T, dt::T, μ::T, ξ::T) where {dim, T}
    iszero(d) && return zero(Vec{dim, T})
    n = normalize(d)
    vᵣ_tan = vᵣ - (vᵣ ⋅ n) * n
    f_nor = (1-ξ) * 2m/dt^2 * d
    f_tan = m * (vᵣ_tan / dt)
    Contact(:friction, μ; sep = true)(f_nor + f_tan, n)
end

##########################
# grid-to-point transfer #
##########################

function G2P!(pointstate::AbstractVector, grid::Grid, cache::MPCache, models::Vector{<: MaterialModel}, materials::Vector{<: InputMaterial}, dt::Real)
    default_grid_to_point!(pointstate, grid, cache, dt)
    @inbounds Threads.@threads for p in eachindex(pointstate)
        matindex = pointstate.matindex[p]
        model = models[matindex]
        params = materials[matindex]
        ρ0 = params.density
        V0 = pointstate.m[p] / ρ0
        J = pointstate.V[p] / V0 # not updated jacobian
        σ, dϵ, J = compute_σ_dϵ_J(model, pointstate.σ[p], pointstate.∇v[p], J, dt)
        pointstate.σ[p] = σ
        pointstate.ϵ[p] += dϵ
        pointstate.V[p] = V0 * J
    end
end

# compute contact force at points
function G2P_contact!(pointstate::AbstractVector, grid::Grid, cache::MPCache, rigidbody::GeometricObject, mask::AbstractVector{Bool})
    Poingr.fillzero!(pointstate.fc)
    grid_to_point!(pointstate.fc, cache, mask) do it, I, p
        @_inline_meta
        @_propagate_inbounds_meta
        fc = grid.state.fc[I]
        -fc * it.N * pointstate.m[p] / grid.state.m_contacted[I]
    end
end

######################
# stress integration #
######################

function compute_σ_dϵ_J(model::VonMises, σ_n::SymmetricSecondOrderTensor, ∇v::SecondOrderTensor, J::Real, dt::Real)
    dϵ = symmetric(∇v) * dt
    σ = matcalc(Val(:stress), model, σ_n, dϵ)
    σ = matcalc(Val(:jaumann_stress), σ, σ_n, ∇v, dt)
    σ, dϵ, J*exp(tr(dϵ))
end

function compute_σ_dϵ_J(model::DruckerPrager, σ_n::SymmetricSecondOrderTensor, ∇v::SecondOrderTensor, J::Real, dt::Real)
    dϵ = symmetric(∇v) * dt
    σ = matcalc(Val(:stress), model, σ_n, dϵ)
    σ = matcalc(Val(:jaumann_stress), σ, σ_n, ∇v, dt)
    if mean(σ) > model.tension_cutoff
        # In this case, since the soil particles are not contacted with
        # each other, soils should not act as continuum.
        # This means that the deformation based on the contitutitive model
        # no longer occurs.
        # So, in this process, we just calculate the elastic strain to keep
        # the consistency with the stress which is on the edge of the yield
        # function, and ignore the plastic strain to prevent excessive generation.
        # If we include this plastic strain, the volume of the material points
        # will continue to increase unexpectedly.
        σ_tr = matcalc(Val(:stress), model.elastic, σ_n, dϵ)
        σ = Poingr.tension_cutoff(model, σ_tr)
        dϵ = model.elastic.Dinv ⊡ (σ - σ_n)
    end
    σ, dϵ, J*exp(tr(dϵ))
end

function compute_σ_dϵ_J(model::NewtonianFluid, σ_n::SymmetricSecondOrderTensor, ∇v::SecondOrderTensor, J::Real, dt::Real)
    J = (1 + tr(∇v)*dt) * J
    σ = matcalc(Val(:stress), model, ∇v, J, dt)
    σ, zero(σ), J # NOTE: dϵ is used
end

######################
# boundary condition #
######################

function boundary_velocity(v::Vec{2}, n::Vec{2}, boundary_contacts)
    n == Vec(-1,  0) && (side = :left)
    n == Vec( 1,  0) && (side = :right)
    n == Vec( 0, -1) && (side = :bottom)
    n == Vec( 0,  1) && (side = :top)
    contact = boundary_contacts[side]
    v + contact(v, n)
end

##########
# output #
##########

function write_vtk_points(vtk, pointstate::AbstractVector)
    σ = pointstate.σ
    ϵ = pointstate.ϵ
    vtk["velocity"] = pointstate.v
    vtk["mean stress"] = @dot_lazy mean(σ)
    vtk["pressure"] = @dot_lazy -mean(σ)
    vtk["deviatoric stress"] = @dot_lazy deviatoric_stress(σ)
    vtk["volumetric strain"] = @dot_lazy volumetric_strain(ϵ)
    vtk["deviatoric strain"] = @dot_lazy deviatoric_strain(ϵ)
    vtk["stress"] = σ
    vtk["strain"] = ϵ
    vtk["density"] = @dot_lazy pointstate.m / pointstate.V
    vtk["material index"] = pointstate.matindex
end
