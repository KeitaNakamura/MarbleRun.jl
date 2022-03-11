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

function timestep(model::DruckerPrager, pt, dx) # currently support only LinearElastic
    ρ = pt.m / pt.V
    v = norm(pt.v)
    vc = matcalc(Val(:sound_speed), model.elastic.K, model.elastic.G, ρ)
    dx / (vc + v)
end

# https://doi.org/10.1016/j.ijnonlinmec.2011.10.007
function timestep(model::NewtonianFluid, pt, dx)
    ρ = pt.m / pt.V
    ν = model.μ / ρ # kinemtatic viscosity
    v = norm(pt.v)
    vc = model.pressure_model.c # speed of sound
    min(dx / (vc + v), dx^2 / ν)
end

##########################################
# generate_pointstate/initialize_stress! #
##########################################

# Since initializing pointstate is dependent on types of simulation,
# `initialize!` phase should be given as argument.
# This `initialize!` function is called on each material to perform
# material-wise initialization.
# After initialization, these `pointstate`s will be concatenated.
# Then, if you have rigid bodies, the invalid points are removed.
function Poingr.generate_pointstate(initialize!::Function, ::Type{PointState}, grid::Grid, input::Input) where {PointState}
    # generate all pointstate first
    Material = input.Material
    pointstates = map(1:length(Material)) do matindex
        material = Material[matindex]
        pointstate′ = generate_pointstate( # call method in `Poingr`
            material.region,
            PointState,
            grid;
            n = input.Advanced.npoints_in_cell,
        )
        initialize!(pointstate′, matindex)
        pointstate′
    end
    pointstate = first(pointstates)
    append!(pointstate, pointstates[2:end]...)

    remove_invalid_pointstate!(pointstate, input)

    pointstate
end

function remove_invalid_pointstate!(pointstate, input::Input)
    α = input.Advanced.contact_threshold_scale
    !isempty(input.RigidBody) && deleteat!(
        pointstate,
        findall(eachindex(pointstate)) do p
            xₚ = pointstate.x[p]
            rₚ = pointstate.r[p]
            any(map(x -> x.model, input.RigidBody)) do rigidbody
                # remove pointstate which is in rigidbody or is in contact with rigidbody
                in(xₚ, rigidbody) || distance(rigidbody, xₚ, α * mean(rₚ)) !== nothing
            end
        end,
    )
    pointstate
end

# This function is basically based on Material.Initialization
function initialize_stress!(pointstate::AbstractVector, material::Input_Material, g)
    init = material.init
    ρ0 = material.density
    if init isa Input_Material_init_K0
        K0 = init.K0
        for p in eachindex(pointstate)
            σ_y = -ρ0 * g * (init.height_ref - pointstate.x[p][2]) # TODO: handle 3D
            σ_x = K0 * σ_y
            pointstate.σ[p] = (@Mat [σ_x 0.0 0.0
                                     0.0 σ_y 0.0
                                     0.0 0.0 σ_x]) |> symmetric
        end
    elseif init isa Input_Material_init_Uniform
        for p in eachindex(pointstate)
            pointstate.σ[p] = init.mean_stress * one(pointstate.σ[p])
        end
    else
        error("unreachable")
    end
end

####################
# advance timestep #
####################

function advancestep!(grid::Grid{dim}, pointstate::AbstractVector, rigidbodies, cache, dt, input::Input, phase::Input_Phase) where {dim}
    g = input.General.gravity
    materials = input.Material
    matmodels = map(x -> x.model, materials)

    if isempty(rigidbodies)
        update!(cache, grid, pointstate)
    else
        spatmasks = [falses(size(grid)) for i in 1:Threads.nthreads()]
        Threads.@threads for rigidbody in rigidbodies
            mask = spatmasks[Threads.threadid()]
            broadcast!(in(rigidbody), mask, grid)
        end
        exclude = broadcast(|, spatmasks...)
        update!(cache, grid, pointstate; exclude)
    end

    # Point-to-grid transfer
    P2G!(grid, pointstate, cache, dt, input)

    # Compute contact forces between rigid bodies
    contact_forces = fill(zero(Vec{dim}), length(rigidbodies))
    contact_moments = fill(zero(Vec{3}), length(rigidbodies))
    compute_contact_forces_moments!(contact_forces, contact_moments, rigidbodies, grid, dt, input)

    # Point-to-grid transfer for contact and update positions of rigid bodies
    for (i, rigidbody, input_rigidbody) in zip(1:length(rigidbodies), rigidbodies, input.RigidBody)
        α = input.Advanced.contact_threshold_scale
        ξ = input.Advanced.contact_penalty_parameter
        mask = P2G_contact!(grid, pointstate, cache, dt, rigidbody, input_rigidbody.frictions, α, ξ)
        if phase.update_motion == true
            # Update rigid bodies
            b = input_rigidbody.body_force
            if input_rigidbody.control
                GeometricObjects.update!(rigidbody, b, zero(Vec{3}), dt)
            else
                G2P_contact!(pointstate, grid, cache, rigidbody, mask)
                inds = findall(mask)
                Fc, Mc = GeometricObjects.compute_force_moment(rigidbody, view(pointstate.fc, inds), view(pointstate.x, inds))
                Fc += rigidbody.m * Vec(0,-g) + b
                GeometricObjects.update!(rigidbody, contact_forces[i]+Fc, contact_moments[i]+Mc, dt)
            end
        end
    end

    # Boundary conditions
    for dirichlet in input.BoundaryCondition.Dirichlet
        dirichlet.displacement += norm(dirichlet.velocity) * dt
        dirichlet.reaction_force = 0.0
    end
    @inbounds for bound in eachboundary(grid)
        I = bound.I
        vᵢ = grid.state.v[I]
        v = vᵢ
        isdirichlet = false
        for dirichlet in input.BoundaryCondition.Dirichlet
            if dirichlet.active_nodes[I]
                vᵢ′ = dirichlet.velocity
                dirichlet.reaction_force += grid.state.m[I] * (norm(vᵢ-vᵢ′) / dt)
                v = vᵢ′
                isdirichlet = true
            end
        end
        if !isdirichlet
            v = boundary_velocity(vᵢ, bound.n, input.BoundaryCondition)
        end
        grid.state.v[I] = v
    end

    # Grid-to-point transfer
    G2P!(pointstate, grid, cache, matmodels, dt, input, phase)
end

##########################
# point-to-grid transfer #
##########################

function P2G!(grid::Grid, pointstate::AbstractVector, cache::MPCache, dt::Real, input::Input)
    input.General.transfer.point_to_grid!(grid, pointstate, cache, dt)
end

function P2G_contact!(grid::Grid, pointstate::AbstractVector, cache::MPCache, dt::Real, rigidbody::GeometricObject, frictions, α::Real, ξ::Real)
    mask = @. distance($Ref(rigidbody), pointstate.x, α * mean(pointstate.r)) !== nothing
    !any(mask) && return mask
    point_to_grid!((grid.state.d, grid.state.vᵣ, grid.state.μ, grid.state.m_contacted), cache, mask) do it, p, i
        @_inline_meta
        @_propagate_inbounds_meta
        N = it.N
        mₚ = pointstate.m[p]
        xₚ = pointstate.x[p]
        vₚ = pointstate.v[p]
        m = N * mₚ
        d₀ = α * mean(pointstate.r[p]) # threshold
        coef = frictions[pointstate.matindex[p]]
        if length(coef) == 1
            μ = only(coef)
            dₚ = distance(rigidbody, xₚ, d₀)
        else
            # friction is interpolated
            dₚ, μ = distance(rigidbody, xₚ, d₀, coef)
        end
        d = d₀*normalize(dₚ) - dₚ
        vᵣ = vₚ - velocityat(rigidbody, xₚ)
        m*d, m*vᵣ, m*μ, m
    end
    @dot_threads grid.state.m′ = grid.state.m / (grid.state.m/rigidbody.m + 1)
    @dot_threads grid.state.d /= grid.state.m
    @dot_threads grid.state.vᵣ /= grid.state.m_contacted
    @dot_threads grid.state.μ /= grid.state.m_contacted
    @dot_threads grid.state.fc = compute_contact_force(grid.state.d, grid.state.vᵣ, grid.state.m′, dt, grid.state.μ, ξ)
    @dot_threads grid.state.v += (grid.state.fc / grid.state.m) * dt
    mask
end

function compute_contact_force(d::Vec{dim, T}, vᵣ::Vec{dim, T}, m::T, dt::T, μ::Vec{2, T}, ξ::T) where {dim, T}
    iszero(d) && return zero(Vec{dim, T})
    n = normalize(d)
    vᵣ_tan = vᵣ - (vᵣ ⋅ n) * n
    f_nor = (1-ξ) * 2m/dt^2 * d
    f_tan = m * (vᵣ_tan / dt)
    Contact(:friction, μ[1]; sep = true, thresh = μ[2])(f_nor + f_tan, n)
end

##########################
# grid-to-point transfer #
##########################

function G2P!(pointstate::AbstractVector, grid::Grid, cache::MPCache, models::Vector{<: MaterialModel}, dt::Real, input::Input, phase::Input_Phase)
    input.General.transfer.grid_to_point!(pointstate, grid, cache, dt)
    @inbounds Threads.@threads for p in eachindex(pointstate)
        matindex = pointstate.matindex[p]
        model = models[matindex]
        σ, dϵ, F, J = compute_σ_dϵ_F_J(model, pointstate.σ[p], pointstate.∇v[p], pointstate.F[p], pointstate.J[p], dt)
        pointstate.σ[p] = σ
        pointstate.ϵ[p] += dϵ
        pointstate.F[p] = F
        pointstate.J[p] = J
        pointstate.V[p] = J * pointstate.V0[p]
        if phase.update_motion == false
            pointstate.x[p] -= pointstate.v[p] * dt
            pointstate.v[p] = zero(pointstate.v[p])
            pointstate.∇v[p] = zero(pointstate.∇v[p])
        end
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

function compute_σ_dϵ_F_J(model::VonMises, σ_n::SymmetricSecondOrderTensor, ∇v::SecondOrderTensor, F::SecondOrderTensor, J::Real, dt::Real)
    dt∇v = dt * ∇v
    F = F + dt∇v ⋅ F
    dϵ = symmetric(dt∇v) # rate of deformation tensor
    σ = matcalc(Val(:stress), model, σ_n, dϵ)
    σ = matcalc(Val(:jaumann_stress), σ, σ_n, ∇v, dt)
    σ, dϵ, F, det(F)
end

function compute_σ_dϵ_F_J(model::DruckerPrager, σ_n::SymmetricSecondOrderTensor, ∇v::SecondOrderTensor, F::SecondOrderTensor, J::Real, dt::Real)
    dt∇v = dt * ∇v
    dϵ = symmetric(dt∇v) # rate of deformation tensor
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
        # modify velocity gradient from re-calculated strain tensor (symmetric component)
        dt∇v = dϵ + skew(dt∇v)
    end
    F = F + dt∇v ⋅ F
    σ, dϵ, F, det(F)
end

function compute_σ_dϵ_F_J(model::NewtonianFluid, σ_n::SymmetricSecondOrderTensor, ∇v::SecondOrderTensor, F::SecondOrderTensor, J::Real, dt::Real)
    J = (1 + tr(∇v)*dt) * J # direct linear integration for `J` is more accurate than using `F` (maybe?)
    σ = matcalc(Val(:stress), model, ∇v, J, dt)
    # σ = matcalc(Val(:jaumann_stress), σ, σ_n, ∇v, dt)
    σ, zero(σ), F, J # NOTE: `dϵ` is not used
end

######################
# boundary condition #
######################

function boundary_velocity(v::Vec{2}, n::Vec{2}, bc::Input_BoundaryCondition)
    n == Vec(-1,  0) && (side = :left)
    n == Vec( 1,  0) && (side = :right)
    n == Vec( 0, -1) && (side = :bottom)
    n == Vec( 0,  1) && (side = :top)
    contact = getproperty(bc, side)
    v + contact(v, n)
end

##################################
# contact forces/moments for DEM #
##################################

function compute_contact_forces_moments!(contact_forces::Vector, contact_moments::Vector, rigidbodies, grid, dt, input)
    @assert length(contact_forces) == length(contact_moments) == length(rigidbodies) == length(input.RigidBody)
    @inbounds for (i, rigidbody1) in enumerate(rigidbodies)
        input.RigidBody[i].control && continue
        for rigidbody2 in rigidbodies
            if rigidbody1 !== rigidbody2
                vals = compute_contactforce_position(rigidbody1, rigidbody2, dt, input)
                if vals !== nothing
                    fc, x = vals
                    contact_forces[i] += fc
                    contact_moments[i] += (x - centroid(rigidbody1)) × fc
                end
            end
        end
        vals = compute_contactforce_position(rigidbody1, grid, dt, input)
        if vals !== nothing
            fc, x = vals
            contact_forces[i] += fc
            contact_moments[i] += (x - centroid(rigidbody1)) × fc
        end
    end
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

function quickview_sparsity_pattern(spat::AbstractMatrix{Bool}; maxwidth::Int = 50, maxheight::Int = 25)
    spat′ = reverse(spat', dims = 1)
    spy(spat′; maxwidth, maxheight).graphics
end
