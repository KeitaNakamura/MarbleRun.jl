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

    # remove invalid pointstate
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
function initialize_stress!(σₚ::AbstractVector, material::Input_Material, g)
    init = material.init
    ρ0 = material.density
    if init isa Input_Material_init_K0
        K0 = init.K0
        for p in eachindex(σₚ)
            σ_y = -ρ0 * g * init.height_ref
            σ_x = K0 * σ_y
            σₚ[p] = (@Mat [σ_x 0.0 0.0
                           0.0 σ_y 0.0
                           0.0 0.0 σ_x]) |> symmetric
        end
    elseif init isa Input_Material_init_Uniform
        for p in eachindex(σₚ)
            σₚ[p] = init.mean_stress * one(σₚ[p])
        end
    else
        error("unreachable")
    end
end

####################
# advance timestep #
####################

function advancestep!(grid::Grid, pointstate::AbstractVector, rigidbodies, cache, dt, input::Input, phase::Input_Phase)
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
    P2G!(grid, pointstate, cache, dt)
    for (i, rigidbody) in enumerate(rigidbodies)
        α = input.Advanced.contact_threshold_scale
        ξ = input.Advanced.contact_penalty_parameter
        mask = P2G_contact!(grid, pointstate, cache, dt, rigidbody, i, materials, α, ξ)
        if phase.update_motion == true
            # Update rigid bodies
            input_rigidbody = input.RigidBody[i]
            b = input_rigidbody.body_force
            if input_rigidbody.control
                GeometricObjects.update!(rigidbody, b, zero(Vec{3}), dt)
            else
                G2P_contact!(pointstate, grid, cache, rigidbody, mask)
                inds = findall(mask)
                Fc, Mc = GeometricObjects.compute_force_moment(rigidbody, view(pointstate.fc, inds), view(pointstate.x, inds))
                Fc += rigidbody.m * Vec(0,-g) + b
                GeometricObjects.update!(rigidbody, Fc, Mc, dt)
            end
        end
    end

    # Boundary conditions
    for bd in eachboundary(grid)
        @inbounds grid.state.v[bd.I] = boundary_velocity(grid.state.v[bd.I], bd.n, input.BoundaryCondition)
    end

    # Grid-to-point transfer
    G2P!(pointstate, grid, cache, matmodels, materials, dt, phase) # `materials` are for densities
end

##########################
# point-to-grid transfer #
##########################

function P2G!(grid::Grid, pointstate::AbstractVector, cache::MPCache, dt::Real)
    default_point_to_grid!(grid, pointstate, cache, dt)
end

function P2G_contact!(grid::Grid, pointstate::AbstractVector, cache::MPCache, dt::Real, rigidbody::GeometricObject, rigidbody_index::Int, materials, α::Real, ξ::Real)
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
        mat = materials[pointstate.matindex[p]]
        coef = mat.friction_with_rigidbodies[rigidbody_index]
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

function G2P!(pointstate::AbstractVector, grid::Grid, cache::MPCache, models::Vector{<: MaterialModel}, materials::Vector{<: Union{Input_Material, Input_SoilLayer}}, dt::Real, phase::Input_Phase)
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

function boundary_velocity(v::Vec{2}, n::Vec{2}, bc::Input_BoundaryCondition)
    n == Vec(-1,  0) && (side = :left)
    n == Vec( 1,  0) && (side = :right)
    n == Vec( 0, -1) && (side = :bottom)
    n == Vec( 0,  1) && (side = :top)
    contact = getproperty(bc, side)
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
