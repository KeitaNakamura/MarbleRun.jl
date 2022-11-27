####################
# grid/point state #
####################

function merge_namedtuple_type(::Type{NamedTuple{names1, types1}}, ::Type{NamedTuple{names2, types2}}) where {names1, types1, names2, types2}
    NamedTuple{(names1..., names2...), Tuple{types1.parameters..., types2.parameters...}}
end

function gridstate_type(input::Input, ::Val{dim}, ::Type{T}) where {dim, T}
    GridState = gridstate_type_base(Val(dim), T)
    if !isempty(input.RigidBody)
        GridState = merge_namedtuple_type(GridState, gridstate_type_contact(Val(dim), T))
    end
    if input.General.v_p_formulation
        GridState = merge_namedtuple_type(GridState, gridstate_type_vpform(Val(dim), T))
    end
    GridState
end

function gridstate_type_base(::Val{dim}, ::Type{T}) where {dim, T}
    @NamedTuple begin
        m   :: T
        v   :: Vec{dim, T}
        v_n :: Vec{dim, T}
    end
end
function gridstate_type_contact(::Val{dim}, ::Type{T}) where {dim, T}
    @NamedTuple begin
        m′          :: T
        m_contacted :: T
        vᵣ          :: Vec{dim, T}
        fc          :: Vec{dim, T}
        d           :: Vec{dim, T}
        μ           :: Vec{dim, T} # [μ, c]
    end
end
function gridstate_type_vpform(::Val{dim}, ::Type{T}) where {dim, T}
    @NamedTuple begin
        poly_coef :: Vec{3, T}
        poly_mat  :: Mat{3, 3, T, 9}
    end
end

function pointstate_type(input::Input, ::Val{dim}, ::Type{T}) where {dim, T}
    PointState = pointstate_type_base(Val(dim), T)
    if input.General.transfer isa APIC
        PointState = merge_namedtuple_type(PointState, pointstate_type_affine(Val(dim), T))
    end
    if input.General.transfer isa Marble.WLSTransfer
        PointState = merge_namedtuple_type(PointState, pointstate_type_linwls(Val(dim), T))
    end
    if !isempty(input.RigidBody)
        PointState = merge_namedtuple_type(PointState, pointstate_type_contact(Val(dim), T))
    end
    if input.General.v_p_formulation
        PointState = merge_namedtuple_type(PointState, pointstate_type_vpform(Val(dim), T))
    end
    PointState
end

function pointstate_type_base(::Val{dim}, ::Type{T}) where {dim, T}
    @NamedTuple begin
        m        :: T
        V        :: T
        x        :: Vec{dim, T}
        x0       :: Vec{dim, T}
        v        :: Vec{dim, T}
        b        :: Vec{dim, T}
        σ        :: SymmetricSecondOrderTensor{3, T, 6}
        ϵ        :: SymmetricSecondOrderTensor{3, T, 6}
        ∇v       :: SecondOrderTensor{3, T, 9}
        r        :: Vec{dim, T}
        index    :: Int
        matindex :: Int
    end
end
function pointstate_type_affine(::Val{dim}, ::Type{T}) where {dim, T}
    @NamedTuple begin
        B :: Mat{dim, dim, T, dim*dim}
    end
end
function pointstate_type_linwls(::Val{dim}, ::Type{T}) where {dim, T}
    @NamedTuple begin
        C :: Mat{dim, dim+1, T, dim*(dim+1)}
    end
end
function pointstate_type_contact(::Val{dim}, ::Type{T}) where {dim, T}
    @NamedTuple begin
        fc :: Vec{dim, T}
    end
end
function pointstate_type_vpform(::Val{dim}, ::Type{T}) where {dim, T}
    @NamedTuple begin
        P :: T
    end
end

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
    vc = @matcalc(:soundspeed; model.elastic.K, model.elastic.G, ρ)
    dx / (vc + v)
end

# https://doi.org/10.1016/j.ijnonlinmec.2011.10.007
function timestep(model::NewtonianFluid, pt, dx)
    ρ = pt.m / pt.V
    ν = model.μ / ρ # kinemtatic viscosity
    v = norm(pt.v)
    vc = model.eos.c # speed of sound
    min(dx / (vc + v), dx^2 / ν)
end

##############################################
# generate_pointstate/initialize_pointstate! #
##############################################

# Since initializing pointstate is dependent on types of simulation,
# `initialize!` phase should be given as argument.
# This `initialize!` function is called on each material to perform
# material-wise initialization.
# After initialization, these `pointstate`s will be concatenated.
# Then, if you have rigid bodies, the invalid points are removed.
function Marble.generate_pointstate(initialize!::Function, ::Type{PointState}, grid::Grid, input::Input) where {PointState}
    # generate all pointstate first
    Material = input.Material
    pointstates = map(1:length(Material)) do matindex
        material = Material[matindex]
        if material.region isa Function
            pointstate′ = generate_pointstate( # call method in `Marble`
                material.region,
                PointState,
                grid;
                n = input.Advanced.npoints_in_cell,
                random = input.Advanced.random_points_generation,
            )
        else # material.region isa String
            pointstate′ = generate_pointstate(PointState, deserialize(joinpath(input.project, material.region)).pointstate)
        end
        initialize!(pointstate′, matindex)
        pointstate′
    end
    pointstate = first(pointstates)
    append!(pointstate, pointstates[2:end]...)

    remove_invalid_pointstate!(pointstate, input)

    pointstate
end

function remove_invalid_pointstate!(pointstate, input::Input)
    α = input.Advanced.contact_threshold_scale_for_initial_state
    !isempty(input.RigidBody) && deleteat!(
        pointstate,
        findall(eachindex(pointstate)) do p
            xₚ = pointstate.x[p]
            rₚ = pointstate.r[p]
            any(input.RigidBody) do input_rigidbody
                rigidbody = input_rigidbody.model
                inverse = input_rigidbody.inverse
                isinbody = in(xₚ, rigidbody)
                # remove pointstate which is inside of rigidbody or is in contact with rigidbody
                (inverse ? !isinbody : isinbody) || distance(rigidbody, xₚ, α * mean(rₚ)) !== nothing
            end
        end,
    )
    pointstate
end

function initialize_pointstate!(pointstate::AbstractVector, material::Input_Material, g)
    init_stress_mass!(pointstate, material.model, material.init, g)
    @. pointstate.x0 = pointstate.x
    @. pointstate.b = Vec(0.0, -g)
end

function init_stress_mass!(pointstate::AbstractVector, ::MaterialModel, init::InitK0, g)
    ρ0 = init.density
    for p in eachindex(pointstate)
        σ_y = -ρ0 * g * (init.height_ref - pointstate.x[p][2]) # TODO: handle 3D
        σ_x = init.K0 * σ_y
        pointstate.σ[p] = (@Mat [σ_x 0.0 0.0
                                 0.0 σ_y 0.0
                                 0.0 0.0 σ_x]) |> symmetric
        pointstate.m[p] = ρ0 * pointstate.V[p]
    end
end

function init_stress_mass!(pointstate::AbstractVector, ::MaterialModel, init::InitUniform, g)
    ρ0 = init.density
    for p in eachindex(pointstate)
        pointstate.σ[p] = init.mean_stress * one(pointstate.σ[p])
        pointstate.m[p] = ρ0 * pointstate.V[p]
    end
end

function init_stress_mass!(pointstate::AbstractVector, model::NewtonianFluid, init::InitUniform, g)
    if init.density !== nothing
        @warn "`density` in `Material.init.Uniform` is overwritten based on the `mean_stress` for `NewtonianFluid` model"
    end
    for p in eachindex(pointstate)
        pointstate.σ[p] = init.mean_stress * one(pointstate.σ[p])
        pointstate.m[p] = @matcalc(:density, model; p = -init.mean_stress) * pointstate.V[p]
    end
end

####################
# advance timestep #
####################

function advancestep!(grid::Grid{dim}, gridstate::AbstractArray{<: Any, dim}, pointstate::AbstractVector, rigidbodies, space::MPSpace, dt, input::Input, phase::Input_Phase) where {dim}
    g = input.General.gravity
    materials = input.Material
    matmodels = map(x -> x.model, materials)

    if isempty(rigidbodies)
        update!(space, pointstate)
    else
        spatmasks = [falses(size(grid)) for i in 1:Threads.nthreads()]
        Threads.@threads for rigidbody in rigidbodies
            mask = spatmasks[Threads.threadid()]
            broadcast!(in(rigidbody), mask, grid)
        end
        exclude = broadcast(|, spatmasks...)
        update!(space, pointstate; exclude)
    end
    update_sparsity_pattern!(gridstate, space)

    # Point-to-grid transfer
    P2G!(gridstate, pointstate, space, dt, input)

    # Compute contact forces between rigid bodies
    contact_forces = fill(zero(Vec{dim}), length(rigidbodies))
    contact_moments = fill(zero(Vec{3}), length(rigidbodies))
    compute_contact_forces_moments!(contact_forces, contact_moments, rigidbodies, grid, dt, input)

    # Point-to-grid transfer for contact and update positions of rigid bodies
    for (i, rigidbody, input_rigidbody) in zip(1:length(rigidbodies), rigidbodies, input.RigidBody)
        α = input.Advanced.contact_threshold_scale
        ξ = input.Advanced.contact_penalty_parameter
        mask = P2G_contact!(gridstate, pointstate, space, dt, rigidbody, input_rigidbody.FrictionWithMaterial, α, ξ)
        if phase.update_motion == true
            # Update rigid bodies
            b = input_rigidbody.body_force
            if input_rigidbody.control
                GeometricObjects.update!(rigidbody, b, zero(Vec{3}), dt)
            else
                G2P_contact!(pointstate, gridstate, space, rigidbody, mask)
                inds = findall(mask)
                Fc, Mc = GeometricObjects.compute_force_moment(rigidbody, view(pointstate.fc, inds), view(pointstate.x, inds))
                Fc += rigidbody.m * Vec(0,-g) + b
                GeometricObjects.update!(rigidbody, contact_forces[i]+Fc, contact_moments[i]+Mc, dt)
            end
        end
    end

    # Boundary conditions
    # for dirichlet boundary condition, Coulomb frictional contact cannot be used for now.
    # the given nodal velocity is directly applied.
    for dirichlet in input.BoundaryCondition.Dirichlet
        fᵢ′ = 0.0
        @inbounds for I in dirichlet.node_indices
            vᵢ = gridstate.v[I]
            vᵢ′ = dirichlet.velocity
            fᵢ′ += gridstate.m[I] * (norm(vᵢ-vᵢ′) / dt)
            gridstate.v[I] = vᵢ′
        end
        dirichlet.displacement += norm(dirichlet.velocity) * dt
        dirichlet.reaction_force = fᵢ′
    end
    for (side, cond) in input.BoundaryCondition.sides
        @inbounds for (I,n) in gridbounds(grid, side)
            gridstate.v[I] += contacted(cond, gridstate.v[I], n)
        end
    end

    # Grid-to-point transfer
    G2P!(pointstate, gridstate, space, matmodels, dt, input, phase)
end

##########################
# point-to-grid transfer #
##########################

function P2G!(gridstate::AbstractArray, pointstate::AbstractVector, space::MPSpace, dt::Real, input::Input)
    point_to_grid!(input.General.transfer, gridstate, pointstate, space, dt)
end

function P2G_contact!(gridstate::AbstractArray, pointstate::AbstractVector, space::MPSpace, dt::Real, rigidbody::GeometricObject, frictions, α::Real, ξ::Real)
    mask = @. distance($Ref(rigidbody), pointstate.x, α * mean(pointstate.r)) !== nothing
    !any(mask) && return mask
    point_to_grid!((gridstate.d, gridstate.vᵣ, gridstate.μ, gridstate.m_contacted), space, mask) do it, p, i
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
    @dot_threads gridstate.m′ = gridstate.m / (gridstate.m/rigidbody.m + 1)
    @dot_threads gridstate.d /= gridstate.m
    @dot_threads gridstate.vᵣ /= gridstate.m_contacted
    @dot_threads gridstate.μ /= gridstate.m_contacted
    @dot_threads gridstate.fc = compute_contact_force(gridstate.d, gridstate.vᵣ, gridstate.m′, dt, gridstate.μ, ξ, $(prod(gridsteps(space.grid))))
    @dot_threads gridstate.v += (gridstate.fc / gridstate.m) * dt
    mask
end

function compute_contact_force(d::Vec{dim, T}, vᵣ::Vec{dim, T}, m::T, dt::T, μ::Vec{2, T}, ξ::T, A::T = one(T)) where {dim, T}
    iszero(d) && return zero(Vec{dim, T})
    n = normalize(d)
    vᵣ_tan = vᵣ - (vᵣ ⋅ n) * n
    f_nor = (1-ξ) * 2m/dt^2 * d
    f_tan = m * (vᵣ_tan / dt)
    # calculate contact force by using stress (F/A) to handle cohesion properly
    contacted(CoulombFriction(; μ = μ[1], c = μ[2], separation = true), (f_nor+f_tan)/A, -n) * A
end

##########################
# grid-to-point transfer #
##########################

function G2P!(pointstate::AbstractVector, gridstate::AbstractArray, space::MPSpace, models::Vector{<: MaterialModel}, dt::Real, input::Input, phase::Input_Phase)
    grid_to_point!(input.General.transfer, pointstate, gridstate, space, dt)
    @inbounds Threads.@threads for p in eachindex(pointstate)
        matindex = pointstate.matindex[p]
        model = models[matindex]
        # currently we don't have constitutive models using deformation gradient `F`.
        # so we calculate the volume `V` from the equation `Vₙ₊₁ = Vₙ * exp(tr(∇v*dt))`
        # instead of using `det(F)` with linear integration of `F` (`Fₙ₊₁ = Fₙ + dt*∇v ⋅ Fₙ`).
        # the equation can be derived by integrating the definition of the rate of
        # the Jacobian `J̇ = J * tr(∇v)`.
        pt = LazyRow(pointstate, p)
        σ, dϵ = compute_σ_dϵ(model, pt, pointstate.∇v[p], dt)
        pointstate.σ[p] = σ
        pointstate.ϵ[p] += dϵ
        pointstate.V[p] *= exp(tr(dϵ))
        if phase.update_motion == false
            pointstate.x[p] -= pointstate.v[p] * dt
            pointstate.v[p] = zero(pointstate.v[p])
            pointstate.∇v[p] = zero(pointstate.∇v[p])
        end
    end

    if input.General.v_p_formulation
        @dot_threads pointstate.P = -mean(pointstate.σ)
        Marble.smooth_pointstate!(pointstate.P, pointstate.V, gridstate, space)
        @inbounds Threads.@threads for p in eachindex(pointstate)
            matindex = pointstate.matindex[p]
            model = models[matindex]
            P = pointstate.P[p]
            σ = pointstate.σ[p]
            ρ = @matcalc(:density, model; p = P)
            pointstate.σ[p] = -P*I + dev(σ)
            pointstate.V[p] = pointstate.m[p] / ρ
        end
    end
end

# compute contact force at points
function G2P_contact!(pointstate::AbstractVector, gridstate::AbstractArray, space::MPSpace, rigidbody::GeometricObject, mask::AbstractVector{Bool})
    Marble.fillzero!(pointstate.fc)
    grid_to_point!(pointstate.fc, space, mask) do it, I, p
        @_inline_meta
        @_propagate_inbounds_meta
        fc = gridstate.fc[I]
        -fc * it.N * pointstate.m[p] / gridstate.m_contacted[I]
    end
end

######################
# stress integration #
######################

@inline function compute_σ_dϵ(model::VonMises, pt, ∇v::SecondOrderTensor{3}, dt::Real)
    @_propagate_inbounds_meta
    σ_n = pt.σ
    dt∇v = dt * ∇v
    dϵ = symmetric(dt∇v)
    σ = @matcalc(:stress, model; σ_n, dϵ)
    dσᴶ = σ - σ_n
    σ = σ_n + @matcalc(:jaumann2caucy; dσ_jaumann=dσᴶ, σ=σ_n, W=skew(dt∇v))
    σ, dϵ
end

@inline function compute_σ_dϵ(model::DruckerPrager, pt, ∇v::SecondOrderTensor{3}, dt::Real)
    @_propagate_inbounds_meta
    σ_n = pt.σ
    dt∇v = dt * ∇v
    dϵ = symmetric(dt∇v)
    ret = @matcalc(:stressall, model; σ=σ_n, dϵ)
    dσᴶ = ret.σ - σ_n
    σ = σ_n + @matcalc(:jaumann2caucy; dσ_jaumann=dσᴶ, σ=σ_n, W=skew(dt∇v))
    if ret.status.tensioncollapse
        # In this case, since the soil particles are not contacted with
        # each other, soils should not act as continuum.
        # This means that the deformation based on the contitutitive model
        # no longer occurs.
        # So, in this process, we just calculate the elastic strain to keep
        # the consistency with the stress which is on the edge of the yield
        # function, and ignore the plastic strain to prevent excessive generation.
        # If we include this plastic strain, the volume of the material points
        # will continue to increase unexpectedly.
        dϵ = @matcalc(:strain, model.elastic; σ=σ-σ_n)
    end
    σ, dϵ
end

@inline function compute_σ_dϵ(model::NewtonianFluid, pt, ∇v::SecondOrderTensor{3}, dt::Real)
    @_propagate_inbounds_meta
    dt∇v = dt * ∇v
    dϵ = symmetric(dt∇v)
    V = pt.V * exp(tr(dϵ)) # need updated volume
    σ = @matcalc(:stress, model; d = symmetric(∇v), ρ = pt.m/V)
    σ, dϵ
end

######################
# boundary condition #
######################

function boundary_velocity(v::Vec{2}, n::Vec{2}, bc::Input_BoundaryCondition)
    n == Vec( 1,  0) && (side = :left)
    n == Vec(-1,  0) && (side = :right)
    n == Vec( 0,  1) && (side = :bottom)
    n == Vec( 0, -1) && (side = :top)
    cond = getproperty(bc, side)
    v + contacted(cond, v, n)
end

##################################
# contact forces/moments for DEM #
##################################

function compute_contact_forces_moments!(contact_forces::Vector, contact_moments::Vector, rigidbodies, grid::Grid, dt, input)
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

# determine dimensions from second argument
vorticity(∇v::Mat{3,3}, ::Vec{2}) = ∇v[2,1] - ∇v[1,2]
vorticity(∇v::Mat{3,3}, ::Vec{3}) = Vec(∇v[3,2]-∇v[2,3], ∇v[1,3]-∇v[3,1], ∇v[2,1]-∇v[1,2])
function writevtk(vtk, list::Input_Paraview_PointState, pointstate::AbstractVector)
    σ = pointstate.σ
    ϵ = pointstate.ϵ
    list.displacement      && (vtk["displacement"]      = @dot_lazy pointstate.x - pointstate.x0)
    list.velocity          && (vtk["velocity"]          = pointstate.v)
    list.mean_stress       && (vtk["mean stress"]       = @dot_lazy mean(σ))
    list.pressure          && (vtk["pressure"]          = @dot_lazy -mean(σ))
    list.deviatoric_stress && (vtk["deviatoric stress"] = @dot_lazy sqrt(3/2 * dev(σ) ⊡ dev(σ)))
    list.volumetric_strain && (vtk["volumetric strain"] = @dot_lazy tr(ϵ))
    list.deviatoric_stress && (vtk["deviatoric strain"] = @dot_lazy sqrt(2/3 * dev(ϵ) ⊡ dev(ϵ)))
    list.stress            && (vtk["stress"]            = σ)
    list.strain            && (vtk["strain"]            = ϵ)
    list.vorticity         && (vtk["vorticity"]         = @dot_lazy vorticity(pointstate.∇v, pointstate.x))
    list.density           && (vtk["density"]           = @dot_lazy pointstate.m / pointstate.V)
    list.material_index    && (vtk["material index"]    = pointstate.matindex)
    vtk
end

function writevtk(vtk, list::Input_Paraview_GridState, gridstate::AbstractArray)
    list.velocity         && (vtk["velocity"]         = vec(gridstate.v))
    list.contact_force    && (vtk["contact_force"]    = vec(gridstate.fc))
    list.contact_distance && (vtk["contact_distance"] = vec(gridstate.d))
    vtk
end

function quickview_sparsity_pattern(sppat::AbstractMatrix{Bool}; maxwidth::Int = 50, maxheight::Int = 25)
    sppat′ = reverse(sppat', dims = 1)
    spy(sppat′; maxwidth, maxheight).graphics
end
