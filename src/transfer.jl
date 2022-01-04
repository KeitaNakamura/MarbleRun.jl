function compute_contact_force(d::Vec{dim, T}, vᵣ::Vec{dim, T}, m::T, dt::T, μ::T, ξ::T) where {dim, T}
    iszero(d) && return zero(Vec{dim, T})
    n = normalize(d)
    vᵣ_tan = vᵣ - (vᵣ ⋅ n) * n
    f_nor = (1-ξ) * 2m/dt^2 * d
    f_tan = m * (vᵣ_tan / dt)
    Contact(:friction, μ; sep = true)(f_nor + f_tan, n)
end

function compute_σ_dϵ(model::VonMises, σ_n::SymmetricSecondOrderTensor, ∇v::SecondOrderTensor, dt::Real)
    dϵ = symmetric(∇v) * dt
    σ = matcalc(Val(:stress), model, σ_n, dϵ)
    σ = matcalc(Val(:jaumann_stress), σ, σ_n, ∇v, dt)
    σ, dϵ
end

function compute_σ_dϵ(model::DruckerPrager, σ_n::SymmetricSecondOrderTensor, ∇v::SecondOrderTensor, dt::Real)
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
    σ, dϵ
end

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
end

function G2P!(pointstate::AbstractVector, grid::Grid, cache::MPCache, models::Vector{<: MaterialModel}, dt::Real)
    default_grid_to_point!(pointstate, grid, cache, dt)
    @inbounds Threads.@threads for p in eachindex(pointstate)
        model = models[pointstate.matindex[p]]
        σ, dϵ = compute_σ_dϵ(model, pointstate.σ[p], pointstate.∇v[p], dt)
        pointstate.σ[p] = σ
        pointstate.ϵ[p] += dϵ
        pointstate.V[p] *= exp(tr(dϵ))
    end
end

function write_vtk_points(vtk, pointstate::AbstractVector)
    ϵ = pointstate.ϵ
    vtk["velocity"] = pointstate.v
    vtk["mean stress"] = @dot_lazy -mean(pointstate.σ)
    vtk["deviatoric stress"] = @dot_lazy deviatoric_stress(pointstate.σ)
    vtk["volumetric strain"] = @dot_lazy volumetric_strain(ϵ)
    vtk["deviatoric strain"] = @dot_lazy deviatoric_strain(ϵ)
    vtk["stress"] = @dot_lazy -pointstate.σ
    vtk["strain"] = ϵ
    vtk["density"] = @dot_lazy pointstate.m / pointstate.V
    vtk["material index"] = pointstate.matindex
end
