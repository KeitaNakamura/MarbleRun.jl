function parseinput(dict::Dict)
    dict2namedtuple(x::Dict) = (; (Symbol(key) => value for (key, value) in x)...)
    list = map(collect(keys(dict))) do section
        content = dict[section]
        if content isa Dict
            Symbol(section) => dict2namedtuple(content)
        elseif content isa Vector
            Symbol(section) => map(dict2namedtuple, content)
        else
            error("unreachable")
        end
    end
    (; list...)
end

function compute_contact_force(f_nor::Vec{dim, T}, vᵣ::Vec{dim, T}, m::T, dt::T, μ::T) where {dim, T}
    iszero(f_nor) && return zero(Vec{dim, T})
    n = f_nor / norm(f_nor)
    vᵣ_tan = vᵣ - (vᵣ ⋅ n) * n
    f_tan = m * (vᵣ_tan / dt)
    Contact(:friction, μ, sep = true)(f_nor + f_tan, n)
end

# don't need to check if `d` is `nothing`?
function compute_contact_force_normal(d::Vec{dim, T}, d₀::T, m::T, dt::T, ξ::T) where {dim, T}
    norm_d = norm(d)
    n = d / norm_d
    k_nor = (1-ξ) * 2m/dt^2
    -k_nor * (norm_d - d₀) * n
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
    default_point_to_grid!(grid, pointstate, cache)
    @dot_threads grid.state.v += (grid.state.f / grid.state.m) * dt
end

function P2G_contact!(grid::Grid, pointstate::AbstractVector, cache::MPCache, dt::Real, rigidbody::Polygon, v_rigidbody::Vec, α::Real, ξ::Real)
    mask = @. distance($Ref(rigidbody), pointstate.x, α * mean(pointstate.r)) !== nothing
    point_to_grid!((grid.state.fc_nor, grid.state.vᵣ, grid.state.μ, grid.state.w_rigidbody), cache, mask) do it, p, i
        @_inline_meta
        @_propagate_inbounds_meta
        N = it.N
        w = it.w
        xₚ = pointstate.x[p]
        vₚ = pointstate.v[p]
        d₀ = α * mean(pointstate.r[p]) # threshold
        if length(pointstate.μ[p]) == 1
            μ = only(pointstate.μ[p])
            d = distance(rigidbody, xₚ, d₀)
        else
            # friction is interpolated
            d, μ = distance(rigidbody, xₚ, d₀, pointstate.μ[p])
        end
        fc_nor = compute_contact_force_normal(d, d₀, pointstate.m[p], dt, ξ)
        vᵣ = vₚ - v_rigidbody
        N*fc_nor, w*vᵣ, w*μ, w
    end
    @dot_threads grid.state.vᵣ /= grid.state.w_rigidbody
    @dot_threads grid.state.μ /= grid.state.w_rigidbody
    @dot_threads grid.state.fc = compute_contact_force(grid.state.fc_nor, grid.state.vᵣ, grid.state.m, dt, grid.state.μ)
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
