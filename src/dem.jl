function _compute_distance_threshold_pointposition(p1::Polygon{dim, T}, p2::Polygon{dim, T}, reverse::Bool) where {dim, T}
    xc = centroid(p1)
    pos = zero(Vec{dim, T})
    dist = zero(Vec{dim, T})
    count = 0
    @inbounds for i in eachindex(p1)
        xₚ = p1[i]
        xisinp = xₚ in p2
        if (!reverse && xisinp) || (reverse && !xisinp)
            d² = T(Inf)
            dₚ = zero(Vec{dim, T})
            for j in eachindex(p2)
                line = GeometricObjects.getline(p2, j)
                d = distance(line, xₚ, T(Inf))
                d === nothing && continue
                dd = dot(d, d)
                d², dₚ = ifelse(dd < d², (dd, d), (d², dₚ))
            end
            pos += xₚ
            dist += dₚ
            count += 1
        end
    end
    count == 0 && return nothing
    d₀ = zero(T)
    (dist/count, d₀, pos/count)
end

function compute_distance_threshold_pointposition(p1::Polygon{dim, T}, p2::Polygon{dim, T}; reverse1::Bool = false, reverse2::Bool = false) where {dim, T}
    res = _compute_distance_threshold_pointposition(p1, p2, reverse2)
    res !== nothing && return (res..., T(1))
    res = _compute_distance_threshold_pointposition(p2, p1, reverse1)
    res !== nothing && return (res..., -T(1))
    nothing
end

function compute_distance_threshold_pointposition(circle::Circle{dim, T}, poly::Polygon{dim, T}) where {dim, T}
    dₚ = distance(poly, centroid(circle), radius(circle))
    dₚ === nothing && return nothing
    d₀ = radius(circle)
    xₚ = centroid(circle)
    dₚ, d₀, xₚ, T(1)
end

function compute_distance_threshold_pointposition(poly::Polygon{dim, T}, circle::Circle{dim, T}) where {dim, T}
    vals = compute_distance_threshold_pointposition(circle, poly)
    vals === nothing && return nothing
    dₚ, d₀, xₚ = vals
    dₚ, d₀, xₚ, -T(1)
end

function compute_distance_threshold_pointposition(c1::Circle{dim, T}, c2::Circle{dim, T}) where {dim, T}
    dₚ = distance(c2, centroid(c1), radius(c1))
    dₚ === nothing && return nothing
    d₀ = radius(c1)
    xₚ = centroid(c1)
    dₚ, d₀, xₚ, T(1)
end

function _compute_contactforce_position(
        input::Input,
        X::GeometricObject{dim, T},
        Y::GeometricObject{dim, T},
        dt::T,
        dₚ::Vec{dim, T},
        d₀::T,
        xₚ::Vec{dim, T},
        sign::T,
    ) where {dim, T}

    d = d₀*normalize(dₚ) - dₚ
    x = xₚ + dₚ + d/2
    vᵣ = velocityat(X, x) - velocityat(Y, x)
    if isinf(X.m)
        m = Y.m
    elseif isinf(Y.m)
        m = X.m
    else
        m = X.m * Y.m / (X.m + Y.m)
    end
    μ = zero(Vec{2, T})
    ξ = T(input.Advanced.dem_contact_penalty_parameter)
    fc = compute_contact_force(d, vᵣ, m, dt, μ, ξ)

    sign*fc, x
end

function compute_contactforce_position(X::GeometricObject, Y::GeometricObject, dt::Real, input)
    vals = compute_distance_threshold_pointposition(X[], Y[])
    vals === nothing && return nothing
    _compute_contactforce_position(input, X, Y, dt, vals...)
end

# for boundaries
function compute_contactforce_position(X::GeometricObject{<: Any, <: Any, <: Polygon}, Y::Grid{2}, dt::Real, input)
    frame = GeometricObject(Polygon(Y[1,1], Y[end,1], Y[end,end], Y[1,end]))
    frame.m = Inf
    vals = compute_distance_threshold_pointposition(X[], frame[]; reverse1 = false, reverse2 = true)
    vals === nothing && return nothing
    _compute_contactforce_position(input, X, frame, dt, vals...)
end
function compute_contactforce_position(X::GeometricObject{<: Any, <: Any, <: Circle}, Y::Grid{2}, dt::Real, input)
    frame = GeometricObject(Polygon(Y[1,1], Y[end,1], Y[end,end], Y[1,end]))
    frame.m = Inf
    vals = compute_distance_threshold_pointposition(X[], frame[])
    vals === nothing && return nothing
    _compute_contactforce_position(input, X, frame, dt, vals...)
end
