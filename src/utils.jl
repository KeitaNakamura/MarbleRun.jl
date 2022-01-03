dict2namedtuple(x) = x
dict2namedtuple(x::Dict) = (; (Symbol(key) => dict2namedtuple(value) for (key, value) in x)...)

function parseinput(dict::Dict)
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

########################
# create_materialmodel #
########################

function create_materialmodel(mat::NamedTuple, coordinate_system)
    create_materialmodel(first(mat), Base.tail(mat), coordinate_system)
end

function create_materialmodel(::Type{DruckerPrager}, params, coordinate_system)
    E = params.youngs_modulus
    ν = params.poissons_ratio
    c = params.cohesion
    ϕ = params.friction_angle
    ψ = params.dilatancy_angle
    tension_cutoff = params.tension_cutoff
    elastic = LinearElastic(; E, ν)
    if coordinate_system isa PlaneStrain
        DruckerPrager(elastic, :plane_strain; c, ϕ, ψ, tension_cutoff)
    else
        DruckerPrager(elastic, :circumscribed; c, ϕ, ψ, tension_cutoff)
    end
end

######################
# initialize_stress! #
######################

function initialize_stress!(σₚ::AbstractVector, material::NamedTuple, g)
    Initialization = material.Initialization
    ρ0 = material.density
    if Initialization.type == "K0"
        for p in eachindex(σₚ)
            σ_y = -ρ0 * g * Initialization.reference_height
            σ_x = Initialization.K0 * σ_y
            σₚ[p] = (@Mat [σ_x 0.0 0.0
                           0.0 σ_y 0.0
                           0.0 0.0 σ_x]) |> symmetric
        end
    else
        throw(ArgumentError("wrong initialization type, got $(condition.type)"))
    end
end

###########
# Outputs #
###########

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
