
@kwdef struct RunParams
    C::Float64 
    σ0::Float64
    w::Float64
    n::Float64
    t_end::Float64
    r_max::Float64
    r_min::Float64
    snaptimes::Vector{Float64}
    snap_loc::String
end


struct EvaporationParams
    # Interpolations for background structure: Input is logm, Output is Float64
    logP::FunctionWrapper{Float64, Tuple{Float64}} 
    logr::FunctionWrapper{Float64, Tuple{Float64}} 
    dPhi::FunctionWrapper{Float64, Tuple{Float64}} 
    # Entropy profile: Input is logr, Output is Float64
    s::FunctionWrapper{Float64, Tuple{Float64}}   

    # Nested parameters and constants
    run_params::RunParams
    G::Float64
    v2_i::Float64
end


struct StructureParams
    s::FunctionWrapper{Float64, Tuple{Float64}}
    G::Float64
end

Base.eltype(::Type{<:StructureParams}) = Float64
