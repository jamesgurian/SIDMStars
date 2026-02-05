
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

