

# Define abstract types for the categories
abstract type ConvectionMode end
abstract type EvaporationMode end

# Define singleton types for the specific options
struct ConvectionOn <: ConvectionMode end
struct ConvectionOff <: ConvectionMode end

struct EvaporationOn <: EvaporationMode end
struct EvaporationOff <: EvaporationMode end
@kwdef struct RunParams{TConv <: ConvectionMode, TEvap <: EvaporationMode}
 
    convection::TConv = ConvectionOn()
    evaporation::TEvap = EvaporationOn()
    G::Float64
    C::Float64 
    σ0::Float64
    w::Float64
    n::Float64
    t_end::Float64
    r_max::Float64
    r_min::Float64
    gridpoints::Int64
    snap_interval_time::Float64
    snap_interval_density::Float64
    snap_loc::String
end

struct CompositeODESolution{T, U, S1, S2}
    t::T
    u::U
    sol_out::S1
    sol_in::S2
    t_fit::Float64
    potential_shift::Float64
end

function CompositeODESolution(sol_out, sol_in; t_fit = sol_out.t[end])
    # We want the zero point defined by sol_out, so we shift sol_in to match at t_fit.
    u_out_fit = sol_out(t_fit)
    u_in_fit  = sol_in(t_fit)
    shift = u_out_fit[3] - u_in_fit[3]
    u_in_shifted = map(sol_in.u) do u
        u_new = copy(u)
        u_new[3] += shift
        u_new
    end

    # Inward integrations have decreasing time arrays. We reverse them 
    # and drop the duplicate fitting point at index 1 of the reversed array.
    t_combined = vcat(sol_out.t, reverse(sol_in.t)[2:end])
    u_combined = vcat(sol_out.u, reverse(u_in_shifted)[2:end])

    # Wrap in VectorOfArray to enable sol.u[2, :] matrix-like slicing
    u_voa = VectorOfArray(u_combined)

    return CompositeODESolution(t_combined, u_voa, sol_out, sol_in, t_fit, shift)
end

# --- DifferentialEquations.jl API Compatibility ---

function (csol::CompositeODESolution)(t::Number)
    if t <= csol.t_fit
        return csol.sol_out(t)
    else
        u_val = csol.sol_in(t)
        return [u_val[1], u_val[2], u_val[3] + csol.potential_shift]
    end
end

function (csol::CompositeODESolution)(t_vec::AbstractArray)
    u_evaluated = [csol(t) for t in t_vec]
    return DiffEqArray(u_evaluated, t_vec)
end

# 1. Scalar time input for derivatives
function (csol::CompositeODESolution)(t::Number, ::Type{Val{N}}) where N
    if t <= csol.t_fit
        return csol.sol_out(t, Val{N})
    else
        return csol.sol_in(t, Val{N})
    end
end

# 2. Vector time input for derivatives
function (csol::CompositeODESolution)(t_vec::AbstractArray, ::Type{Val{N}}) where N
    # Evaluate the derivative interpolant at each time point
    u_evaluated = [csol(t, Val{N}) for t in t_vec]
    
    # Wrap in DiffEqArray so the derivative output also supports matrix slicing
    return DiffEqArray(u_evaluated, t_vec)
end

Base.getindex(csol::CompositeODESolution, i, j) = csol.u[i, j]
Base.getindex(csol::CompositeODESolution, i) = csol.u[i]
Base.length(csol::CompositeODESolution) = length(csol.t)
Base.firstindex(csol::CompositeODESolution) = 1
Base.lastindex(csol::CompositeODESolution) = length(csol.t)