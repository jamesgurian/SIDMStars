

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
    # LMFP conductivity calibration constant, measured jointly with the
    # crossover bump below (fit_kappa_local.jl).
    C::Float64 = 0.8620
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
    collapse_density_factor::Float64 = Inf
    collapse_mass_frac::Float64 = 0.0
    # Conductivity crossover: log-Gaussian resistivity bump on the sum of the
    # branch resistivities (the single form; see calc_κ in transport.jl),
    #   1/κ = (Kinv_smfp + Kinv_lmfp) · [1 + A exp(−u²/2W²)],
    #   u = ln(Kinv_smfp/Kinv_lmfp).
    # The bump dies super-polynomially, so both tails are exactly the
    # asymptotic branches (SMFP side Burnett-consistent O(Kn_loc²)).
    # Defaults are the LOCAL values measured from the DSMC heating
    # trajectories (fit_kappa_local.jl, jointly with C above):
    # (A, W, C) = (1.2046 ± 0.039, 1.7499 ± 0.050, 0.8620 ± 0.026).
    bump_A::Float64 = 1.2046
    bump_W::Float64 = 1.7499
    # K_p prefactors (Eq. K_p): K_p = K_p_prefactor * v2^(n/2), with
    # K_p ∝ Γ((n+p+3)/2), for the power-law cross section σ = σ0 (v/w)^n.
    # K5 reproduces the previous `f_k5_prefactor` exactly.
    K5_prefactor::Float64 = ((1/3)  * 2^(n-1) * gamma(n/2 + 4)) / (w^n)
    K7_prefactor::Float64 = ((1/12) * 2^(n-1) * gamma(n/2 + 5)) / (w^n)
    K9_prefactor::Float64 = ((1/60) * 2^(n-1) * gamma(n/2 + 6)) / (w^n)
    # Derived constants for the log-space calc_κ hot path. All K_p share the
    # v2^(n/2) scaling, so K_eff^(2)(v2) = Keff2_prefactor · v2^(n/2) exactly
    # (Eq. k_secondorder evaluated on the prefactors), and the branch
    # resistivities are pure power laws:
    #   Kinv_smfp = exp(lci_s + (n−1)/2 · ln v²)            b = 25√π/32
    #   Kinv_lmfp = exp(lci_l − ln ρ − (n+3)/2 · ln v²)     a = 4/√π
    Keff2_prefactor::Float64 = (28*K5_prefactor^2 + 80*K5_prefactor*K9_prefactor -
                                64*K7_prefactor^2) /
                               (77*K5_prefactor - 112*K7_prefactor + 80*K9_prefactor)
    lci_s::Float64 = log(σ0 * Keff2_prefactor * 32 / (25 * sqrt(pi)))
    lci_l::Float64 = log(pi^1.5 * G / (C * σ0 * K5_prefactor))
    inv2W2::Float64 = 1 / (2 * bump_W^2)
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
        u_new = SVector{3}(u[1], u[2], u[3] + shift)
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