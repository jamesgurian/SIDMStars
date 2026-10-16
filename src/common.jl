

# Define abstract types for the categories
abstract type ConvectionMode end
abstract type EvaporationMode end

# Define singleton types for the specific options
"""
    ConvectionOn()
    ConvectionOff()

Convection mode for [`RunParams`](@ref) (`convection` keyword). Default: `ConvectionOn()`.
"""
struct ConvectionOn <: ConvectionMode end
struct ConvectionOff <: ConvectionMode end

"""
    EvaporationOn()
    EvaporationOff()

Evaporation mode for [`RunParams`](@ref) (`evaporation` keyword). Default: `EvaporationOn()`.
"""
struct EvaporationOn <: EvaporationMode end
struct EvaporationOff <: EvaporationMode end

# --- Cross sections ----------------------------------------------------------
# σ(v) is carried into RunParams by dispatch. InterpolatedCS handles an arbitrary
# shape σ̂(v²) = σ(v)/σ0 by tabulating the K5, K7 and K_eff^(2) moments on a
# uniform grid in ln v² at construction.
"""
    CrossSection

Abstract supertype for velocity-dependent cross sections passed to
[`RunParams`](@ref) via the `cs` keyword. Concrete types: [`PowerLawCS`](@ref),
[`InterpolatedCS`](@ref) (and its convenience constructor [`RegularizedCS`](@ref)).
"""
abstract type CrossSection end

"""
    PowerLawCS(; σ0, w, n)

Power-law cross section per unit mass, `σ(v) = σ0 * (v/w)^n`. This is what
[`RunParams`](@ref) builds by default from its `σ0`, `w`, `n` fields.
"""
@kwdef struct PowerLawCS <: CrossSection
    σ0::Float64
    w::Float64
    n::Float64
end

# Uniform grid in ln v² with O(1) lookup. Linear interpolation between nodes and
# linear extrapolation beyond them, which is exact for power-law tails.
struct LogTable
    lv2_min::Float64
    inv_dl::Float64
    vals::Vector{Float64}
end

function lookup(tbl::LogTable, lv2)
    N = length(tbl.vals)
    t = (lv2 - tbl.lv2_min) * tbl.inv_dl
    i = clamp(floor(Int, t) + 1, 1, N - 1)
    frac = t - (i - 1)
    @inbounds return tbl.vals[i] + frac * (tbl.vals[i+1] - tbl.vals[i])
end

"""
    InterpolatedCS(σ̂, σ0; v2_center, halfwidth_lv2=45.0, npoints=4096)

Cross section of arbitrary shape, `σ(v) = σ0 * σ̂(v^2)`. Note that `σ̂` takes
the *squared* velocity `v^2` and returns the dimensionless shape `σ/σ0`.

The transport moments are tabulated once at construction on `npoints` points
uniform in `ln v^2`, spanning `ln(v2_center) ± halfwidth_lv2`; outside that range
they are extrapolated linearly in log-log. Construction takes a few seconds.

Pass the result to [`RunParams`](@ref) as `cs=...`, together with `n = NaN`.
"""
struct InterpolatedCS{F} <: CrossSection
    σ̂::F                  # normalized shape σ̂(v²) = σ(v)/σ0, called with v² directly
    σ0::Float64
    logK5::LogTable        # ln K5        vs ln v²
    logK7::LogTable        # ln K7        vs ln v²
    logKeff2::LogTable     # ln K_eff^(2) vs ln v²
end

# Maxwellian moments (Eq. K_p):
#   K_p(v2) = c_p ∫₀^∞ y^(p+2) e^(−y²) σ̂(4 v2 y²) dy,   (c5,c7,c9) = (1/3,1/12,1/60)
function _moment(σ̂, v2, p_pow, c_p)
    integral, _ = quadgk(y -> y^p_pow * exp(-y^2) * σ̂(4 * v2 * y^2), 0, 8; rtol=1e-10)
    return c_p * integral
end

# Tables span ln(v2_center) ± halfwidth_lv2 (default ±45, i.e. ±19.5 decades).
function InterpolatedCS(σ̂::F, σ0::Real; v2_center::Real, halfwidth_lv2::Real=45.0,
                         npoints::Integer=4096) where F
    lv2_min = log(v2_center) - halfwidth_lv2
    lv2_max = log(v2_center) + halfwidth_lv2
    dl = (lv2_max - lv2_min) / (npoints - 1)
    logK5_vals    = Vector{Float64}(undef, npoints)
    logK7_vals    = Vector{Float64}(undef, npoints)
    logKeff2_vals = Vector{Float64}(undef, npoints)
    for i in 1:npoints
        lv2 = lv2_min + (i - 1) * dl
        v2  = exp(lv2)
        k5 = _moment(σ̂, v2, 7,  1/3)
        k7 = _moment(σ̂, v2, 9,  1/12)
        k9 = _moment(σ̂, v2, 11, 1/60)
        # Eq. k_secondorder; a non-positive denominator means the second-Sonine
        # approximation has broken down for this σ̂.
        denom = 77*k5 - 112*k7 + 80*k9
        @assert denom > 0 "Keff2 denominator non-positive at lv2=$lv2 (second-Sonine breakdown)"
        keff2 = (28*k5^2 + 80*k5*k9 - 64*k7^2) / denom
        logK5_vals[i]    = log(k5)
        logK7_vals[i]    = log(k7)
        logKeff2_vals[i] = log(keff2)
    end
    inv_dl = 1 / dl
    logK5    = LogTable(lv2_min, inv_dl, logK5_vals)
    logK7    = LogTable(lv2_min, inv_dl, logK7_vals)
    logKeff2 = LogTable(lv2_min, inv_dl, logKeff2_vals)
    return InterpolatedCS{F}(σ̂, Float64(σ0), logK5, logK7, logKeff2)
end

# Power-law tails (n = 0 as v² → 0, n = −4 as v² → ∞), so the table
# extrapolation is asymptotically exact.
"""
    RegularizedCS(σ0, w)

[`InterpolatedCS`](@ref) for `σ(v) = σ0 / (1 + v^2/w^2)^2`.
"""
RegularizedCS(σ0::Real, w::Real) =
    InterpolatedCS(v2 -> 1 / (1 + v2/w^2)^2, σ0; v2_center = w^2)

# Direct σ(v), used by the evaporation integral.
σ_at_v2(cs::PowerLawCS, v2) = cs.σ0 * (sqrt(v2) / cs.w)^cs.n
σ_at_v2(cs::InterpolatedCS, v2) = cs.σ0 * cs.σ̂(v2)

# Log-space constants for the InterpolatedCS conductivity (calc_κ_log):
#   ln Kinv_smfp = lb_s + ln K_eff^(2)(v²) − (1/2) ln v²        b = 25√π/32
#   ln Kinv_lmfp = lb_l − ln ρ − (3/2) ln v² − ln K_lmfp(v²)    a = 4/√π
# PowerLawCS uses lci_s/lci_l instead, so these are NaN for it.
smfp_logconst(cs::PowerLawCS) = NaN
lmfp_logconst(cs::PowerLawCS, G, C) = NaN
smfp_logconst(cs::InterpolatedCS) = log(cs.σ0 * 32 / (25 * sqrt(pi)))
lmfp_logconst(cs::InterpolatedCS, G, C) = log(pi^1.5 * G / (C * cs.σ0))

"""
    RunParams(; G, σ0, w, n, t_end, r_max, r_min, gridpoints,
                snap_interval_time, snap_interval_density, snap_loc, kwargs...)

Parameters for a run of [`evolve`](@ref). Units are whatever system `G` is given in.

Required keywords:
- `G`: gravitational constant.
- `σ0`, `w`, `n`: cross section per unit mass `σ(v) = σ0 (v/w)^n`.
- `t_end`: stop time.
- `r_max`, `r_min`: outer boundary radius, and the inner radius used to set up the grid.
- `gridpoints`: number of mass-grid faces (snapshots have `gridpoints - 1` cells).
- `snap_interval_time`: write a snapshot at least this often in time (`Inf` disables).
- `snap_interval_density`: also write one when the central density has changed
  by this fractional amount since the last snapshot (`Inf` disables).
- `snap_loc`: output HDF5 file path (overwritten).

Optional keywords:
- `convection = ConvectionOn()`: or `ConvectionOff()`.
- `evaporation = EvaporationOn()`: or `EvaporationOff()`.
- `cs = PowerLawCS(σ0=σ0, w=w, n=n)`: any [`CrossSection`](@ref). With a
  non-power-law `cs`, set `n = NaN`.
- `C = 0.8505`, `bump_A = 1.31`, `bump_W = 1.68`: conductivity calibration constants.
- `collapse_density_factor = Inf`, `collapse_mass_frac = 0.0`: stop early once
  the density at enclosed mass `collapse_mass_frac * M_tot` has grown by
  `collapse_density_factor` over its initial value.
- `regrid = false`: adaptively redistribute mass shells during the run.
- `lmfp_K7_frac = 0.0`: K5/K7 mixture in the long-mean-free-path conductivity.
  The LMFP velocity moment K5 is replaced by `K_lmfp = (1 − a) K5 + a K7` with
  `a = lmfp_K7_frac`. `0` is the standard K5 closure. Constant σ is unaffected
  (K7 = K5); for a power law `K7/K5 = 1 + n/8`, so `a > 0` lowers the LMFP
  conductivity of a falling σ(v). `a = 0.73` fits deep-LMFP KiSS-SIDM runs
  with σ ∝ v⁻⁴ at large v.

The remaining fields are derived from the above; do not set them.
"""
@kwdef struct RunParams{TConv <: ConvectionMode, TEvap <: EvaporationMode, TCS <: CrossSection}

    convection::TConv = ConvectionOn()
    evaporation::TEvap = EvaporationOn()
    G::Float64
    # LMFP conductivity calibration constant (Eq. kappa_lmfp), the deep-LMFP
    # anchor that goes with the (A, W) bump below.
    C::Float64 = 0.8505
    σ0::Float64
    w::Float64
    n::Float64
    # For a non-power-law shape, pass cs explicitly and set n = NaN; the
    # power-law fields below are then unused.
    cs::TCS = PowerLawCS(σ0=σ0, w=w, n=n)
    lb_s::Float64 = smfp_logconst(cs)
    lb_l::Float64 = lmfp_logconst(cs, G, C)
    t_end::Float64
    r_max::Float64
    r_min::Float64
    gridpoints::Int64
    # Redistribute mass shells toward steep gradients each outer step, with a
    # v²-weighted donor-cell remap that conserves ∫v²s dm locally.
    regrid::Bool = false
    snap_interval_time::Float64
    snap_interval_density::Float64
    snap_loc::String
    collapse_density_factor::Float64 = Inf
    collapse_mass_frac::Float64 = 0.0
    # SMFP–LMFP crossover: a log-Gaussian bump on the summed resistivities,
    #   1/κ = (Kinv_smfp + Kinv_lmfp) · [1 + A exp(−u²/2W²)],  u = ln(Kinv_smfp/Kinv_lmfp).
    # The bump decays faster than any power, so both limits are the pure branches.
    # (A, W) are measured from the heat-flux moment of KiSS-SIDM runs that span
    # the crossover.
    bump_A::Float64 = 1.31
    bump_W::Float64 = 1.68
    # LMFP velocity weighting K_lmfp = (1 − a) K5 + a K7, a = lmfp_K7_frac.
    # K7 weights pair velocities as y⁹e^{−y²} rather than y⁷e^{−y²}, so a > 0
    # lowers κ_lmfp for a falling σ(v) (K7/K5 = 1 + n/8 for a power law) and leaves
    # constant σ, hence C, unchanged. Not gated by collisionality: it acts wherever
    # the LMFP branch matters.
    lmfp_K7_frac::Float64 = 0.0
    # Power law σ = σ0 (v/w)^n (Eq. K_p): K_p = K_p_prefactor · v2^(n/2), K_p ∝ Γ((n+p+3)/2).
    K5_prefactor::Float64 = ((1/3)  * 2^(n-1) * gamma(n/2 + 4)) / (w^n)
    K7_prefactor::Float64 = ((1/12) * 2^(n-1) * gamma(n/2 + 5)) / (w^n)
    K9_prefactor::Float64 = ((1/60) * 2^(n-1) * gamma(n/2 + 6)) / (w^n)
    # All K_p share the v2^(n/2) scaling, so both branch resistivities are power laws:
    #   Kinv_smfp = exp(lci_s + (n−1)/2 · ln v²)            b = 25√π/32
    #   Kinv_lmfp = exp(lci_l − ln ρ − (n+3)/2 · ln v²)     a = 4/√π
    Keff2_prefactor::Float64 = (28*K5_prefactor^2 + 80*K5_prefactor*K9_prefactor -
                                64*K7_prefactor^2) /
                               (77*K5_prefactor - 112*K7_prefactor + 80*K9_prefactor)
    lci_s::Float64 = log(σ0 * Keff2_prefactor * 32 / (25 * sqrt(pi)))
    Klmfp_prefactor::Float64 = lmfp_K7_frac == 0 ? K5_prefactor :
                               (1 - lmfp_K7_frac) * K5_prefactor + lmfp_K7_frac * K7_prefactor
    lci_l::Float64 = log(pi^1.5 * G / (C * σ0 * Klmfp_prefactor))
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