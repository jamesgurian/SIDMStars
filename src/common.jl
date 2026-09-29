

# Define abstract types for the categories
abstract type ConvectionMode end
abstract type EvaporationMode end

# Define singleton types for the specific options
struct ConvectionOn <: ConvectionMode end
struct ConvectionOff <: ConvectionMode end

struct EvaporationOn <: EvaporationMode end
struct EvaporationOff <: EvaporationMode end

# --- Cross-section type hierarchy -------------------------------------------
# Carries σ(v) into RunParams via multiple dispatch. PowerLawCS reproduces the
# original hardwired σ = σ0 (v/w)^n exactly (all existing calibrations —
# C=0.8620, bump (A,W)=(1.2046,1.7499), the Ψ model — were measured against
# it). InterpolatedCS supports an arbitrary normalized shape σ̂(v²)=σ(v)/σ0 by
# precomputing log-log tables of the K5 and K_eff^(2) moments (the only two
# the conductivity hot path needs) on a uniform grid in ln v², built once at
# construction via QuadGK.
abstract type CrossSection end

@kwdef struct PowerLawCS <: CrossSection
    σ0::Float64
    w::Float64
    n::Float64
end

# Uniform grid in ln v² with allocation-free O(1) lookup. Linear interpolation
# between nodes; outside the grid range the last segment's slope is used
# unmodified, i.e. plain linear extrapolation — exact for a shape whose tail
# is locally a power law (as for the regularized cross section below).
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

struct InterpolatedCS{F} <: CrossSection
    σ̂::F                  # normalized shape σ̂(v²) = σ(v)/σ0, called with v² directly
    σ0::Float64
    logK5::LogTable        # ln K5        vs ln v²
    logKeff2::LogTable     # ln K_eff^(2) vs ln v²
end

# Generic Maxwellian moments (Eq. K_p), verified analytically to reproduce the
# power-law prefactors above exactly (Γ-functions and all):
#   K_p(v2) = c_p ∫₀^∞ y^(p+2) e^(−y²) σ̂(4 v2 y²) dy,   (c5,c7,c9) = (1/3,1/12,1/60)
# (the integrand argument 4·v2·y² is v² for v = 2√v2·y). K_eff^(2) is the same
# rational combination of K5/K7/K9 as below (Eq. k_secondorder), valid for any σ.
function _moment(σ̂, v2, p_pow, c_p)
    integral, _ = quadgk(y -> y^p_pow * exp(-y^2) * σ̂(4 * v2 * y^2), 0, 8; rtol=1e-10)
    return c_p * integral
end

# Constructor: builds ln K5 and ln K_eff^(2) tables on a uniform grid in ln v²
# spanning `v2_center` ± `halfwidth_lv2` nats (default ±45 nats ≈ ±19.5 decades
# each side). Cost is ~3·npoints smooth quadratures — a few seconds, once per run.
function InterpolatedCS(σ̂::F, σ0::Real; v2_center::Real, halfwidth_lv2::Real=45.0,
                         npoints::Integer=4096) where F
    lv2_min = log(v2_center) - halfwidth_lv2
    lv2_max = log(v2_center) + halfwidth_lv2
    dl = (lv2_max - lv2_min) / (npoints - 1)
    logK5_vals    = Vector{Float64}(undef, npoints)
    logKeff2_vals = Vector{Float64}(undef, npoints)
    for i in 1:npoints
        lv2 = lv2_min + (i - 1) * dl
        v2  = exp(lv2)
        k5 = _moment(σ̂, v2, 7,  1/3)
        k7 = _moment(σ̂, v2, 9,  1/12)
        k9 = _moment(σ̂, v2, 11, 1/60)
        # Second-order Chapman-Enskog denominator (Eq. k_secondorder); a
        # non-positive value signals the second-Sonine approximation has
        # broken down for this σ̂ at this v2.
        denom = 77*k5 - 112*k7 + 80*k9
        @assert denom > 0 "Keff2 denominator non-positive at lv2=$lv2 (second-Sonine breakdown)"
        keff2 = (28*k5^2 + 80*k5*k9 - 64*k7^2) / denom
        logK5_vals[i]    = log(k5)
        logKeff2_vals[i] = log(keff2)
    end
    inv_dl = 1 / dl
    logK5    = LogTable(lv2_min, inv_dl, logK5_vals)
    logKeff2 = LogTable(lv2_min, inv_dl, logKeff2_vals)
    return InterpolatedCS{F}(σ̂, Float64(σ0), logK5, logKeff2)
end

# Convenience: the regularized shape used by the companion DSMC simulations,
#   σ̂(v²) = 1/(1+v²/w²)²,
# with exact power-law tails (n=0 as v²→0, n=−4 as v²→∞), so log-log linear
# extrapolation beyond the table is asymptotically exact.
RegularizedCS(σ0::Real, w::Real) =
    InterpolatedCS(v2 -> 1 / (1 + v2/w^2)^2, σ0; v2_center = w^2)

# Direct σ(v) evaluation — used only by the evaporation integral (transport.jl).
σ_at_v2(cs::PowerLawCS, v2) = cs.σ0 * (sqrt(v2) / cs.w)^cs.n
σ_at_v2(cs::InterpolatedCS, v2) = cs.σ0 * cs.σ̂(v2)

# Derived log-space constants for the InterpolatedCS calc_κ_log hot path
# (transport.jl) — mirror lci_s/lci_l below, but with K5/K_eff^(2) read from
# tables instead of closed-form power-law prefactors:
#   ln Kinv_smfp = lb_s + ln K_eff^(2)(v²) − (1/2) ln v²        b = 25√π/32
#   ln Kinv_lmfp = lb_l − ln ρ − (3/2) ln v² − ln K5(v²)        a = 4/√π
# PowerLawCS doesn't use these (the power-law hot path reads lci_s/lci_l
# directly); NaN here is unused, never read.
smfp_logconst(cs::PowerLawCS) = NaN
lmfp_logconst(cs::PowerLawCS, G, C) = NaN
smfp_logconst(cs::InterpolatedCS) = log(cs.σ0 * 32 / (25 * sqrt(pi)))
lmfp_logconst(cs::InterpolatedCS, G, C) = log(pi^1.5 * G / (C * cs.σ0))

@kwdef struct RunParams{TConv <: ConvectionMode, TEvap <: EvaporationMode, TCS <: CrossSection}

    convection::TConv = ConvectionOn()
    evaporation::TEvap = EvaporationOn()
    G::Float64
    # LMFP conductivity calibration constant, measured jointly with the
    # crossover bump below (fit_kappa_local.jl).
    C::Float64 = 0.8620
    σ0::Float64
    w::Float64
    n::Float64
    # Cross section carried via multiple dispatch (see CrossSection hierarchy
    # above). Defaults to the power law σ = σ0 (v/w)^n built from σ0/w/n
    # above, so existing keyword construction is unaffected. To use another
    # shape (e.g. RegularizedCS(σ0, w)), pass cs= explicitly and set n = NaN
    # (n then only propagates into the unused legacy prefactor/lci fields
    # below; the interpolated path never reads them).
    cs::TCS = PowerLawCS(σ0=σ0, w=w, n=n)
    # Derived log-space constants for the InterpolatedCS calc_κ_log hot path
    # (unused/NaN for PowerLawCS — see smfp_logconst/lmfp_logconst above).
    lb_s::Float64 = smfp_logconst(cs)
    lb_l::Float64 = lmfp_logconst(cs, G, C)
    t_end::Float64
    r_max::Float64
    r_min::Float64
    gridpoints::Int64
    # Adaptive-regrid toggle. When true, evolve calls update_grid! each outer
    # step to redistribute mass shells toward steep-gradient regions (e.g. the
    # migrating convective boundary), using a v²-weighted donor-cell remap that
    # conserves the heat integral ∫v²s dm (and ∫v² dm) locally. Default false so
    # the calibrated deep-SMFP/LMFP runs stay on the static grid, unchanged.
    regrid::Bool = false
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
    # Defaults are the DIRECT heat-flux-moment measurement (notes/flux_moment.tex,
    # flux_cross_joint.jl), the joint fit restricted to the three DSMC runs that
    # straddle the crossover peak (sigma = 7482, 25000, 83333 cm^2/g):
    #   (A, W) = (1.31 ± 0.05, 1.68 ± 0.09),   5167 bins, chi2/dof = 1.17.
    # Verified against fluid collapse times (rms ln(fluid/DSMC) = 0.044, vs 0.039
    # for the superseded set below and 0.085 for the all-six-run variant).
    # SUPERSEDES (A, W) = (1.2046 ± 0.039, 1.7499 ± 0.050) from the DSMC heating
    # trajectories (fit_kappa_local.jl), which were fit jointly with C = 0.8620.
    # C is deliberately NOT moved: the flux-moment reduction assumed the deep-LMFP
    # anchor C = 0.8505, but ln(0.8620/0.8505) = 0.013 shifts u by 0.8% of W, so
    # (A, W) are insensitive to the choice and every other calibration keeps 0.8620.
    bump_A::Float64 = 1.31
    bump_W::Float64 = 1.68
    # Diagnostic toggle: drop the LMFP branch (and the crossover bump, which
    # dies as Kinv_lmfp→0) so κ = 1/Kinv_smfp is the pure 2-Sonine SMFP
    # Chapman-Enskog conductivity everywhere. For testing how far the SMFP
    # closure alone tracks the collapse (see calc_κ_log in transport.jl).
    lmfp_off::Bool = false
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