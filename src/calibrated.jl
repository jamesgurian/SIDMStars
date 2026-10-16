# Simplified interface: RunParams and evolve() with the calibrated conductivity
# and the numerical settings used for the production runs filled in.

# C, bump_A and bump_W are the RunParams defaults. The Γ₅/Γ₇ mixture fitted to
# deep-LMFP KiSS-SIDM runs with σ ∝ v⁻⁴ tails is used when calib = :lmfp (see
# readme, "Conductivity and the Γ₅/Γ₇ mixture").
const CALIBRATED_LMFP_K7_FRAC = 0.73

# G in kpc (km/s)² / M_sun.
const G_KPC_KMS_MSUN = 4.30091727003628e-6

as_cross_section(cs::CrossSection) = cs
as_cross_section(σ0::Real) = PowerLawCS(σ0 = σ0, w = 1.0, n = 0.0)

# Γ₅/Γ₇ mixture for each calibration choice. It only matters for velocity-dependent
# σ: for constant σ, Γ₇ = Γ₅ and both choices give the same conductivity.
function calib_lmfp_K7_frac(calib::Symbol)
    calib === :smfp && return 0.0
    calib === :lmfp && return CALIBRATED_LMFP_K7_FRAC
    throw(ArgumentError("calib must be :smfp or :lmfp, got :$calib"))
end

# The (σ0, w, n) fields RunParams requires alongside cs.
power_law_fields(cs::PowerLawCS) = (σ0 = cs.σ0, w = cs.w, n = cs.n)
power_law_fields(cs::InterpolatedCS) = (σ0 = cs.σ0, w = NaN, n = NaN)

"""
    calibrated_params(cs; calib=:smfp, t_end, r_max, snap_loc, kwargs...)

[`RunParams`](@ref) with the calibrated conductivity. `cs` is any
[`CrossSection`](@ref), or a number for a constant cross section `σ0`. `σ0`,
`w` and `n` are taken from `cs`.

Calibration:
- `C`, `bump_A`, `bump_W`: the [`RunParams`](@ref) defaults (0.8505, 1.31,
  1.68), for every cross section.
- `calib` sets the LMFP velocity weighting for a velocity-dependent σ:
  - `:smfp` (default): standard Γ₅ closure (`lmfp_K7_frac = 0`). Use it for
    halos in the SMFP regime or near the SMFP–LMFP crossover.
  - `:lmfp`: Γ₅/Γ₇ mixture (`lmfp_K7_frac = 0.73`), fitted to deep-LMFP
    KiSS-SIDM runs. Near the crossover it gives collapse times 20–30% too long.
  For constant σ, Γ₇ = Γ₅, so `calib` has no effect. An explicit
  `lmfp_K7_frac` overrides `calib`.

Other defaults: `G` in kpc (km/s)²/M_sun, `r_min = 1e-6 r_max`,
`gridpoints = 200`, `snap_interval_time = t_end/1000`,
`snap_interval_density = 0.1`, `convection = ConvectionOn()`,
`evaporation = EvaporationOff()`.

Every default can be overridden by keyword. Any other `RunParams` keyword
(e.g. `collapse_density_factor`, `regrid`) is passed through.

```julia
params = calibrated_params(RegularizedCS(1e-6, 10.0); calib = :lmfp,
                           t_end = 6000.0, r_max = 40.0, snap_loc = "out.hdf5")
```
"""
function calibrated_params(cs; t_end, r_max, snap_loc,
                           G = G_KPC_KMS_MSUN,
                           r_min = 1e-6 * r_max,
                           gridpoints = 200,
                           snap_interval_time = t_end / 1000,
                           snap_interval_density = 0.1,
                           calib::Symbol = :smfp,
                           lmfp_K7_frac = nothing,
                           convection = ConvectionOn(),
                           evaporation = EvaporationOff(),
                           kwargs...)
    cs = as_cross_section(cs)
    a7 = something(lmfp_K7_frac, calib_lmfp_K7_frac(calib))
    return RunParams(; power_law_fields(cs)..., cs, G,
                     lmfp_K7_frac = a7, t_end, r_max, r_min, gridpoints,
                     snap_interval_time, snap_interval_density, snap_loc,
                     convection, evaporation, kwargs...)
end

"""
    run_halo(rho_initial, cs; calib=:smfp, t_end, r_max, snap_loc, kwargs...)

Evolve the halo with density profile `rho_initial(r)` and cross section `cs`
(a [`CrossSection`](@ref), or a number for constant σ) using
[`calibrated_params`](@ref), and write snapshots to `snap_loc`. Returns the
`RunParams` used.

`p_outer` defaults to [`calc_psurf`](@ref)`(rho_initial, r_max, G)`. The
[`evolve`](@ref) numerics default to the production settings `eps_kh = 0.01`,
`eps_diff = 0.5`, `mass_min_frac = 1e-7`. All other keywords, including
`calib = :smfp | :lmfp`, go to [`calibrated_params`](@ref).

```julia
ρ(r) = nfw_profile(1.18, 2.73e7, r)
run_halo(ρ, 50 * 2.09e-10; t_end = 6000.0, r_max = 40.0, snap_loc = "iso.hdf5")
run_halo(ρ, RegularizedCS(1e-6, 10.0); calib = :lmfp, t_end = 6000.0,
         r_max = 40.0, snap_loc = "vdep_lmfp.hdf5")
```
"""
function run_halo(rho_initial, cs; r_max, G = G_KPC_KMS_MSUN,
                  p_outer = calc_psurf(rho_initial, r_max, G),
                  eps_kh = 0.01, eps_diff = 0.5, mass_min_frac = 1e-7,
                  velocity_bins = 50, abstol = 1e-5, reltol = 1e-5,
                  kwargs...)
    params = calibrated_params(cs; r_max, G, kwargs...)
    evolve(params, rho_initial; p_outer, eps_kh, eps_diff, mass_min_frac,
           velocity_bins, abstol, reltol)
    return params
end
