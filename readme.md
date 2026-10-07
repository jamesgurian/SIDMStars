# SIDMStars

A 1D Lagrangian gravothermal fluid code for self-interacting dark matter (SIDM) halos,
written in Julia. It evolves a spherical halo in hydrostatic equilibrium under heat
conduction by dark matter self-interactions, with optional convection and evaporation,
and writes profiles to HDF5.

## Installation

From the Julia REPL, press `]` to enter the package manager, then:

```
pkg> dev /path/to/SIDMStars
```

## Quick start

```julia
using SIDMStars

G     = 4.30091727003628e-6          # kpc (km/s)^2 / M_sun
r_max = 30.0                         # kpc
ρ(r)  = nfw_profile(1.18, 2.73e7, r) # r_s [kpc], rho_s [M_sun/kpc^3]

# Pressure at r_max from the weight of the halo outside it.
p_outer = calc_psurf(ρ, r_max, G)

params = RunParams(
    G          = G,
    σ0         = 2.09e-10 * 50,      # 50 cm^2/g in kpc^2/M_sun
    w          = 10.0,               # km/s
    n          = 0.0,                # velocity-independent
    t_end      = 25.0,               # kpc/(km/s) ≈ 0.978 Gyr
    r_max      = r_max,
    r_min      = 1e-3,
    gridpoints = 200,
    snap_interval_time    = 0.1,
    snap_interval_density = 0.1,
    snap_loc   = "nfw_out.hdf5",
    convection  = ConvectionOff(),
    evaporation = EvaporationOff(),
)

evolve(params, ρ; p_outer = p_outer, mass_min_frac = 1e-6, eps_diff = 0.2, eps_kh = 0.05)

snap = read_snapshot_file("nfw_out.hdf5")
snap["t"]            # snapshot times
snap["rho"][1, :]    # central density vs. time
```

More examples are in `test/` (`run_nfw.jl`, `run_diemer_kravtsov.jl`, `run_isothermal.jl`).
Each exported name has a docstring, so `?evolve` or `?RunParams` in the REPL shows the
details.

## Units

The code has no built-in units. You pick them by the value of `G`, and every other input
must be in the same system. The examples use kpc, M_sun and km/s:

| Quantity      | Unit                     |
|---------------|--------------------------|
| `G`           | 4.30091727e-6 kpc (km/s)² / M_sun |
| length        | kpc                      |
| mass          | M_sun                    |
| velocity, `w` | km/s                     |
| time          | kpc/(km/s) ≈ 0.978 Gyr   |
| `σ0`          | kpc²/M_sun (1 cm²/g ≈ 2.09e-10 kpc²/M_sun) |

`setup_dk_profile` is the exception: its inputs are hardwired to kpc and M_sun.

## API

### Running

| Function | Description |
|----------|-------------|
| `evolve(run_params, rho_initial; kwargs...)` | Build the initial state from `rho_initial(r)`, evolve it, and write snapshots to `run_params.snap_loc`. Returns `nothing`. |
| `initialize_profile(run_params, rho_initial; p_outer, mass_min_frac)` | Only the initial hydrostatic-equilibrium solve. Returns `(sol, s_func, M_tot)`. `evolve` calls it for you. |
| `calc_psurf(rho_func, r_max, G)` | Pressure at `r_max` from the weight of the material outside it (hydrostatic equilibrium). Pass it as `p_outer`. |

Keyword arguments to `evolve`:

| Keyword | Default | Meaning |
|---------|---------|---------|
| `p_outer` | `0` | Pressure at `r_max`. Set it to the weight of the material outside `r_max`, e.g. `calc_psurf(rho_initial, r_max, G)`; `0` can destabilize the outer shells. |
| `mass_min_frac` | `1e-6` | Inner mass boundary, as a fraction of the mass inside `r_max`. Inside it the density is held flat. |
| `eps_diff` | `1` | Timestep as a fraction of the explicit diffusion stability limit. |
| `eps_kh` | `0.05` | Timestep as a fraction of the shortest entropy-change timescale. |
| `velocity_bins` | `50` | Velocity resolution of the evaporation calculation. |
| `abstol`, `reltol` | `1e-5` | Starting tolerances for the structure solves. They are tightened automatically if a solve fails. |

### `RunParams`

All run settings, passed as keyword arguments.

**Required**

| Field | Meaning |
|-------|---------|
| `G` | Gravitational constant (fixes the unit system). |
| `σ0`, `w`, `n` | Cross section per unit mass, `σ(v) = σ0 (v/w)^n`. |
| `t_end` | Stop time. |
| `r_max` | Outer boundary radius. |
| `r_min` | Inner radius used to set up the grid. |
| `gridpoints` | Number of grid faces. Snapshots have `gridpoints - 1` cells. |
| `snap_interval_time` | Maximum time between snapshots (`Inf` disables). |
| `snap_interval_density` | Also write a snapshot when the central density has changed by this fraction (`Inf` disables). |
| `snap_loc` | Output HDF5 path. An existing file there is overwritten. |

**Optional**

| Field | Default | Meaning |
|-------|---------|---------|
| `convection` | `ConvectionOn()` | Or `ConvectionOff()`. |
| `evaporation` | `EvaporationOn()` | Or `EvaporationOff()`. |
| `cs` | `PowerLawCS(σ0, w, n)` | Cross-section object (see below). |
| `C` | `0.8620` | Conductivity calibration constant. |
| `bump_A`, `bump_W` | `1.31`, `1.68` | Conductivity calibration constants. |
| `collapse_density_factor` | `Inf` | Stop early once the density at enclosed mass `collapse_mass_frac * M_tot` reaches this multiple of its initial value. |
| `collapse_mass_frac` | `0.0` | See above. `0` tracks the innermost cell. |
| `regrid` | `false` | Redistribute mass shells toward steep gradients during the run. |
| `lmfp_K7_frac` | `0.0` | Γ₅/Γ₇ mixture in the LMFP conductivity (see below). `0` is the standard Γ₅ closure. |

The remaining fields of `RunParams` are derived from these; do not set them.

### Conductivity and the Γ₅/Γ₇ mixture

The conductivity joins the short- and long-mean-free-path (SMFP, LMFP) limits through a
log-Gaussian bump in the summed resistivities,

```
1/κ = (1/κ_SMFP + 1/κ_LMFP) · [1 + A exp(−u²/2W²)],   u = ln(κ_LMFP/κ_SMFP),
```

with `A = bump_A` and `W = bump_W`. The LMFP branch has the standard calibration constant `C`.

The cross section enters through Maxwellian velocity moments
`Γ_p = σ0 K_p`, with `K_p ∝ ∫ y^(p+2) e^(−y²) σ̂ dy` (`K5`, `K7`, `K9` in the code). The SMFP
conductivity uses the second-order Chapman–Enskog combination of K5, K7 and K9. The LMFP conductivity
uses Γ₅ by default, but this may underestimate the collapse time for velocity-dependent cross sections in the deep LMFP. 

`lmfp_K7_frac = a` replaces Γ₅ in the LMFP branch by the mixture

```
Γ_LMFP = (1 − a) Γ₅ + a Γ₇.
```

- Γ₇ weights the pair velocities more heavily toward fast pairs. For a falling σ(v) it is
  smaller than Γ₅, so `a > 0` lowers the LMFP conductivity. For a power law,
  `Γ₇/Γ₅ = 1 + n/8`.
- For constant σ, Γ₇ = Γ₅, so the switch has no effect and `C` keeps its calibration.
- The mixture applies wherever the LMFP branch contributes. It is not gated by the local
  Knudsen number.
- `a = 0` (the default) is the standard closure. `a = 0.73` fits deep-LMFP KiSS-SIDM runs of
  cross sections that fall as v⁻⁴ at high velocity. This also
  slows collapse near the SMFP–LMFP crossover. In those v⁻⁴ tests, crossover collapse times came
  out 20–30% longer than KiSS-SIDM. For that reason it is off by default.
- It works with every cross-section type. `InterpolatedCS` tabulates K7 alongside K5.

```julia
params = RunParams(; G, σ0, w, n = -4.0, lmfp_K7_frac = 0.73, ...)
```

### Cross sections

`RunParams` builds a power law from `σ0`, `w`, `n` by default. For any other velocity
dependence, pass a cross-section object as `cs` **and set `n = NaN`**. `σ0` and `w` are
still required fields, so pass the same values used to build `cs`.

| Type | Cross section |
|------|---------------|
| `PowerLawCS(; σ0, w, n)` | `σ0 (v/w)^n` (the default) |
| `RegularizedCS(σ0, w)` | `σ0 / (1 + v²/w²)²` |
| `InterpolatedCS(σ̂, σ0; v2_center, halfwidth_lv2=45, npoints=4096)` | `σ0 σ̂(v²)` for any shape `σ̂`. Note that **`σ̂` takes `v²`**, not `v`. |

`InterpolatedCS` tabulates what the solver needs on a log grid centered on `v2_center`
when it is constructed, which takes a few seconds, and extrapolates beyond the table as
a power law.

```julia
cs = RegularizedCS(1e-6, 10.0)
params = RunParams(; G, σ0 = 1e-6, w = 10.0, n = NaN, cs = cs, ...)
```

### Initial profiles

| Function | Description |
|----------|-------------|
| `nfw_profile(r_s, rho_0, r)` | NFW density at `r`. |
| `dk_profile(r_s, rho_0, rvir, nu, r)` | NFW density times the truncation factor `(1 + (r/rvir)^4/5)^(-nu)`. |
| `setup_dk_profile(M_vir, c_vir, nu; z=0, delta_vir=200, H0=67.74, Om0=0.3089)` | Returns `r -> ρ(r)` for a truncated NFW halo of virial mass `M_vir` [M_sun] and concentration `c_vir`. The virial overdensity is relative to the mean matter density. |

`rho_initial` can be any function of radius; these are just conveniences.

### Output

`read_snapshot_file(filename)` returns a `Dict`:

| Key | Contents |
|-----|----------|
| `"t"` | Snapshot times, as a vector. |
| `"r"`, `"m"`, `"P"`, `"rho"`, `"s"`, `"L"` | Matrices of size `(gridpoints - 1, n_snapshots)` at cell centers: radius, enclosed mass, pressure, density, entropy `log(P^(3/2)/ρ^(5/2))`, and luminosity. Index 1 is the innermost cell. |
| `"metadata"` | `Dict` of the run's `RunParams` fields, plus `git_commit`. |

Snapshots are appended to the file as the run proceeds.

### Other exports

- `Constants`: the CODATA 2018 physical constants, with units, from PhysicalConstants.jl.
  For example, `Constants.G` and `Constants.k_B`.
