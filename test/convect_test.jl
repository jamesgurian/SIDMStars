using PyPlot
using DataInterpolations
using Printf
using SIDMStars

# convect! dispatches on run_params.convection (ConvectionOn/Off); only the
# convection mode matters here, the rest of RunParams is unused by convect!.
run_params = RunParams(
    convection = ConvectionOn(),
    G = 4.30091727003628e-6,
    σ0 = 1e-6,
    w = 10.,
    n = 0.,
    t_end = 1.,
    r_max = 1.,
    r_min = 1e-3,
    gridpoints = 10,
    snap_interval_time = Inf,
    snap_interval_density = Inf,
    snap_loc = "",
)

N = 200
m = range(0.1, 1.0, length=N) # Mass coordinate
s_rad = collect(m)            # Base stable profile (linear rise)

# B. Create Instability
# Add a "Dip" (entropy inversion). In this context, s dropping is unstable.
# We create a Gaussian dip centered at m=0.5
instability = -0.15 * exp.(-(m .- 0.5).^2 ./ (2*0.05^2))
s_initial = s_rad + instability

s_interp = AkimaInterpolation(s_initial, log.(m))
# C. Create Non-Isothermal Weighting (Temperature)
# T drops as 1/sqrt(m) (hotter at center)
v2 = 1.0 ./ sqrt.(m)

# D. Run Simulation
s_final = copy(s_initial)
SIDMStars.convect!(s_interp, s_final, v2, m, run_params)

# E. Verification (Conservation of Energy)
# We compute \int (s_final - s_initial) * T * dm over the whole domain
# This should be practically zero.
energy_diff = sum((s_final .- s_initial) .* v2 .* (vcat(0.0, diff(m))))
@printf "Net Weighted Entropy Change: %.2e (Should be close to 0)\n" energy_diff
semilogx(m, s_initial)
semilogx(m, s_final,"--")


##Now the case where the zone extends to the core
instability_core = -0.2 * exp.(-(m .- 0.2).^2 ./ (2*0.1^2)) 
s_initial_core = s_rad + instability_core

s_interp_core = AkimaInterpolation(s_initial_core, log.(m))
s_final_core = copy(s_initial_core)
semilogx(m, s_initial_core)

SIDMStars.convect!(s_interp_core, s_final_core, v2, m, run_params)

# Verify
energy_diff_core = sum((s_final_core .- s_initial_core) .* v2 .* (vcat(0.0, diff(m))))
@printf "Core Net Change: %.2e\n" energy_diff_core

semilogx(m, s_final_core, "--")


##Now the case where the zone extends to the surface
instability_surf = -0.2 * exp.(-(m .- 1.1).^2 ./ (2*0.1^2)) 
s_initial_surf = s_rad + instability_surf

s_interp_surf = AkimaInterpolation(s_initial_surf, log.(m))
s_final_surf = copy(s_initial_surf)
semilogx(m, s_initial_surf)

SIDMStars.convect!(s_interp_surf, s_final_surf, v2, m, run_params)

# Verify
energy_diff_surf = sum((s_final_surf .- s_initial_surf) .* v2 .* (vcat(0.0, diff(m))))
@printf "Surface Net Change: %.2e\n" energy_diff_surf

semilogx(m, s_final_surf, "--")