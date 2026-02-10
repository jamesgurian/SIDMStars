
include("../src/SIDMStars.jl")

using .SIDMStars
using PyPlot

#using PyPlot


# 3. Setup Params
dm_params = RunParams(
    C = 1.0,
    σ0 = 1.04e-6, #kpc^2/M_sun
    w = 10, # km/s
    n = -4.,
    t_end = 7,
    r_max = 2, #kpc
    r_min = 1e-3, #kpc
    snap_interval_time = 0.1,
    snap_interval_density = 0.1,
    snap_loc="init_out.hdf5"
)

profile = setup_dk_profile(1e6, 7, 1, z=10)
sol_init = initialize_profile(dm_params,profile, p_outer=0);
sol = evolve(dm_params,profile,gridpoints=100,eps=0.1, mass_min_frac=1e-4);
