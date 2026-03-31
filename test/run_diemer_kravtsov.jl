
using Revise
using SIDMStars
using PyPlot

#using PyPlot

# 3. Setup Params
dm_params = RunParams(
    C = 1.0,
    G = 4.30091727003628e-6,
    σ0 = 1.04e-6, #kpc^2/M_sun
    w = 10., # km/s
    n = -4.,
    t_end = 10.,
    r_max = 2., #kpc
    r_min = 1e-3, #kpc
    snap_interval_time = 0.1,
    snap_interval_density = 0.1,
    gridpoints=150,
    snap_loc="dk_out.hdf5"
)

profile = setup_dk_profile(5e6, 7, 1.0, z=10)
sol_init = initialize_profile(dm_params,profile, p_outer=0);
#eps_diff: inverse number of timesteps per von Neumann constraint on dt
#eps_kh: inverse number of timesteps per 1/\dot s
evolve(dm_params,profile,eps_diff=0.2, eps_kh=.05, mass_min_frac=1e-4);
snap = read_snapshot_file("dk_out.hdf5")
#loglog(snap["r"], snap["rho"])