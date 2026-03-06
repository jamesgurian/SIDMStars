
using SIDMStars
using PyPlot


    


# 3. Setup Params
dm_params = RunParams(
    C = 1.0,
    σ0 = 1e-8, #kpc^2/msun
    w = 10, # km/s
    n = 0.0,
    t_end = 10.0,
    r_max = 30, #kpc
    r_min = 1e-3, #kpc
    snap_interval_time = 0.1,
    snap_interval_density = 0.1,
    snap_loc= "nfw_out.hdf5"
)
sol_init = initialize_profile(dm_params, x-> nfw_profile(1, 1e7, x), p_outer=1e5);
#semilogx(sol[1].t, sol[2].(sol[1].t))
#savefig("out.png")
sol = evolve(dm_params, x -> nfw_profile(1, 1.0e7,  x),gridpoints=150,eps_kh=.05, eps_diff=0.2, mass_min_frac=1e-4, p_outer=1e5);