
include("../src/SIDMStars.jl")

using .SIDMStars
using PyPlot


    


# 3. Setup Params
dm_params = RunParams(
    C = 1.0,
    σ0 = 1e-8, #kpc^2/msun
    w = 10, # km/s
    n = 0.0,
    t_end = 1.0,
    r_max = 40, #kpc
    r_min = 1e-3, #kpc
    snap_interval_time = 0.1,
    snap_interval_density = 0.1
)
sol_init = initialize_profile(dm_params, x-> nfw_profile(1, 1e7, x), p_outer=1e-20);
#semilogx(sol[1].t, sol[2].(sol[1].t))
#savefig("out.png")
sol = evolve(dm_params, x -> nfw_profile(1, 1.0e7,  x),gridpoints=50,eps=0.001, mass_min_frac=1e-4);