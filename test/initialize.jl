
include("../src/SIDMStars.jl")

using .SIDMStars
using PyPlot

using ProfileView
    


# 3. Setup Params
dm_params = RunParams(
    C = 1.0,
    σ0 = 1e-5, #kpc^2/M_sun
    w = 10, # km/s
    n = -4.,
    t_end = 1,
    r_max = 100, #kpc
    r_min = 1e-3, #kpc
    snaptimes = range(.01, stop = 1, length=30)
)
sol_init = initialize_profile(dm_params, x-> nfw_profile(1, 1e7, x), p_outer=1e-20)
#semilogx(sol[1].t, sol[2].(sol[1].t))
#savefig("out.png")
sol = evolve(dm_params, x -> nfw_profile(1, 1.0e7,  x),gridpoints=100,eps=0.001, mass_min_frac=1e-4);
#@profview sol = evolve(dm_params, x -> nfw_profile(1, 1.0e7,  x),gridpoints=50,eps=0.001, mass_min_frac=1e-4);
