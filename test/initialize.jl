
include("../src/SIDMStars.jl")

using .SIDMStars
using PyPlot

    


# 3. Setup Params
dm_params = RunParams(
    C = 1.0,
    σ0 = 5e-5, #kpc^2/M_sun
    w = 10, # km/s
    n = -4.,
    t_end = 1,
    r_max = 10, #kpc
    r_min = 1e-3, #kpc
    snaptimes = range(.01, stop = 1, length=30)
)

profile = setup_dk_profile(1e6, 7, 1, z=10)
sol_init = initialize_profile(dm_params,profile, p_outer=0)
#semilogx(sol[1].t, sol[2].(sol[1].t))
#savefig("out.png")
sol = evolve(dm_params,profile,gridpoints=200,eps=0.001, mass_min_frac=1e-4);
#@profview sol = evolve(dm_params, x -> nfw_profile(1, 1.0e7,  x),gridpoints=50,eps=0.001, mass_min_frac=1e-4);
