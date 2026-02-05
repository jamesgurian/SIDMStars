
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
    t_end = 4.5,
    r_max = 2, #kpc
    r_min = 1e-3, #kpc
    snaptimes = range(.01, stop = 4.5, length=30),
    snap_loc="init_out.hdf5"
)

profile = setup_dk_profile(1e6, 7, 1, z=10)
sol_init = initialize_profile(dm_params,profile, p_outer=0);
#semilogx(sol[1].t, sol[2].(sol[1].t))
#savefig("out.png")
sol = evolve(dm_params,profile,gridpoints=100,eps=0.1, mass_min_frac=1e-4);

#@profview sol = evolve(dm_params,profile,gridpoints=200,eps=0.1, mass_min_frac=1e-4);

#@profview sol = evolve(dm_params, x -> nfw_profile(1, 1.0e7,  x),gridpoints=50,eps=0.001, mass_min_frac=1e-4);
