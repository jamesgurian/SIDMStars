
using Unitful
using UnitfulAstro
include("../src/SIDMStars.jl")
using Revise
using .SIDMStars
using PyPlot


    


# 3. Setup Params
dm_params = RunParams(
    C = 1.0,
    σ0 = ustrip(uconvert(u"kpc^2/Msun", 100.0u"cm^2/g")),
    w = ustrip(10.0u"km/s"), # Example velocity scale
    n = 1.0,
    t_end = 10.0,
    r_max = ustrip(100.0u"kpc"),
    r_min = ustrip(1.0e-3u"kpc")
)
sol_init = initialize_profile(dm_params, x-> nfw_profile(1, 1e8, x), p_outer=0.)
#semilogx(sol[1].t, sol[2].(sol[1].t))
savefig("out.png")
sol = evolve(dm_params, x -> nfw_profile(1, 1.0e8, x));
