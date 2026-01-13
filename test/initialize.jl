using FunctionWrappers
import FunctionWrappers:FunctionWrapper
using Unitful
using UnitfulAstro
include("../src/SIDMStars.jl")
using Revise
using .SIDMStars
rho_0 = 1.0e7#u"Msun/kpc^3"
R_s = 1.0#u"kpc"
  

    


# 3. Setup Params
dm_params = RunParams(
    C = 1.0,
    σ0 = ustrip(uconvert(u"kpc^2/Msun", 100.0u"cm^2/g")),
    w = ustrip(10.0u"km/s"), # Example velocity scale
    n = 1.0,
    t_end = 10.0,
    r_max = ustrip(100.0u"kpc"),
    r_min = ustrip(1.0e-3u"kpc"),
    r_s=R_s,
    rho_s=1e7
)

sol = initialize_profile(dm_params);