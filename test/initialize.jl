using FunctionWrappers
import FunctionWrappers:FunctionWrapper
using Unitful
using UnitfulAstro
include("../src/SIDMStars.jl")
using .SIDMStars
rho_0 = 1.0e7u"Msun/kpc^3"
R_s = 1.0u"kpc"
  
function nfw_profile(r)
    x = r / R_s
    return rho_0 / (x * (1 + x)^2)
end
    
# Wrapper for the struct
rho_wrapper = FunctionWrapper{typeof(rho_0), Tuple{typeof(R_s)}}(nfw_profile)

# 3. Setup Params
dm_params = RunParams(
    C = 1.0,
    σ0 = 100.0u"km/s",
    w = 10.0u"km/s", # Example velocity scale
    n = 1.0,
    t_end = 10.0u"Gyr",
    rho_initial = rho_wrapper,
    r_max = 10.0u"kpc",
    r_min = 1.0e-3u"kpc"
)

sol = initialize_profile(dm_params)