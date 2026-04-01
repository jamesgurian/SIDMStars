
using Revise
using SIDMStars
#using PyPlot

using Unitful
using UnitfulAstro
using DifferentialEquations
using PyPlot
m_particle = 1u"GeV/c^2"
T = 100.0u"K"      # Temperature
n_c = 1e3u"cm^-3"  # Central number density
μ = 1.0            # Mean molecular weight
m_gas = μ * m_particle

const G_cgs = ustrip(u"cm^3 * g^-1 * s^-2", Constants.G)
const k_B_cgs = ustrip(u"erg/K", Constants.k_B)
const T_cgs = ustrip(u"K", T)
const m_gas_cgs = ustrip(u"g", m_gas)
function hse_eq!(du, u, p, r)
    m, ρ = u
    dpdrho = k_B_cgs * T_cgs / m_gas_cgs  # Isothermal equation of state: P = n k_B T
    du[1] = 4π * ρ * r^2                     # dm/dr = 4πr²ρ
    du[2] = -G_cgs * m * ρ / (r^2 * dpdrho)  # dρ/dr from hydrostatic equilibrium
end
# Initial conditions
u0 = [0.0, ustrip(u"g/cm^3", m_gas * n_c)]
r_span = (1e10, 1e20)  # r/a (start just above 0 to avoid singularity)

# Solve the ODE
prob = ODEProblem(hse_eq!, u0, r_span)
sol = solve(prob; abstol=1e-8, reltol=1e-8)


# 3. Setup Params
dm_params = RunParams(
    convection=ConvectionOff(),
    evaporation=EvaporationOff(),
    G = 4.30091727003628e-6,
    C = 0.84,
    σ0 = 1e-2, #kpc^2/msun
    w = 10., # km/s
    n = 0.0,
    t_end = 25.,
    r_max = 0.01, #kpc
    r_min = 1e-4, #kpc
    snap_interval_time = 0.1,
    snap_interval_density = 0.1,
    gridpoints=200,
    snap_loc= "run_isothermal.hdf5"
)
prof_init(x) = ustrip(u"Msun/kpc^3", sol(ustrip(u"cm", x*u"kpc"))[2]*u"g/cm^3")
P_surf = sol(ustrip(u"cm", 0.01u"kpc"))[2]*u"g/cm^3" * Constants.k_B * T /m_gas
println(uconvert(u"km^2/s^2 * Msun/kpc^3", P_surf))
rvals = 10 .^ range(log10(1e-4), log10(0.01),length=40)
loglog(rvals, prof_init.(rvals))
sol_init = initialize_profile(dm_params,prof_init, p_outer=ustrip(u"km^2/s^2 * Msun/kpc^3", P_surf));
sol = evolve(dm_params, prof_init ,eps_kh=.05, eps_diff=0.2, mass_min_frac=1e-5, p_outer=ustrip(u"km^2/s^2 * Msun/kpc^3", P_surf));