
using SIDMStars
#using PyPlot

using QuadGK
    
r_max = 100.
prof = x -> nfw_profile(1.18, 2.73e7,  x)


function find_psurf(rho_func, r_max,G)
    M_max, _ = quadgk(r -> 4 * π * r^2 * rho_func(r), 0.0, r_max)
    P_integrand = function(r)
        m_tail, _ = quadgk(x -> 4 * π * x^2 * rho_func(x), r_max, r)
        M_enclosed = M_max + m_tail
        return rho_func(r) * G * M_enclosed / r^2
    end
    
    P_max, _ = quadgk(P_integrand, r_max, Inf)
    return P_max
end

P_surf = find_psurf(prof, r_max, 4.30091727003628e-6)


# 3. Setup Params
dm_params = RunParams(
    convection=ConvectionOff(),
    evaporation=EvaporationOff(),
    G = 4.30091727003628e-6,
    C = 0.85,
    σ0 = 2.044e-10*50, #kpc^2/msun
    w = 10., # km/s
    n = 0.0,
    t_end = 25.,
    r_max = r_max, #kpc
    r_min = 1e-3, #kpc
    snap_interval_time = 0.1,
    snap_interval_density = 0.1,
    gridpoints=200,
    snap_loc= "nfw_out.hdf5"
)
#semilogx(sol[1].t, sol[2].(sol[1].t))
#savefig("out.png")
sol = evolve(dm_params, prof,eps_kh=.05, eps_diff=0.2, mass_min_frac=1e-4, p_outer=P_surf);