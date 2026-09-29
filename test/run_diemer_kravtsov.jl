
using Revise
using SIDMStars
using PyPlot

using QuadGK


# Surface pressure from the weight of the (near-zero-density) material
# overlying r_max, exactly as in run_nfw.jl. Leaving p_outer at its default
# (0 -> log(0) = -Inf) starves the outer few mass shells of any confining
# pressure and produces a non-monotonic entropy profile right at the surface,
# which is what was breaking the HSE shooting solve.
function find_psurf(rho_func, r_max, G)
    M_max, _ = quadgk(r -> 4 * π * r^2 * rho_func(r), 0.0, r_max)
    P_integrand = function(r)
        m_tail, _ = quadgk(x -> 4 * π * x^2 * rho_func(x), r_max, r)
        M_enclosed = M_max + m_tail
        return rho_func(r) * G * M_enclosed / r^2
    end

    P_max, _ = quadgk(P_integrand, r_max, Inf)
    return P_max
end

# 3. Setup Params
dm_params = RunParams(
    C = 1.0,
    G = 4.30091727003628e-6,
    σ0 = 1.04e-6, #kpc^2/M_sun
    w = 10., # km/s
    n = -4.,
    t_end = 10.,
    r_max = 2., #kpc
    r_min = 1e-3, #kpc
    snap_interval_time = 0.1,
    snap_interval_density = 0.1,
    gridpoints=150,
    snap_loc="dk_out.hdf5"
)

profile = setup_dk_profile(5e6, 7, 1.0, z=10)
P_surf = find_psurf(profile, dm_params.r_max, dm_params.G)
#eps_diff: inverse number of timesteps per von Neumann constraint on dt
#eps_kh: inverse number of timesteps per 1/\dot s
evolve(dm_params,profile,eps_diff=0.2, eps_kh=.001, mass_min_frac=1e-5, p_outer=P_surf);
snap = read_snapshot_file("dk_out.hdf5")
