
using SIDMStars

    
r_max = 30.
prof = x -> nfw_profile(1.18, 2.73e7,  x)


P_surf = calc_psurf(prof, r_max, 4.30091727003628e-6)


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
sol = evolve(dm_params, prof,eps_kh=.05, eps_diff=0.2, mass_min_frac=1e-6, p_outer=P_surf);