



# 3. Setup Params
dm_params = RunParams(
    C = 1.0,
    σ0 = 1.04e-6, #kpc^2/M_sun
    w = 10, # km/s
    n = -4.,
    t_end = .01,
    r_max = 2, #kpc
    r_min = 1e-3, #kpc
    snap_interval_time = Inf,
    snap_interval_density = Inf,
    snap_loc="init_out.hdf5"
)

profile = setup_dk_profile(1e6, 7, 1, z=10)
evolve(dm_params,profile,gridpoints=100,eps_diff=1., eps_kh=0.1, mass_min_frac=1e-4);
