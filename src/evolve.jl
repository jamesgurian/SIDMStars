

function evolve(run_params, rho_initial;eps_diff=1, eps_kh=0.05,p_outer=0,mass_min_frac=1e-6,gridpoints=50,velocity_bins=50)
    sol, s, M_tot = initialize_profile(run_params, rho_initial,p_outer=p_outer,mass_min_frac=mass_min_frac)
   
    
    M_min = mass_min_frac * M_tot
    r_min = sol(M_min)[1]
    G = ustrip(uconvert(u"kpc*km^2/s^2/Msun", Constants.G))

    M_tot -= 1e-4*M_tot #bvp has a hard time with the true zero pressure BC, because we log transform P
    p_surf = sol(M_tot)[2]
    

    function κ(P, s)
        ρ = P^(3/5) *exp(-2s/5)
        v2 = P/ρ
        k_tot=  1/ (Kinv_smfp(v2, run_params) + 
            Kinv_lmfp(v2, ρ, run_params, G))
        return k_tot
    end
    
    #initialize the mechanical and thermodynamic state
    t=0.
    t_sub = 0.
    dt = 0.
    dt_kh = 0.
    m_face = collect(range((M_min), stop = (M_tot), length=gridpoints))
    m_center = 0.5 .* (m_face[1:end-1] .+ m_face[2:end])
   # m_center = sqrt.(m_face[1:end-1] .* m_face[2:end])
    logm_center = log.(m_center)
    logm_face = log.(m_face)
    r_center = sol(m_center)[1,:]
    r_face = sol(m_face)[1,:]
    P_center = sol(m_center)[2,:]


    logP_core = log(sol(M_min)[2])
    logP_core_old = logP_core
    dm_cell  = diff(m_face)
    s_center = s.(m_center) 
    ρ_center  = P_center.^(3/5) .* exp.(-2 .* s_center./ 5)


    #may as well preallocate all the internal arrays for the entropy update
    v2_center = zeros(length(m_center))
    κ_center = zeros(length(m_center))
    κ_face = zeros(length(κ_center)-1)
    prefactor_face = zeros(length(κ_face))
    L_face = zeros(length(m_face))
    ρ_face = zeros(length(ρ_center)-1)
    dm_between_centers = zeros(length(m_center)-1)
    dL_dm = zeros(length(m_center))
    dsdt_conduct = zeros(length(m_center))
    dsdt_evap = zeros(length(m_center))
    mdot_cum = zeros(length(m_center))
    dsdt = zeros(length(m_center))
    dt_diff_array = zeros(length(prefactor_face))
    rspan = (log(M_min), log(M_tot))

    ode_params = StructureParams(s, G)
    ρ_core = M_min/(4/3 * pi * r_min^3)
    ρ_core_prev = ρ_core
    rspan = (log(M_min), log(M_tot))
    init_hdf5(run_params.snap_loc, gridpoints-1)
    sol = nothing #we'll just want this in scope later
    s_interp = nothing
    t_prev_snap = 0.
    snap_idx = 1
    write_snapshot!(run_params.snap_loc, snap_idx, t, P_center, r_center,m_center,  ρ_center, s_center, dsdt)


    #Note that this mutates the entire state, not just dsdt. We should probably make this a function rather than a closure
    #and pass everything explicitly
    function explicit_dsdt_conduction!(dsdt, s_interp, s_center, log_P_guess)
        ode_params = StructureParams(s_interp, G)
        sol = shoot_1p((log(M_min), log(M_tot)), log_P_guess, log(r_min), log(p_surf), ode_params)
        m_face .= collect(range(M_min, stop = M_tot, length = gridpoints))
        m_center .= 0.5 .* (m_face[1:end-1] .+ m_face[2:end])

        logm_face .= log.(m_face)
    
        #logm_face .= (range(log(M_min), stop = log(M_tot), length=gridpoints))
        #m_face .= exp.(logm_face)
        #m_center .= sqrt.(m_face[1:end-1] .* m_face[2:end])
        logm_center .= log.(m_center)

        P_center .= exp.(sol(logm_center)[2,:]) 

        # Calculate derived center variables
        ρ_center  .= P_center.^(3/5) .* exp.(-2 .* s_center./ 5)
        v2_center .= P_center ./ ρ_center
        κ_center  .= κ.(P_center, s_center) 
        r_face .= exp.(sol(logm_face)[1,:])
        r_center .= exp.(sol(logm_center)[1,:])

        #set r_min by assuming the mass-radius relationship is a power law and extrapolating
        slope = log(r_center[2]/r_center[1])/log(m_center[2]/m_center[1])
        r_min = r_center[1] * (M_min/m_center[1])^slope
        dm_between_centers .= diff(m_center)


        dv2_dm_face = diff(v2_center) ./ dm_between_centers
        κ_face .= 2*(κ_center[1:end-1] .* κ_center[2:end]) ./(κ_center[1:end-1] .+ κ_center[2:end]) #harmonic mean
    
        ρ_face .= sqrt.(ρ_center[1:end-1] .* ρ_center[2:end])
        prefactor_face .= 16 * pi^2 .* r_face[2:end-1].^4 .* ρ_face .* κ_face
        L_face[2:end-1] .= -prefactor_face .* dv2_dm_face

        dL_dm .= diff(L_face) ./ dm_cell
        dsdt_conduct .= -dL_dm ./ v2_center
        slope = (dsdt_conduct[3] - dsdt_conduct[2])/log(dm_between_centers[2])
        dsdt_conduct[1] = dsdt_conduct[2] - slope * log(dm_between_centers[1])

        slope = (dsdt_conduct[end] - dsdt_conduct[end-1])/log(dm_between_centers[end])
        dsdt_conduct[end] = dsdt_conduct[end-1] + slope * log(dm_between_centers[end])
        
        ρ_core = 4/3 * pi * M_min/r_min^3
        P_core = sol(log(M_min))[2]
        
        dt_diff_array .= dm_between_centers .^2 ./prefactor_face
        dt_diffusion = eps_diff *  minimum(dt_diff_array)
        return sol, P_core, ρ_core, dt_diffusion
    end


    snap_idx=1
   while t < run_params.t_end
    println("t:", t) 
        #The von Neumann Stability condition is very strict compared to the thermal timescale
        #and calculating the evaporation rate is relatively expensive
        #so we subscycle
        while t_sub <= dt_kh &&  t < run_params.t_end
            dlogp_core = (logP_core) - logP_core_old
            logP_guess = logP_core + dlogp_core
            logP_core_old = logP_core
            #we need a monotonicity preserving interpolant for convection...
            s_interp = PCHIPInterpolation(s_center, log.(m_center),extrapolation = ExtrapolationType.Linear)
            sol, logP_core, ρ_core, dt_diffusion = explicit_dsdt_conduction!(dsdt, s_interp, s_center, logP_guess)
            s_center .+= dt_diffusion .* (dsdt_conduct .+ dsdt_evap)  
            m_center .+= mdot_cum .* dt_diffusion
            M_tot += mdot_cum[end] .* dt_diffusion
            t += dt_diffusion
            t_sub += dt_diffusion

            
        end
        convect!(s_interp, s_center, v2_center, m_center)

        v2_esc_surface = 2 * (G * M_tot)/r_face[end] #neglecting the mass exterior to the integration region
        v2_esc_core =  2 * sol[end][3]
        vrange = exp.(range(log(sqrt(v2_esc_surface)), stop=log(sqrt(2*v2_esc_core)),length=velocity_bins))
        dsdt_evap, mdot_cum = calc_evaporation(run_params, G, sol, s_interp, m_face, m_center, v2_center, vrange)
        t_sub = 0.
        dt_kh = eps_kh *  minimum((abs.(1 ./ (dsdt_evap .+ dsdt_conduct))))



        if t - t_prev_snap > run_params.snap_interval_time || abs((ρ_core - ρ_core_prev) / ρ_core_prev) > run_params.snap_interval_density
            ρ_core_prev = ρ_core

            println("writing snapshot $snap_idx")
            write_snapshot!(run_params.snap_loc, snap_idx, t, P_center, r_center,m_center,  ρ_center, s_center, dsdt)
            snap_idx += 1;
            t_prev_snap = t
        end
        
  
       
   end
end
