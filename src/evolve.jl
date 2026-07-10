
struct GridState
    M_min::Float64
    M_tot::Ref{Float64}
    p_surf::Float64
    logM_fit::Ref{Float64}

    m_face::Vector{Float64}
    m_center::Vector{Float64}
    logm_face::Vector{Float64}
    logm_center::Vector{Float64}
    dm_cell::Vector{Float64}
    dm_between_centers::Vector{Float64}

    r_face::Vector{Float64}
    r_center::Vector{Float64}
    P_center::Vector{Float64}
    s_center::Vector{Float64}
    ρ_center::Vector{Float64}
    ρ_face::Vector{Float64}
    v2_center::Vector{Float64}
    κ_center::Vector{Float64}

    prefactor_face::Vector{Float64}
    L_face::Vector{Float64}
    dL_dm::Vector{Float64}
    dsdt_conduct::Vector{Float64}
    dsdt_evap::Vector{Float64}
    mdot_cum::Vector{Float64}
    dt_diff_array::Vector{Float64}
end

function initialize_grid(gridpoints, rho_initial, run_params, p_outer, mass_min_frac)
    M_tot = quadgk(x -> 4 * pi *x^2 * rho_initial(x), run_params.r_min, run_params.r_max)[1]
    M_min = mass_min_frac*M_tot
    x = (0.5 .* (cos.(range(pi, stop=0, length = gridpoints)) .+ 1)).^(2)
    
    m_face =  M_min .+ (M_tot - M_min) .* x
    m_center = 0.5 .* (m_face[1:end-1] .+ m_face[2:end])

    sol, s = initialize_profile(run_params, rho_initial,m_center[1], M_tot,p_outer=p_outer)

    m_center = 0.5 .* (m_face[1:end-1] .+ m_face[2:end])   
    logm_center = log.(m_center)
    logm_face = log.(m_face)
    r_center = Vector{Float64}(undef, length(m_center))
    P_center = Vector{Float64}(undef, length(m_center))
   

    for i in eachindex(logm_center)
        r_center[i] = (sol(m_center[i])[1])
        P_center[i] = (sol(m_center[i])[2])
    end
    r_face = Vector{Float64}(undef,gridpoints)


    for i in 2:gridpoints
        r_face[i] = (sol(m_face[i])[1])
    end
    m_c1 = m_center[1]
    m_f1 = m_face[1]
    r_c1 = r_center[1]
   
    # Solve for r_face[1]
    r_face[1] =  r_c1 * (m_f1 / m_c1)^(1/3)
    
    dm_cell = diff(m_face)
    dm_between_centers = diff(m_center)

    s_center = s.(m_center) 
    ρ_center = P_center.^(3/5) .* exp.(-2 .* s_center./ 5)


    # Preallocate the internal arrays
    v2_center = zeros(length(m_center))
    κ_center = zeros(length(m_center))
    prefactor_face = zeros(length(κ_center)-1)
    L_face = zeros(length(m_face))
    ρ_face = zeros(length(ρ_center)-1)
    dL_dm = zeros(length(m_center))
    dsdt_conduct = zeros(length(m_center))
    dsdt_evap = zeros(length(m_center))
    mdot_cum = zeros(length(m_center))
    dt_diff_array = zeros(length(prefactor_face))

    p_surf = Float64(p_outer)
    fit_idx = searchsortedlast(P_center, sqrt(P_center[1]), rev=true)
    logM_fit = log(m_center[fit_idx])

    return GridState(M_min, Ref(M_tot), p_surf, Ref(logM_fit),
                     m_face, m_center, logm_face, logm_center, dm_cell, dm_between_centers,
                     r_face, r_center, P_center, s_center, ρ_center, ρ_face, v2_center,
                     κ_center, prefactor_face, L_face, dL_dm, dsdt_conduct,
                     dsdt_evap, mdot_cum, dt_diff_array), r_face[1], sol, s
end


#we want to construct a grid which has resolution where log(r) or s are varying rapidly

function calc_Q!(state, run_params::RunParams{ConvectionOn, T}) where T

    r_norm = (log.(state.r_center) .- log.(state.r_center[1])) ./ (log.(state.r_center[end]) - log.(state.r_center[1]))

    invP = 1 ./ state.P_center

    invP_norm = (log.(invP) .- log.(invP[1])) ./ (log.(invP[end]) - log.(invP[1]))

    s_norm = (state.s_center .- state.s_center[1]) ./ (state.s_center[end] - state.s_center[1])



    Q_center = (0.45) .* (r_norm .+ s_norm  ) .+ 0.1 .* invP_norm

end



#s is only guaranteed monotonic with convection on 

#so we need a different Q which doesn't  rely on s

function calc_Q!(state, run_params::RunParams{ConvectionOff, T}) where T

    r_norm = (log.(state.r_center) .- log.(state.r_center[1])) ./ (log.(state.r_center[end]) - log.(state.r_center[1]))

    m_norm = (log.(state.m_center) .- log.(state.m_center[1])) ./ (log.(state.m_center[end]) - log.(state.m_center[1]))



    invP = 1 ./ state.P_center

    invP_norm = (log.(invP) .- log.(invP[1])) ./ (log.(invP[end]) - log.(invP[1]))



    Q_center = (0.4) .* (r_norm .+ m_norm ) .+ 0.2 .* invP_norm

end





function update_grid!(state::GridState, s_interp, sol, gridpoints, run_params)

 

    Q_center = calc_Q!(state, run_params)

    #so we interpolate in normalized logr and s to get our new mass grid

    logm_of_Q = PCHIPInterpolation(state.logm_center, Q_center)

    Q_range = collect(range(Q_center[1], stop = Q_center[end], length=gridpoints-1))

    logm_center_new = logm_of_Q(Q_range)



    #and then update the state

    state.s_center .= s_interp(logm_center_new)

    state.logm_center .= logm_center_new

    state.m_center .= exp.(state.logm_center)



    state.m_face[2:end-1] .= 0.5 .* (state.m_center[2:end] .+ state.m_center[1:end-1])

    state.logm_face .= log.(state.m_face)



    state.dm_cell .= diff(state.m_face)

    state.dm_between_centers .= diff(state.m_center)



    for i in eachindex(state.logm_center)
        
        state.r_center[i] = exp(sol(state.logm_center[i])[1])
        state.P_center[i] = exp(sol(state.logm_center[i])[2])
    end
    for i in eachindex(state.logm_face)[2:end]
        state.r_face[i] = exp(sol(state.logm_face[i])[1])
    end

    state.logm_center .= log.(state.m_center)

    state.logm_face .= log.(state.m_face)



    state.ρ_center  .= state.P_center.^(3/5) .* exp.(-2 .* state.s_center./ 5)

    state.v2_center .= state.P_center ./ state.ρ_center







end

function core_energy(P_c, ρ_c, r_c, r_min, m_f1, G)
    P_face = P_c + (2/3)*pi*G*ρ_c^2 * (r_c^2 - r_min^2)
    V1 = (4/3)*pi*r_min^3
    U1 = 1.5*P_face*V1 + 1.5*(3/10)*G*m_f1^2/r_min
    W1 = -(3/5)*G*m_f1^2/r_min
    return U1 + W1
end

function hse_update!(state::GridState, logP_guess, s_interp, E_old, G, abstol, reltol)
    ode_params = StructureParams(s_interp, G)
    log_rmax = log(state.r_face[end])

    logP_core_old = log(state.P_center[1])
    sol, abstol, reltol, logM_fit_new = shoot_2p_fit((state.logm_center[1], log(state.M_tot[])), logP_guess, logP_core_old, log_rmax, log(state.r_center[1]), log(state.p_surf), state.logM_fit[], ode_params, abstol, reltol)
    state.logM_fit[] = logM_fit_new
    

    for i in eachindex(state.logm_center) 
        state.r_center[i], state.P_center[i]  = exp.(sol(state.logm_center[i]))
    end
    for i in eachindex(state.logm_face)
        state.r_face[i] = exp(sol(state.logm_face[i])[1])
    end
    state.ρ_center  .= state.P_center.^(3/5) .* exp.(-2 .* state.s_center./ 5)
    state.v2_center .= state.P_center ./ state.ρ_center
  
    m_c1 = state.m_center[1]
    m_f1 = state.m_face[1]
    r_c1 = state.r_center[1]
    P_c1 = state.P_center[1]

    r_min_new = r_c1 * (m_f1 / m_c1)^(1/3)
    E_new = core_energy(P_c1, state.ρ_center[1], r_c1, r_min_new, m_f1, G)
    dE_implied = E_new - E_old
    state.r_face[1] = r_min_new

    ρ_core = state.M_min/(4/3 * pi * r_min_new^3)
    logP_core = log(state.P_center[1])
    return sol, logP_core, ρ_core, abstol, reltol, E_new, dE_implied
end

function explicit_dsdt_conduction!(state::GridState, G, eps_diff, run_params)
    
    # Calculate derived center variables; one log(P) feeds ρ, v² and κ
    @inbounds @fastmath @simd for i in eachindex(state.P_center)
        P_c = state.P_center[i]
        s_c = state.s_center[i]

        lP  = log(P_c)
        lρ  = 0.6 * lP - 0.4 * s_c
        lv2 = 0.4 * (lP + s_c)
        ρ_c = exp(lρ)

        state.ρ_center[i] = ρ_c
        state.v2_center[i] = P_c / ρ_c
        state.κ_center[i] = calc_κ_log(lρ, lv2, run_params)
    end
   
   
    @turbo for j in 1:(length(state.κ_center) - 1)
        
        dv2_dm = (state.v2_center[j+1] - state.v2_center[j]) / state.dm_between_centers[j]
        
        # Averages across the face
        κ_f = 0.5 * (state.κ_center[j] + state.κ_center[j+1])
        ρ_f = sqrt(state.ρ_center[j] * state.ρ_center[j+1])
     
        
        r_f = state.r_face[j+1] 
        
        pref = 16 * pi^2 * r_f^4 * ρ_f * κ_f
        state.prefactor_face[j] = pref
        
        state.L_face[j+1] = -pref * dv2_dm * 1.5
    end
    
    @turbo for i in eachindex(state.dL_dm)
        state.dL_dm[i] = (state.L_face[i+1] - state.L_face[i]) / state.dm_cell[i]
        state.dsdt_conduct[i] = - state.dL_dm[i] / state.v2_center[i]
    end
        

    slope_conduct2 = (state.dsdt_conduct[end-1] - state.dsdt_conduct[end-2])/log(state.dm_between_centers[end-1])
    state.dsdt_conduct[end-1] + slope_conduct2 * log(state.dm_between_centers[end])

    state.dt_diff_array .= state.dm_between_centers .^2 ./ state.prefactor_face

    dt_diffusion = eps_diff * minimum(state.dt_diff_array)
    
    return dt_diffusion
end

function evolve(run_params, rho_initial;eps_diff=1, eps_kh=0.05,p_outer=0,mass_min_frac=1e-6,velocity_bins=50, abstol=1e-5, reltol=1e-5)
    gridpoints = run_params.gridpoints

    G = run_params.G

    
    
    state, r_min, sol, s = initialize_grid(gridpoints, rho_initial, run_params, p_outer, mass_min_frac)

    M_tot_init = state.M_tot[]
    collapse_m_target = run_params.collapse_mass_frac * M_tot_init
    collapse_idx = argmin(abs.(state.m_center .- collapse_m_target))
    ρ_collapse_init = state.ρ_center[collapse_idx]

    logP_core = log(sol(state.m_center[1])[2])
    logP_core_old = logP_core
    ρ_core = state.M_min/(4/3 * pi * r_min^3)
    ρ_core_prev = ρ_core
    init_hdf5(run_params.snap_loc, gridpoints-1, run_params)
    t=0.
    t_sub = 0.
    dt_kh = 0.
    t_prev_snap = 0.
    dt_diffusion= 0.

    s_interp = PCHIPInterpolation(state.s_center, state.logm_center, extrapolation = ExtrapolationType.Constant)
    E_core = core_energy(state.P_center[1], state.ρ_center[1], state.r_center[1], state.r_face[1], state.m_face[1], G)
    sol, logP_core, ρ_core, abstol, reltol, E_core, dE_implied =
        hse_update!(state, logP_core, s_interp, E_core, G, abstol, reltol)
    
    
        dt_diffusion = explicit_dsdt_conduction!(state,  G, eps_diff, run_params)
    #update_grid!(state, s_interp, sol, gridpoints, run_params)

        snap_idx = 1
    if !isinf(run_params.snap_interval_time)
        write_snapshot!(run_params.snap_loc, snap_idx, t, state)
    end
    snap_idx = 2
    logP_guess = log(state.P_center[1])

    while t < run_params.t_end
    println("t:", t) 
        #The von Neumann Stability condition is very strict compared to the thermal timescale
        #and calculating the evaporation rate is relatively expensive
        #so we subscycle

        while t_sub <= dt_kh &&  t < run_params.t_end && 
            (t - t_prev_snap <= run_params.snap_interval_time && 
            abs((ρ_core - ρ_core_prev) / ρ_core_prev) <= run_params.snap_interval_density)
            
            t += dt_diffusion
            t_sub += dt_diffusion
            state.s_center .+= dt_diffusion .* (state.dsdt_conduct .+ state.dsdt_evap)  
            state.m_center .+= state.mdot_cum .* dt_diffusion
            state.M_tot[] += state.mdot_cum[end] * dt_diffusion

            state.m_face[2:end-1] .= 0.5 .* (state.m_center[2:end] .+ state.m_center[1:end-1])
            state.m_face[end] = state.M_tot[]
            state.dm_cell  .= diff(state.m_face)
            state.dm_between_centers .= diff(state.m_center)

            state.logm_face .= log.(state.m_face)  
            state.logm_center .= log.(state.m_center) 
            dt_diffusion = explicit_dsdt_conduction!(state,G, eps_diff, run_params)
            
            

            

            
        end

        #at this point, the conduction and evaporation rates are on different grids
        #so we can't add them
        dt_kh = eps_kh *  min(minimum((abs.(1 ./ (state.dsdt_evap)))), minimum(abs.(1 ./state.dsdt_conduct)))
        dlogp_core = logP_core - logP_core_old
        logP_guess = log(state.P_center[1]) + dlogp_core
        logP_core_old = logP_core 
        sol, logP_core, ρ_core, abstol, reltol, E_core, dE_implied =
            hse_update!(state, logP_guess, s_interp, E_core, G, abstol, reltol)
        state.L_face[1] =  dE_implied/t_sub
        convect!(s_interp, state.s_center, state.v2_center, state.m_center, run_params)
       # update_grid!(state, s_interp, sol, gridpoints, run_params)

        v2_esc_surface = 2 * (G * state.M_tot[]) / state.r_face[end] # neglecting the mass exterior to the integration region
        v2_esc_core = 2 * sol[end][3]
        vrange = exp.(range(log(sqrt(v2_esc_surface)), stop=log(sqrt(2*v2_esc_core)),length=velocity_bins))
        dsdt_evap_new, mdot_cum_new = calc_evaporation(run_params, G, sol, s_interp, state.m_face, state.m_center, state.v2_center, vrange)
        
        state.dsdt_evap .= dsdt_evap_new
        state.mdot_cum .= mdot_cum_new
        t_sub = 0.

        if t - t_prev_snap > run_params.snap_interval_time || abs((ρ_core - ρ_core_prev) / ρ_core_prev) > run_params.snap_interval_density
            ρ_core_prev = ρ_core
            println("writing snapshot $snap_idx")
            write_snapshot!(run_params.snap_loc, snap_idx, t, state)
            snap_idx += 1;
            t_prev_snap = t
        end

        if !isinf(run_params.collapse_density_factor)
            collapse_idx = argmin(abs.(state.m_center .- collapse_m_target))
            if state.ρ_center[collapse_idx] / ρ_collapse_init >= run_params.collapse_density_factor
                println("collapse termination: density factor $(run_params.collapse_density_factor) reached at t=$t")
                write_snapshot!(run_params.snap_loc, snap_idx, t, state)
                snap_idx += 1
                break
            end
        end

   end
end
