
struct StructureParams
    s::FunctionWrapper{Float64, Tuple{Float64}}
    G::Float64
end

Base.eltype(::Type{<:StructureParams}) = Float64


function evolve(run_params, rho_initial;eps=.005,p_outer=0,mass_min_frac=1e-6,gridpoints=50)
    sol, L, s, dsdt, t_kh, M_tot = initialize_profile(run_params, rho_initial,p_outer=p_outer,mass_min_frac=mass_min_frac)
    dt = 0;
    
    M_min = mass_min_frac * M_tot
    r_min = sol(M_min)[1]
    G = ustrip(uconvert(u"kpc*km^2/s^2/Msun", Constants.G))

    M_tot -= 1e-4*M_tot #bvp has a hard time with the true zero pressure BC, because we log transform P
    p_surf = sol(M_tot)[2]
    println("surf", p_surf)
    

    function κ(P, s)
        ρ = P^(3/5) *exp(-2s/5)
        v2 = P/ρ
        k_tot=  1/ (Kinv_smfp(v2, run_params) + 
            Kinv_lmfp(v2, ρ, run_params, G))
        return k_tot
    end
    
    t=0.
    m_face = 10 .^ range(log10(M_min), stop = log10(M_tot), length=gridpoints)

    m_center = sqrt.(m_face[1:end-1] .* m_face[2:end])
    r_center = sol(m_center)[1,:]
    logP_core = log(sol(M_min)[2])
    logP_core_old = logP_core
    dm_cell  = diff(m_face)
    s_center = s.(m_center) 
    rspan = (log(M_min), log(M_tot))
    ode_params = StructureParams(s, G)
    ρ_core = M_min/(4/3 * pi * r_min^3)
    logrmax = log(sol(M_tot)[1])
    println("rtot, ", sol(M_tot)[1])
    println("rmin", sol(M_min)[1])
    println("mtot", log(M_tot))
    logrmax_old = logrmax
    rspan = (log(M_min), log(M_tot))
    snap_t = zeros(length(run_params.snaptimes))
    snap_p = zeros(length(run_params.snaptimes), length(m_center))
    snap_ρ = zeros(length(run_params.snaptimes),length(m_center))
    snap_r = zeros(length(run_params.snaptimes),length(m_center))
    snap_s = zeros(length(run_params.snaptimes),length(m_center))
    snap_dict = Dict(["t"=>snap_t, "P"=>snap_p, "ρ" => snap_ρ, "r"=>snap_r, "s"=>snap_s])
    function explicit_dsdt(s_center, r_center, ρ_core, log_P_guess, log_R_guess)
        s_interp = AkimaInterpolation(s_center, log.(m_center),extrapolation = ExtrapolationType.Linear)
        ode_params = StructureParams(s_interp, G)
       # loglog(m_center,r_center)
        #@time sol = shoot((log(M_min), log(M_tot)), log_P_guess, log_R_guess, log(r_min), log(p_surf), ode_params)
        sol = shoot_1p((log(M_min), log(M_tot)), log_P_guess, log(r_min), log(p_surf), ode_params)
        
        #  scatter(exp(sol.t[end]), exp(sol[end][1]))
       
        P_center = exp.(sol(log.(m_center))[2,:]) 

        # Calculate derived center variables
        ρ_center  = P_center.^(3/5) .* exp.(-2 .* s_center./ 5)
        v2_center = P_center ./ ρ_center
        κ_center  = κ.(P_center, s_center) 
        r_face = exp.(sol(log.(m_face))[1,:])
       #loglog(r_center, ρ_center)
      # semilogx(exp.(sol.t), (sol[1,:]))
        #println("logrmin")
   #     println("rtot", log(r_face[end]))
        r_center = exp.(sol(log.(m_center))[1,:])
               # loglog(m_center,r_center)

        slope = log(r_center[2]/r_center[1])/log(m_center[2]/m_center[1])
        r_min = r_center[1] * (M_min/m_center[1])^slope
        dm_between_centers = diff(m_center)

        dv2_dm_face = diff(v2_center) ./ dm_between_centers
        κ_face = 2*(κ_center[1:end-1] .* κ_center[2:end]) ./(κ_center[1:end-1] .+ κ_center[2:end]) # mean
    
        ρ_face  = sqrt.(ρ_center[1:end-1] .* ρ_center[2:end])
        L_face = zeros(length(m_face))
        prefactor_face = 16 * pi^2 .* r_face[2:end-1].^4 .* ρ_face .* κ_face
        L_face[2:end-1] .= -prefactor_face .* dv2_dm_face
              #  semilogx(r_face, L_face)

        dL_dm = diff(L_face) ./ dm_cell
        dsdt_conduct = -dL_dm ./ v2_center
        slope = (dsdt_conduct[3] - dsdt_conduct[2])/log(dm_between_centers[2])
        dsdt_conduct[1] = dsdt_conduct[2] - slope * log(dm_between_centers[1])

        slope = (dsdt_conduct[end] - dsdt_conduct[end-1])/log(dm_between_centers[end])
        dsdt_conduct[end] = dsdt_conduct[end-1] + slope * log(dm_between_centers[end])
        slope = log(ρ_center[2]/ρ_center[1])/log(r_center[2]/r_center[1])
       # v2_esc_surface = 2 * (G * M_tot)/r_face[end] #neglectting the mass exterior to the integration region
       # v2_esc_core =  2 * sol[end][3]
       # println("v2_surf", v2_esc_surface )
        #grid is an odd size so we'll get errors if we broadcast something we shouldn't 
        #vrange = exp.(range(log(sqrt(v2_esc_surface)), stop=log(sqrt.(2*v2_esc_core)),length=50))
        #
        #Φ = calc_Φ(G, m_center, M_tot, exp(sol(log(M_tot))[1]),r_center, diff(r_face))
       # v2_esc_surface = 2 * (G * M_tot)/r_face[end] #neglectting the mass exterior to the integration region
        #v2_esc_core = -2 * Φ[1]
        #grid is an odd size so we'll get errors if we broadcast something we shouldn't 
        
        #vrange = exp.(range(log(.001), stop=log(sqrt.(v2_esc_core)),length=2*gridpoints+3))
        #calc_τ(run_params,sol, s_interp,m_face, m_center, sqrt(v2_esc_surface))
    

        #sdot_evap = calc_evaporation(run_params, Φ, diff(m_face), r_center, diff(r_face), ρ_center, v2_center, vrange)
        #dsdt = dsdt_conduct #.+ sdot_evap
        dsdt_tot =  dsdt_conduct
        ρ_core =  ρ_center[1] * (r_min[1]/r_center[1])^slope #extrapolate to r_min
        return dsdt_tot, r_center, sol(log(M_min))[2], sol(log(M_tot))[1], ρ_core, sol
    end


    snap_idx=1
  while t < run_params.t_end
    #
           println("t:", t) 
        
        dlogp_core = (logP_core) - logP_core_old
        logP_guess = logP_core + dlogp_core
        dlogrmax = logrmax - logrmax_old
        logP_core_old = logP_core

        logRguess = logrmax 
        logrmax_old = logrmax
        dsdt,r_center, logP_core, logrmax, ρ_core, sol = explicit_dsdt(s_center, r_center, ρ_core, logP_guess, logRguess)
        s_center .+= dt .* (dsdt)       
        
        t += dt
        if t > run_params.snaptimes[snap_idx] && snap_idx <= length(run_params.snaptimes)
            pvals = exp.(sol(log.(m_center))[2,:])
            rvals = exp.(sol(log.(m_center))[1,:])
            ρ_vals =  pvals.^(3/5) .* exp.(-2 .* s_center./ 5)
            write_snapshot!(snap_dict, snap_idx,t, pvals, rvals, ρ_vals, s_center)
            snap_idx += 1;
        end
        dt = eps * (1/maximum(abs.(dsdt ./ s_center)))
        println(dt)
  
       
   end
   return snap_dict
end

function write_snapshot!(snap_dict, i,t,p, r, ρ,s)
    snap_dict["t"][i,:] .= t
    snap_dict["P"][i,:] .= p
    snap_dict["r"][i,:] .= r
    snap_dict["ρ"][i,:] .= ρ
    snap_dict["s"][i,:] .= s
end