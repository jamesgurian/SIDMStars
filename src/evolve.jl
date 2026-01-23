
struct StructureParams
    s::FunctionWrapper{Float64, Tuple{Float64}}
    G::Float64
end

Base.eltype(::Type{<:StructureParams}) = Float64


function evolve(run_params, rho_initial;eps=.005,p_outer=0,mass_min_frac=1e-6,gridpoints=50)
    sol, L, s, dsdt, t_kh, M_tot = initialize_profile(run_params, rho_initial,p_outer=p_outer,mass_min_frac=mass_min_frac)
    dt = eps*minimum(t_kh)
    
    M_min = mass_min_frac * M_tot
    r_min = sol(M_min)[1]
    G = ustrip(uconvert(u"kpc*km^2/s^2/Msun", Constants.G))

    M_tot *= 0.7 #bvp has a hard time with the true zero pressure BC
    p_surf = sol(M_tot)[2]



    function κ(P, s)
        ρ = P^(3/5) *exp(-2s/5)
        v2 = P/ρ
        return 1/ (Kinv_smfp(v2, run_params) + 
            Kinv_lmfp(v2, ρ, run_params, G))
    end
    
    t=0.
    m_face = ((10 .^ range(log10(M_min), stop=log10(M_tot), length=gridpoints+1)))
    m_center = 0.5 .* (m_face[1:end-1] .+ m_face[2:end])
    logP_core = log(sol(M_min)[2])
    logP_core_old = logP_core
    dm_cell  = diff(m_face)
    s_center = s.(m_center) 
    rspan = (log(M_min), log(M_tot))
    ode_params = StructureParams(s, G)
    ρ_core = M_min/(4/3 * pi * r_min^3)
    logrmax = log(sol(M_tot)[1])
    logrmax_old = logrmax
    #println("rminit", r_min)
    rspan = (log(M_min), log(M_tot))

    function explicit_dsdt(s_center, ρ_core, log_P_guess, log_R_guess)
        s_interp = AkimaInterpolation(s_center, log.(m_center),extrapolation = ExtrapolationType.Constant)
        #s_func(x) = s_interp(log(x))
        
        ode_params = StructureParams(s_interp, G)
        #r_min = ((M_min)/ (4/3 * pi * ρ_core))^(1/3)
        #println("rm", r_min)
        @time sol = shoot((log(M_min), log(M_tot)), log_P_guess, log_R_guess, log(r_min), log(p_surf), ode_params)
        P_center = exp.(sol(log.(m_center))[2,:]) 
        
        # Calculate derived center variables
        ρ_center  = P_center.^(3/5) .* exp.(-2 .* s_center./ 5)
        v2_center = P_center ./ ρ_center
        κ_center  = κ.(P_center, s_center) # Assuming kappa depends on local state
        r_face = exp.(sol(log.(m_face))[1,:])
        dm_between_centers = diff(m_center)

        dv2_dm_face = diff(v2_center) ./ dm_between_centers
        κ_face = 2*(κ_center[1:end-1] .* κ_center[2:end]) ./(κ_center[1:end-1] .+ κ_center[2:end]) # mean
    
        ρ_face  = sqrt.(ρ_center[1:end-1] .* ρ_center[2:end])
        L_face = zeros(length(m_face))
        prefactor_face = 16 * pi^2 .* r_face[2:end-1].^4 .* ρ_face .* κ_face
        L_face[2:end-1] .= -prefactor_face .* dv2_dm_face
        dL_dm = diff(L_face) ./ dm_cell
        #dL_dm[1] = 0
        #extrapolate for the first grid point
        #slope =((dL_dm[3] - dL_dm[2])/log(dm_cell[2]))# + (dL_dm[4] - dL_dm[3])/log(dm_cell[3]))
        #dL_dm[1] = dL_dm[2] - slope * log(dm_cell[1])    
        dsdt = -dL_dm ./ v2_center
        slope = (dsdt[3] - dsdt[2])/log(dm_between_centers[2])
        dsdt[1] = dsdt[2] - slope * log(dm_between_centers[1])
        loglog(exp.(sol(log.(m_center))[1,:]),ρ_center)
        return dsdt, sol(log(M_min))[2], sol(log(M_tot))[1], ρ_center[1]
    end


    #
   while t < run_params.t_end
   # for i in 1:10
           println("t:", t) 
        
        dlogp_core = (logP_core) - logP_core_old
        logP_guess = logP_core + dlogp_core
        dlogrmax = logrmax - logrmax_old
        logP_core_old = logP_core

        logRguess = logrmax + dlogrmax
        logrmax_old = logrmax
        dsdt, logP_core, logrmax, ρ_core = explicit_dsdt(s_center, ρ_core, logP_guess, logRguess)
        s_center .+= dt .* (dsdt)       #dsdt_predictor, logP_core, logrmax= explicit_dsdt(s_center, logP_guess, logrmax)
        #s_predictor = s_center .+ dsdt_predictor .* dt

        #dsdt_corrector, _, _ = explicit_dsdt(s_predictor, logP_guess, logrmax)
        
       # s_center .+= dt .* 0.5*(dsdt_predictor .+ dsdt_corrector)
        #semilogx(m_center,0.5*(dsdt_predictor .+ dsdt_corrector))
        #=s_func = x->AkimaInterpolation(s_center, log.(m_center),extrapolation = ExtrapolationType.Constant)(log(x))
        

        
        ode_params = StructureParams(s_func, G)
        
        @time sol = shoot((log(M_min), log(M_tot)), logP_guess, logrmax, log(r_min), log(p_surf), ode_params)
        P_center = exp.(sol(log.(m_center))[2,:]) 
        #P_core_old = P_core
        logP_core = log(P_center[1])
        logrmax_old = logrmax
        logrmax = sol(log(M_tot))[1]

        # Calculate derived center variables
        ρ_center  = P_center.^(3/5) .* exp.(-2 .* s_center./ 5)
        v2_center = P_center ./ ρ_center
        κ_center  = κ.(P_center, s_center) 
        r_face = exp.(sol(log.(m_face))[1,:])
        dm_between_centers = diff(m_center)

        dv2_dm_face = diff(v2_center) ./ dm_between_centers
        κ_face = 2*(κ_center[1:end-1] .* κ_center[2:end]) ./(κ_center[1:end-1] .+ κ_center[2:end]) # mean
      
        ρ_face  = sqrt.(ρ_center[1:end-1] .* ρ_center[2:end])
        L_face = zeros(length(m_face))
        prefactor_face = 16 * pi^2 .* r_face[2:end-1].^4 .* ρ_face .* κ_face
        L_face[2:end-1] .= -prefactor_face .* dv2_dm_face
        
         L_face[end] = L_face[end-1]
        dL_dm = diff(L_face) ./ dm_cell

        slope = (dL_dm[3] - dL_dm[2])/log(dm_cell[2])
        dL_dm[1] = dL_dm[2] - slope * log(dm_cell[1])    
        dsdt = -dL_dm ./ v2_center 
        s_center .+= dt .* (dsdt)
        #dsdt = explicit_dsdt(s_center, dt, log_P_guess, log_R_guess)
        #dt=1e-4
        #dt = minimum(0.4 .* dm_between_centers .^2 ./ prefactor_face)
        #dt_prev = dt=#
        t += dt

        dt = eps * minimum(abs.(s_center ./ dsdt))
        
        #loglog(m_center, v2_center)
        #loglog(exp.(sol(log.(m_center))[1,:]),ρ_center)
        #semilogx(m_center, s_center)
       # semilogx(m_center, dL_dm)
        println("dt", dt)
       
   end
    
end
#=

function evolve(run_params, rho_initial;eps=.1,p_outer=0,mass_min_frac=1e-4)
    sol_init, L_init, s_init, dsdt, t_kh, M_tot = initialize_profile(run_params, rho_initial,p_outer=p_outer,mass_min_frac=mass_min_frac)
    dt = eps*minimum(t_kh)
    println("m_tot", M_tot)
    M_min = mass_min_frac * M_tot
    r_min = sol_init(M_min)[1]
    L_min = L_init(M_min)

    G = ustrip(uconvert(u"kpc*km^2/s^2/Msun", Constants.G))
    println("L0", L_init(0))
    println("dt", dt)
    κ(v2, ρ) = 1/ (Kinv_smfp(v2, run_params) + 
            Kinv_lmfp(v2, ρ, run_params, G))
    M_tot *= 0.1 #bvp has a hard time with the true zero pressure BC
    p_surf = sol_init(M_tot)[2]
    v2_surf = p_surf/rho_initial(sol_init(M_tot)[1])
    v2_scale = sol_init(M_min)[2]/rho_initial(r_min)
    L_scale = maximum(L_init.(sol_init[1]))
    #=function bc!(residual, sol, p, t)
        P = exp(sol(M_tot)[2])
        r = exp(sol(M_min)[1])
        residual[1] = log(r/r_min)
        residual[2] = log(P/p_surf)
        residual[3] = (sol(M_tot)[3] - v2_surf)/v2_scale
        residual[4] = (sol(M_min)[4] - L_min)/L_min
        println(residual)
    end=#
    L_scale = L_min

    function bc_l!(resid_a, u_a, p)
        resid_a[1] = u_a[1] - log(r_min)
        println("residl", resid_a)
    end

    function bc_h!(resid_b, u_b, p)
        resid_b[1] = u_b[2] - log(p_surf)
        resid_b[2] = (u_b[3] - v2_surf)/v2_surf
        println("residh", resid_b)
    end
    ode_params = StructureParams(κ,s_init, L_scale, G, dt)
    rspan = (log(M_min), log(M_tot))
    println("range", rspan)
    println("dt", dt)
   # guess(m) = [[0.0,sol(0)[2], s_init(0), 0.0]
    #guess = [sol(M_min)[1],log(sol(M_min)[2]), s_init(M_min), L_init(M_min)]
    guess_T =  1.00 .*[log(sol_init(M_min)[1]),log(sol_init(M_min)[2]), sol_init(M_min)[2]/rho_initial(r_min), L_init(M_min)/L_scale]
    println(guess_T)
    println("Time to solve")
    bvp = BVP.TwoPointBVProblem(structure_system_T!, (bc_l!, bc_h!), guess_T, rspan, ode_params, bcresid_prototype=(zeros(2),zeros(2)),
        )
    #bvp = BVP.TwoPointBVProblem(structure_system_T!, 
    #    (bc_inner!, bc_outer!), 
     #   guess_T, 
      #  rspan, ode_params,bcresid_prototype = (zeros(2), zeros(2)))
    #prob = ODEProblem(structure_system_T!,guess_T, rspan)
    #sol = solve(prob)
    #bc!(resd, sol_init, ode_params, 0)
    #condition(u, t, integrator) = exp(u[1]) < r_min || exp(u[2]) < p_surf || u[3] < v2_surf
    ##affect!(integrator) = terminate!(integrator)
    #cb = DiscreteCallback(condition, affect!)
    @time sol = BVP.solve(bvp,Shooting(Rosenbrock23(autodiff=false), RobustMultiNewton(;autodiff=AutoFiniteDiff())),abstol=1e-3, reltol=1e-3)

    return sol
end

=#
