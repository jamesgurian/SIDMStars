
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

    M_tot *= 0.8 #bvp has a hard time with the true zero pressure BC
    p_surf = sol(M_tot)[2]

    


    function bc_l!(resid_a, u_a, p)
        resid_a[1] = u_a[1] - log(r_min)
        #println(resid_a)
    end

    function bc_h!(resid_b, u_b, p)
        resid_b[1] = u_b[2] - log(p_surf)
        #println(resid_b)
    end
    function κ(P, s)
        ρ = P^(3/5) *exp(-2s/5)
        v2 = P/ρ
        return 1/ (Kinv_smfp(v2, run_params) + 
            Kinv_lmfp(v2, ρ, run_params, G))
    end
    t=0.
    m_face = 10 .^ range(log10(M_min), stop=log10(M_tot), length=gridpoints+1)
    m_center = 0.5 .* (m_face[1:end-1] .+ m_face[2:end])
    P_core = sol(M_min)[2]
    P_core_old = P_core
    dm_cell  = diff(m_face) # Mass of each shell (for the final energy step)
    s_center = s.(m_center) .+ dsdt.(m_center) .* dt 
    #loglog(m_face,abs.(L.(m_face)))
    guess =  [log(r_min),log(P_core)]
    rspan = (log(M_min), log(M_tot))
    ode_params = StructureParams(s, G)
    #prob = ODEProblem(hse_system!, guess, rspan, ode_params)
    bvp = BVP.TwoPointBVProblem(hse_system!, (bc_l!, bc_h!), guess, rspan, ode_params, bcresid_prototype=(zeros(1),zeros(1)),
            )
    println(M_tot)
    for i in 1:3
        println(t) 
        t += dt
        
        s_func = x->CubicSpline(s_center, log.(m_center),extrapolation = ExtrapolationType.Constant)(log(x))
       

        
        ode_params = StructureParams(s_func, G)
        rspan = (log(M_min), log(M_tot))
        dlogp_dt_core = (log(P_core) - log(P_core_old))
        guess =  [log(r_min),log(P_core) + dlogp_dt_core]
      #  bvp = BVP.TwoPointBVProblem(hse_system!, (bc_l!, bc_h!), guess, rspan, ode_params, bcresid_prototype=(zeros(1),zeros(1)),
       #     )
        bvp = remake(bvp, u0=guess, p=ode_params)
        #@time sol = solve(prob, alg=Rosenbrock23(autodiff=false),abstol=1e-4)
        #@time sol = BVP.solve(bvp,MultipleShooting(4,ode_alg=Rosenbrock23(autodiff=false)),odesolve_kwargs=(;verbose=false),abstol=1e-4,reltol=1e-4)
    
        @time sol = BVP.solve(bvp,Shooting(ode_alg=Rosenbrock23(autodiff=false),nlsolve=NewtonRaphson(autodiff=AutoFiniteDiff()),jac_alg=SimpleTrustRegion()),odesolve_kwargs=(;verbose=false),abstol=1e-5,reltol=1e-5)
        P_center = exp.(sol(log.(m_center))[2,:]) 
        P_core_old = P_core
        P_core = P_center[1]

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
        

        # Calculate the local power-law index `beta` between the first two calculated points
        # L2 is at M_face[2], L3 is at M_face[3]
      
        # Extrapolate to M_min (M_face[1])
        # We preserve the sign of L_face[2] (inflow or outflow)
        #L_face[1] = L_face[2] * (m_face[1] / m_face[2])^β
        L_face[end] = L_face[end-1]
        dL_dm = diff(L_face) ./ dm_cell
        # log-log slope between cell 2 and 3
        dL_dm[1] = dL_dm[2]    
        dsdt = -dL_dm ./ v2_center
        s_center .+= dsdt*dt
        #dt=1e-4
        #dt = minimum(0.4 .* dm_between_centers .^2 ./ prefactor_face)
        dt_prev = dt
        dt = eps * minimum(abs.(s_center ./ dsdt))
        
        #loglog(m_center, v2_center)
        loglog(exp.(sol(log.(m_center))[1,:]),ρ_center)
        #scatter(m_center, s_center)
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
