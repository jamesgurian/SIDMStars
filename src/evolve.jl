
struct StructureParams{T,U}
    κ::T
    s_prev::U
    L_scale::Float64
    G::Float64
    dt::Float64
end

Base.eltype(::Type{<:StructureParams}) = Float64


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
        resid_a[2] = u_a[4] - L_min/L_scale
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


