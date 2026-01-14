function evolve(run_params, rho_initial;eps=.001,p_outer=0,mass_min_frac=1e-4)
    sol_init, L_init, s_init, dsdt, t_kh, M_tot = initialize_profile(run_params, rho_initial,p_outer=p_outer,mass_min_frac=mass_min_frac)
    dt = eps*minimum(t_kh)
    #s_outer = log(p_outer^(3/2)/rho_initial(run_params.r_max)^(5/2))
    println("P_out", p_outer)
    M_min = mass_min_frac * M_tot
    r_min = sol_init(M_min)[1]
    L_min = L_init(M_min)

    G = ustrip(uconvert(u"kpc*km^2/s^2/Msun", Constants.G))
   # println("s error", s_init(M_tot)- s_outer)
    println("L0", L_init(0))
    println()
    κ(v2, ρ) = 1/ (Kinv_smfp(v2, run_params) + 
            Kinv_lmfp(v2, ρ, run_params, G))
    M_tot /= 2
    p_surf = sol_init(M_tot)[2]
    v2_surf = p_surf/rho_initial(sol_init(M_tot)[1])

    function bc!(residual, sol, p, t)
        residual[1] = (sol(M_min)[1] - r_min)/r_min
        residual[2] = (sol(M_tot)[2] - p_surf)/sol_init(M_min)[2]
        residual[3] = (sol(M_tot)[3] - v2_surf)/(sol(M_min)[2]/rho_initial(r_min))
        residual[4] = (sol(M_min)[4] - L_min)/L_min
        println(residual)
    end
    ode_params = (dsdt= dsdt, s_prev=s_init, 
                    G = G,κ=κ,dt=dt)
    rspan = (M_min, M_tot)
    println("range", rspan)
    println("dt", dt)
   # guess(m) = [[0.0,sol(0)[2], s_init(0), 0.0]
    #guess = [sol(M_min)[1],log(sol(M_min)[2]), s_init(M_min), L_init(M_min)]
    guess_T =  1.01 .*[sol_init(M_min)[1],(sol_init(M_min)[2]), sol_init(M_min)[2]/rho_initial(r_min), L_init(M_min)]

    println("Time to solve")
    bvp = BVP.BVProblem(structure_system_T!, bc!, guess_T, rspan, ode_params)
    #bvp = BVP.TwoPointBVProblem(structure_system_T!, 
    #    (bc_inner!, bc_outer!), 
     #   guess_T, 
      #  rspan, ode_params,bcresid_prototype = (zeros(2), zeros(2)))
    prob = ODEProblem(structure_system_T!,guess_T, rspan, ode_params)
   # sol = solve(prob)
    resd = zeros(4)
    #bc!(resd, sol_init, ode_params, 0)
    println("resid", resd)
    sol = BVP.solve(bvp,Shooting(Rodas4P2(autodiff=false)))
    return sol
end


