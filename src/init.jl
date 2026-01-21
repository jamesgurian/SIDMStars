

function nfw_profile(r_s, rho_0, r)
    x = r / r_s
    return rho_0 / (x * (1 + x)^2)
end



function initialize_profile(run_params, rho_initial; p_outer=1e-10, mass_min_frac=1e-7)
    # Setup a Radial Grid 
    # Integrate Mass Outwards 
    # dM = 4 * pi * r^2 * ρ * dr
    safety_factor=1 #integrate a bit further, in case the bc was weird
    M_tot, _ = quadgk(r -> 4 * pi * r^2 * rho_initial(r), 1e-3run_params.r_min, run_params.r_max*safety_factor)
    M_span = ((M_tot),(mass_min_frac*M_tot/2)) 

    u0 = [(run_params.r_max*safety_factor), (p_outer)]

    G = ustrip(uconvert(u"kpc*km^2/s^2/Msun", Constants.G))
    ode_params = (ρ_func = rho_initial, 
                    G = G)
   
    prob = ODEProblem(hse_system_rho!,u0, M_span, ode_params)

    sol = solve(prob,Rodas4P2(autodiff=false),reltol=1e-12)
    v2_func(m) = sol(m)[2]/(rho_initial(sol(m)[1]))

    dv2(m) = 1/rho_initial(sol(m)[1]) * sol(m, Val{1})[2] - 
            sol(m)[2]/(4 * pi * sol(m)[1]^2 * rho_initial(sol(m)[1])^3) * ForwardDiff.derivative(rho_initial, sol(m)[1])
    κ_func(m) = 1/ (Kinv_smfp(v2_func(m), run_params) + 
                    Kinv_lmfp(v2_func(m), rho_initial(sol(m)[1]), run_params,G))
    L_func(m) = -(16 * pi^2 * (sol(m)[1])^4 * 
                    rho_initial(sol(m)[1])* 
                    κ_func(m) *dv2(m)) 

    s_func(m) = log(sol(m)[2]^(3/2)/rho_initial(sol(m)[1])^(5/2)) 
    dsdt(m) = -ForwardDiff.derivative(L_func,m)/(v2_func(m))
    t_kh = abs.(sol.t[2:end] .* v2_func.(sol.t[2:end]) ./ L_func.(sol.t[2:end]))
    return sol, L_func, s_func, dsdt, t_kh, M_tot
end

