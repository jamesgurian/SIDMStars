

function nfw_profile(r_s, rho_0, r)
    x = r / r_s
    return rho_0 / (x * (1 + x)^2)
end
# The ODE system: du/dr = f(u, p, r)
# u[1] = radius (r)
# u[2] = Pressure (P)
# u[3] = entropy (s)
# u[4] = luminosity (L)
# u[5] = ds/dt
function hse_inward!(du, u, p, M)
    r = u[1]
    P = u[2]
    #s = u[3]
    #L = u[4]
    
    # Get density at current radius r
    ρ = p.ρ_func(r)
    #v2 = P/ρ
    # 1. Mass Conservation: dr/dm = 1/ (4πr^2ρ)
    du[1] = 1/(4 * pi * r^2 * ρ )
    # 2. Hydrostatic Equilibrium: dP/dm = -GMρ/r^2 (G=1)
    du[2] = -(p.G * M  / (4*pi * r^4)) 

   # du[3] - 2.5/(p.κ(v2) * P* 16*pi^2*r^4) - du[2] / (P)
   # du[4] = -v2 * u[5]

   # du[5] = s - log(P/ρ^(5/2))

end

function initialize_profile(run_params; mass_min_frac=1e-7)
    # Setup a Radial Grid (Log-spaced is best for cores)
    rho_initial(r) = nfw_profile(run_params.r_s, run_params.rho_s,r)
    # Integrate Mass Outwards 
    # dM = 4 * pi * r^2 * ρ * dr
    M_tot, _ = quadgk(r -> 4 * pi * r^2 * rho_initial(r), 1e-3run_params.r_min, run_params.r_max)
    M_span = (M_tot,mass_min_frac * M_tot) 

    u0 = [run_params.r_max, 0.0]

    ρ_func(x) = rho_initial(x)
    G = ustrip(uconvert(u"kpc*km^2/s^2/Msun", Constants.G))
     ode_params = (ρ_func = ρ_func, 
                    G = G)
    

    f = ODEFunction(hse_inward!)
    prob = ODEProblem(hse_inward!,u0, M_span,ode_params)

    sol = solve(prob,Rodas5P())
    v2_func(m) = sol(m)[2]/(rho_initial(sol(m)[1]))
    println("v2",v2_func(1e6))
    κ_func(m) = 1/ (Kinv_smfp(v2_func(m), run_params) + 
                    Kinv_lmfp(v2_func(m), rho_initial(sol(m)[1]), run_params,G))
    L_func(m) = -(16 * pi^2 * (sol(m)[1])^4 * 
                    rho_initial(sol(m)[1])* 
                    κ_func(m) *ForwardDiff.derivative(v2_func,m)) 
                
    return sol, v2_func, κ_func, L_func
    #v2_func(x) = x-> sol[2](x)/ρ_func(x)
    #L(r) = (-4 * pi * x^2 * 
      #  dm_params.kappa(v2_func(r)) * 
     #   ForwardDiff.derivative(v2_func, r))
    #dLdm =
    #dLdm = ForwardDiff.derivative(x -> -κ * rho_initial * )

    #return r_grid, P_grid, s_grid, L_grid, M_grid, s_interp
end