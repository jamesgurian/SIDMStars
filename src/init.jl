
# The ODE system: du/dr = f(u, p, r)
# u[1] = Mass (M)
# u[2] = Pressure (P)
function hse_inward!(du, u, ρ_func, r)
    M = u[1]
    
    # Get density at current radius r
    ρ = ρ_func(r)
    
    # 1. Mass Conservation: dM/dr = 4πr^2ρ
    du[1] = 4 * pi * r^2 * ρ
    
    # 2. Hydrostatic Equilibrium: dP/dr = -GMρ/r^2
    du[2] = -uconvert( u"pc*km^2/s^2/Msun",Constants.G) * M * ρ / r^2
end

function initialize_profile(run_params)
    # 1. Setup a Radial Grid (Log-spaced is best for cores)
    # Avoid r=0 to prevent singularities. Start at a tiny fraction of r_max.
    r_grid = exp10.(range(log10(ustrip(run_params.r_min)), log10(ustrip(run_params.r_max)), length=1000)) * unit(run_params.r_max)
    
    # 2. Evaluate Density on the grid
    # rho_initial is a FunctionWrapper, call it element-wise
    ρ_grid = run_params.rho_initial.(r_grid)
    mass_unit = unit(ρ_grid[1]*run_params.r_min^3)

    # 3. Integrate Mass Outwards (Trapezoidal Rule)
    # dM = 4 * pi * r^2 * ρ * dr
    # We do a cumulative sum.
    
    
    M_tot, _ = quadgk(r -> 4 * pi * r^2 * run_params.rho_initial(r), 1e-3run_params.r_min, run_params.r_max)
    r_span = (run_params.r_max, run_params.r_min)
    P_surf = 0.0* unit(ρ_grid[1]*run_params.w^2)
    u0 = [M_tot, P_surf]
    prob = ODEProblem(hse_inward!,u0, r_span,run_params.rho_initial)
    # 4. Integrate Pressure Inwards (Jeans Equation)
    # dP/dr = - G * M * ρ / r^2
    # Integrate from R_max (surface) down to center.
    #P_grid = similar(ρ_grid)
    sol = solve(prob,alg=Tsit5())
    #return r_grid, P_grid, s_grid, L_grid, M_grid, s_interp
end