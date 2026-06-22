

function nfw_profile(r_s, rho_0, r)
    x = r / r_s
    return rho_0 / (x * (1 + x)^2)
end

function dk_profile(r_s, rho_0, rvir, nu, r)
    f_trans = (1+1/5*(r/rvir)^4)^(-nu)
    nfw_profile(r_s, rho_0, r) * f_trans
end
using Printf

function get_scale_rho_and_r_from_Mvir_cvir(M_vir, c_vir; 
                                            z=0.0, 
                                            delta_vir=200.0, 
                                            H0=67.74,    # km/s/Mpc
                                            Om0=0.3089)  # Omega matter
    
    # --- Constants & Conversions ---
    # G in (kpc * km^2) / (s^2 * M_sun)
    # Derived from G = 6.674e-11 m^3/kg/s^2
    G_astro = 4.301e-6 
    
    # Convert H0 from km/s/Mpc to km/s/kpc 
    H0_kpc = H0 / 1000.0
    
    # --- Cosmology Calculations ---
    # Critical density at z=0: rho_crit_0 = 3H^2 / 8πG
    # Resulting units: M_sun / kpc^3
    rho_crit_0 = (3 * H0_kpc^2) / (8 * π * G_astro)
    
    # Physical Matter Density at redshift z: rho_m(z)
    # rho_m(z) = rho_crit_0 * Om0 * (1+z)^3
    rho_m_kpc = rho_crit_0 * Om0 * (1 + z)^3
    
    # Convert rho_m to M_sun / pc^3 for final output/consistency with source
    # 1 kpc^3 = 10^9 pc^3

    # --- Step 1: Calculate Virial Radius (R_vir) ---
    # R_vir = (3 M / 4 π Δ ρ)^(1/3)
    # Using rho_m_kpc yields R_vir in kpc
    R_vir = (3 * M_vir / (4 * π * delta_vir * rho_m_kpc))^(1/3)

    # --- Step 2: Calculate Scale Radius (r_s) ---
    # Units: kpc
    r_s = R_vir / c_vir

    # --- Step 3: Calculate Characteristic Density (rho_s) ---
    bracket_term = log(1 + c_vir) - c_vir / (1 + c_vir)
    
    # Calculate rho_s in M_sun / kpc^3 first (since r_s is in kpc)
    rho_s_kpc = M_vir / (4 * π * r_s^3 * bracket_term)
    
    # Convert to M_sun / pc^3

    # --- Output ---
    @printf("M_vir: %0.1e, c_vir: %0.2f, z: %0.2f\n", M_vir, c_vir, z)
    println("----------------------------------------")
    println("rho_s: $rho_s_kpc")
    println("r_s: $r_s")
    println("R_vir: $R_vir\n")

    return rho_s_kpc, r_s, R_vir
end

function setup_dk_profile(M_vir, c_vir, nu; kwargs...)
    rho_s, r_s, R_vir = get_scale_rho_and_r_from_Mvir_cvir(M_vir, c_vir; kwargs...)
    x -> dk_profile(r_s, rho_s, R_vir, nu, x)
end

function initialize_profile(run_params, rho_initial; p_outer=1e-10, mass_min_frac=1e-7)
    # Setup a Radial Grid 
    # Integrate Mass Outwards 
    # dM = 4 * pi * r^2 * ρ * dr
    safety_factor=1 #integrate a bit further, in case the bc was weird
        
    r_start = 1e-3 * run_params.r_min
    r_end = run_params.r_max * safety_factor

    # 2. Calculate total mass
    M_tot, _ = quadgk(r -> 4 * pi * r^2 * rho_initial(r), r_start, r_end)

    # 3. Find r_core using Roots.jl
    M_target = mass_min_frac * M_tot

    # Objective function: difference between enclosed mass at r and the target mass
    mass_diff(r) = quadgk(x -> 4 * pi * x^2 * rho_initial(x), r_start, r)[1] - M_target

    # find_zero will quickly find the root within the bracketing interval
    #r_core = find_zero(mass_diff, (r_start, r_end))

    # 4. Create the modified, continuous density profile
   # rho_core = rho_initial(r_core)
   # rho_modified(r) = r > r_core ? rho_initial(r) : rho_core

    # 5. Set up and solve the ODE using the modified density
    M_span = (M_tot, M_target) 
    u0 = [r_end, p_outer] 
    G = run_params.G

    # Pass our new rho_modified function instead of rho_initial
    ode_params = (ρ_func = rho_initial, G = G)

    prob = ODEProblem(hse_system_rho!, u0, M_span, ode_params)
    sol = solve(prob, reltol=1e-8)
    v2_func(m) = sol(m)[2]/(rho_modified(sol(m)[1]))

    dv2(m) = 1/rho_modified(sol(m)[1]) * sol(m, Val{1})[2] - 
            sol(m)[2]/(4 * pi * sol(m)[1]^2 *rho_modified(sol(m)[1])^3) * ForwardDiff.derivative(rho_initial, sol(m)[1])
    κ_func(m) = 1/ (Kinv_smfp(v2_func(m), run_params)^(1. ./ 3.) + 
                    Kinv_lmfp(v2_func(m), rho_modified(sol(m)[1]), run_params,G)^(1. ./ 3.))^(3.)
    L_func(m) = -(16 * pi^2 * (sol(m)[1])^4 * 
                    rho_modified(sol(m)[1])* 
                    κ_func(m) *dv2(m)) 

    s_func(m) = log(sol(m)[2]^(3/2)/rho_modified(sol(m)[1])^(5/2)) 
    dsdt(m) = -ForwardDiff.derivative(L_func,m)/(v2_func(m))
    return sol, s_func, M_tot
end


function initialize_profile(run_params, rho_initial, M_min, M_tot; p_outer=1e-10, mass_min_frac=1e-7)
    # Setup a Radial Grid 
    # Integrate Mass Outwards 
    # dM = 4 * pi * r^2 * ρ * dr
    safety_factor=1 #integrate a bit further, in case the bc was weird
        
    r_start = 1e-3 * run_params.r_min
    r_end = run_params.r_max * safety_factor

    # 2. Calculate total mass
   # M_tot, _ = quadgk(r -> 4 * pi * r^2 * rho_initial(r), r_start, r_end)

    # 3. Find r_core using Roots.jl
   # M_target = mass_min_frac * M_tot

    # Objective function: difference between enclosed mass at r and the target mass
    # quadgk returns a tuple (integral, error), so we grab [1]
    mass_diff(r) = quadgk(x -> 4 * pi * x^2 * rho_initial(x), r_start, r)[1] - M_min

    # find_zero will quickly find the root within the bracketing interval
    r_core = find_zero(mass_diff, (r_start, r_end))

    # 4. Create the modified, continuous density profile
    rho_core = rho_initial(r_core)
    rho_modified(r) = r > r_core ? rho_initial(r) : rho_core

    # 5. Set up and solve the ODE using the modified density
    M_span = (M_tot, M_min)
    u0 = SVector(r_end, Float64(p_outer))
    G = run_params.G

    ode_params = (ρ_func = rho_modified, G = G)

    prob = ODEProblem(hse_system_rho, u0, M_span, ode_params)
    sol = solve(prob, reltol=1e-8)
    v2_func(m) = sol(m)[2]/(rho_modified(sol(m)[1]))

    dv2(m) = 1/rho_modified(sol(m)[1]) * sol(m, Val{1})[2] - 
            sol(m)[2]/(4 * pi * sol(m)[1]^2 *rho_modified(sol(m)[1])^3) * ForwardDiff.derivative(rho_initial, sol(m)[1])
    κ_func(m) = 1/ (Kinv_smfp(v2_func(m), run_params)^(1. ./ 3.) + 
                    Kinv_lmfp(v2_func(m), rho_modified(sol(m)[1]), run_params,G)^(1. ./ 3.))^(3.)
    L_func(m) = -(16 * pi^2 * (sol(m)[1])^4 * 
                    rho_modified(sol(m)[1])* 
                    κ_func(m) *dv2(m)) 

    s_func(m) = log(sol(m)[2]^(3/2)/rho_modified(sol(m)[1])^(5/2)) 
    dsdt(m) = -ForwardDiff.derivative(L_func,m)/(v2_func(m))
    return sol, s_func
end
#=
function initialize_profile(run_params, rho_initial, M_tot; p_outer=1e-10, mass_min_frac=1e-7)
    # Setup a Radial Grid 
    # Integrate Mass Outwards 
    # dM = 4 * pi * r^2 * ρ * dr
    M_span = ((M_tot),(mass_min_frac*M_tot/2)) 
    function enclosed_mass(R)
        integrand(r) = 4π * r^2 * rho_initial(r)
        mass, error = quadgk(integrand, 0.0, R)
        return mass
    end
    r_max = find_zero(x -> enclosed_mass(x) - M_tot, (run_params.r_min, run_params.r_max))
    u0 = [r_max, (p_outer)]

    G = run_params.G
    ode_params = (ρ_func = rho_initial, 
                    G = G)
   
    prob = ODEProblem(hse_system_rho!,u0, M_span, ode_params)

    sol = solve(prob,reltol=1e-8)
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
    
    return sol, s_func
end=#