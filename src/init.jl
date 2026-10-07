

"""
    nfw_profile(r_s, rho_0, r)

NFW density `rho_0 / (x (1+x)^2)` at radius `r`, with `x = r/r_s`.
"""
function nfw_profile(r_s, rho_0, r)
    x = r / r_s
    return rho_0 / (x * (1 + x)^2)
end

"""
    dk_profile(r_s, rho_0, rvir, nu, r)

NFW density times the truncation factor `(1 + (r/rvir)^4 / 5)^(-nu)`.
"""
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

    # --- Step 1: Calculate Virial Radius (R_vir) ---
    # R_vir = (3 M / 4 π Δ ρ)^(1/3)
    # Using rho_m_kpc yields R_vir in kpc
    R_vir = (3 * M_vir / (4 * π * delta_vir * rho_m_kpc))^(1/3)

    # --- Step 2: Calculate Scale Radius (r_s) ---
    # Units: kpc
    r_s = R_vir / c_vir

    # --- Step 3: Calculate Characteristic Density (rho_s) ---
    bracket_term = log(1 + c_vir) - c_vir / (1 + c_vir)
    
    # rho_s in M_sun / kpc^3 (r_s is in kpc)
    rho_s_kpc = M_vir / (4 * π * r_s^3 * bracket_term)

    # --- Output ---
    @printf("M_vir: %0.1e, c_vir: %0.2f, z: %0.2f\n", M_vir, c_vir, z)
    println("----------------------------------------")
    println("rho_s: $rho_s_kpc")
    println("r_s: $r_s")
    println("R_vir: $R_vir\n")

    return rho_s_kpc, r_s, R_vir
end

"""
    setup_dk_profile(M_vir, c_vir, nu; z=0.0, delta_vir=200.0, H0=67.74, Om0=0.3089)

Return a density function `r -> ρ(r)` ([`dk_profile`](@ref)) for a halo of
virial mass `M_vir` [M_sun] and concentration `c_vir`, with `r` in kpc and `ρ`
in M_sun/kpc^3. The virial radius encloses `delta_vir` times the mean matter
density at redshift `z`. Prints the derived `rho_s`, `r_s` and `R_vir`.
"""
function setup_dk_profile(M_vir, c_vir, nu; kwargs...)
    rho_s, r_s, R_vir = get_scale_rho_and_r_from_Mvir_cvir(M_vir, c_vir; kwargs...)
    x -> dk_profile(r_s, rho_s, R_vir, nu, x)
end

"""
    calc_psurf(rho_func, r_max, G)

Pressure at `r_max` from the weight of the material outside it, assuming
hydrostatic equilibrium: `P(r_max) = ∫_{r_max}^∞ ρ(r) G M(r) / r^2 dr`, with
`M(r)` the mass enclosed by `r`. Pass the result as `p_outer` to
[`evolve`](@ref). `rho_func(r)` must fall off fast enough for the integral to
converge.
"""
function calc_psurf(rho_func, r_max, G)
    M_max, _ = quadgk(r -> 4 * π * r^2 * rho_func(r), 0.0, r_max)
    P_integrand = function(r)
        m_tail, _ = quadgk(x -> 4 * π * x^2 * rho_func(x), r_max, r)
        M_enclosed = M_max + m_tail
        return rho_func(r) * G * M_enclosed / r^2
    end

    P_max, _ = quadgk(P_integrand, r_max, Inf)
    return P_max
end

"""
    initialize_profile(run_params, rho_initial; p_outer=1e-10, mass_min_frac=1e-7)

Integrate hydrostatic equilibrium inward in enclosed mass for the density
profile `rho_initial(r)`, from `run_params.r_max` (where the pressure is
`p_outer`) down to `mass_min_frac` times the total mass inside `r_max`.

Returns `(sol, s_func, M_tot)`: `sol(m)` gives `[r, P]` at enclosed mass `m`,
`s_func(m) = log(P^(3/2) / ρ^(5/2))`, and `M_tot` is the total mass.
[`evolve`](@ref) does this internally; call it directly only to inspect the
initial state.
"""
function initialize_profile(run_params, rho_initial; p_outer=1e-10, mass_min_frac=1e-7)
    # Integrate HSE inward in mass from r_max, with dr/dm = 1/(4 pi r^2 ρ)
    r_start = 1e-3 * run_params.r_min
    r_end = run_params.r_max

    # 1. Calculate total mass
    M_tot, _ = quadgk(r -> 4 * pi * r^2 * rho_initial(r), r_start, r_end)

    # 2. Inner mass boundary
    M_target = mass_min_frac * M_tot

    # 3. Set up and solve the ODE
    M_span = (M_tot, M_target)
    u0 = SVector(r_end, Float64(p_outer))
    G = run_params.G

    ode_params = (ρ_func = rho_initial, G = G)

    prob = ODEProblem(hse_system_rho, u0, M_span, ode_params)
    sol = solve(prob, reltol=1e-8)

    s_func(m) = log(sol(m)[2]^(3/2)/rho_initial(sol(m)[1])^(5/2))
    return sol, s_func, M_tot
end


function initialize_profile(run_params, rho_initial, M_min, M_tot; p_outer=1e-10, mass_min_frac=1e-7)
    # Same as above, with the density held flat inside the radius enclosing M_min
    r_start = 1e-3 * run_params.r_min
    r_end = run_params.r_max

    # 1. Radius enclosing M_min
    mass_diff(r) = quadgk(x -> 4 * pi * x^2 * rho_initial(x), r_start, r)[1] - M_min

    r_core = find_zero(mass_diff, (r_start, r_end))

    # 2. Create the modified, continuous density profile
    rho_core = rho_initial(r_core)
    rho_modified(r) = r > r_core ? rho_initial(r) : rho_core

    # 3. Set up and solve the ODE using the modified density
    M_span = (M_tot, M_min)
    u0 = SVector(r_end, Float64(p_outer))
    G = run_params.G

    ode_params = (ρ_func = rho_modified, G = G)

    prob = ODEProblem(hse_system_rho, u0, M_span, ode_params)
    sol = solve(prob, reltol=1e-8)

    s_func(m) = log(sol(m)[2]^(3/2)/rho_modified(sol(m)[1])^(5/2))
    return sol, s_func
end