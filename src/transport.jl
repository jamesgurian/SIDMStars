f_k5(v_over_w,n) = 1/3 * 2^(n-1) * v_over_w^n * gamma(n/2 + 4)

function Kinv_smfp(v2, dm_params) 
    b = 25 * sqrt(pi)/32 #smfp coefficient

    4 * pi * dm_params.σ0 * f_k5(sqrt(v2)/dm_params.w, dm_params.n)/(b * sqrt(v2)) 
end

function Kinv_lmfp(v2, ρ, dm_params,G)
    a = 4/sqrt(pi)

    4 * pi * G / ( a * dm_params.C * ρ * sqrt(v2)^3 * dm_params.σ0 * f_k5(sqrt(v2)/dm_params.w,dm_params.n))
end

function K_eff(v2,ρ, dm_params)
    1/(Kinv_smfp(v2, dm_params) + Kinv_lmfp(v2, ρ, dm_params, G))
end

struct EvaporationParams
    # Interpolations for background structure: Input is logm, Output is Float64
    logP::FunctionWrapper{Float64, Tuple{Float64}} 
    logr::FunctionWrapper{Float64, Tuple{Float64}} 
    dPhi::FunctionWrapper{Float64, Tuple{Float64}} 
    # Entropy profile: Input is logr, Output is Float64
    s::FunctionWrapper{Float64, Tuple{Float64}}    

    # Nested parameters and constants
    run_params::RunParams
    G::Float64
    v2_i::Float64
end

#tau of surface velocity v (v should be scalar or row vector)
function evaporation_integral!(du, u, p, logm)
    τ = u[1]
   #Edot = u[2]
    v2_part = p.v2_i + 2*p.dPhi(logm)
    P = exp(p.logP(logm))
    r = exp(p.logr(logm))
    m = exp(logm)
    ρ = P^(3/5) *exp(-2p.s(p.logr(logm))/5)
    v2_local = P/ρ
    v2_rel = v2_part + v2_local
    dr_dlogm = m/(4 * pi * r^2 * ρ)
    sigma = p.run_params.σ0 * (sqrt(v2_rel) / p.run_params.w)^p.run_params.n

    du[1] = -sigma * ρ * dr_dlogm
    # \dot E(v) = \int  4 pi r^2 \rho v  (0.5 v^2) e^-(\tau) d\tau/dlogm dlogm
    du[2] = 2 * pi * r^2 * ρ * v2_part^(3/2) * exp(-τ) * du[1]
end

condition_τ(u,t, integrator) =  integrator.u[1] > 10
            
cb_τ = DiscreteCallback(condition_τ, terminate!)
function calc_evaporation(run_params,G, sol, s, m_face, m_center, v2s, v_vals)
    #each column is one geodesic
    #with its own initial (surface) velocity
    udot_cell = zeros(length(m_face)-1,length(v_vals))
    v_array = zeros(length(m_center),length(v_vals))

    #v_array = zeros(length(m_center), length(v_vals ))
    Phi_surf =  sol(log(m_face[end]))[3] 
    dPhi(x) = Phi_surf - sol(x)[3] 
   
    log_m = log.(m_face)

    log_m_center = log.(m_center)
    figure()
    
    u0 = [ 0., 0.]
    m_span = (log_m[end], log_m[1])

    prob= ODEProblem(evaporation_integral!,u0, m_span)

    for (i, v_val) in enumerate((v_vals))
        pars = EvaporationParams(x->sol(x)[2], x->sol(x)[1], dPhi, s, run_params, G, v_val^2)

        
        prob= remake(prob, p=pars)#ODEProblem(evaporation_integral!,u0, m_span, pars)
        sol_tau = solve(prob, callback = cb_τ)
        for j in 1:(length(m_face)-1)
            if log(m_face[j+1]) < sol_tau.t[end] #if it's after we terminated (high tau)
                udot_cell[j,i] = 0.
            else
                udot_cell[j,i] = sol_tau(log(m_face[j+1]))[2] - sol_tau(log(m_face[j]))[2]
            end
        end
        v_array[:, i] .= sqrt.(v_val^2 .+ 2 .* dPhi.(log_m_center))
       # loglog(exp.(sol(sol_tau.t)[1,:]), sol_tau[1,:])
    end
    udot_cell =  udot_cell ./diff(m_face)
    dv = diff(v_array, dims=2)
    v_integrand = maxwellian.(v2s, v_array) .* udot_cell
   # loglog(m_center, abs.(udot_cell))
  #  loglog(exp.(sol(log_m_center)[1,:]), v_array)
    v_integrand_centered = sqrt.(v_integrand[:, 1:end-1] .* v_integrand[:, 2:end])
    udot = -sum(v_integrand_centered .* dv, dims=2)
    #v2 = 2/3 u
    #s = (3/2)log(v2) - log(rho)
    
    sdot = udot ./ v2s 


end

function calc_τ(run_params, dr::Vector{Float64}, Φ_vals, ρ_center,v2_center, v)
    v_vals = sqrt.(v.^2 .- 2 * (Φ_vals .- Φ_vals[end]))

    integrand = ρ_center .* run_params.σ0 .* (sqrt.(v_vals.^2 .+ v2_center) ./ run_params.w).^(run_params.n) .* dr
        
    # Integrate along columns
    τ_vals = reverse(cumsum(reverse(integrand, dims=1), dims=1), dims=1)    
    return v_vals, τ_vals
end

#=
function calc_τ(run_params,sol, s, m_face, m_center,  v_surf)
    #integrate from surface
    logm_face_grid = reverse(log.(m_face))
    logm_center_grid = reverse(log.(m_center))
    #println("m", logm_grid)
    Phi_surf = sol(logm_face_grid[1])[3]
    
   #semilogx(logm_grid, sol(logm_grid)[3,:])
    # dτ/d(logm) = σ * ρ * (dr/d(logm))
    function integrand(logm)
        # --- Retrieve State ---
        u = sol(logm)
        logr, logP, Phi = u
        r = exp(logr)
        P = exp(logP)
        
        rho = P^(3/5) * exp(-2*s(logr)/5)
        
        v2_part = v_surf^2 .+ 2*(Phi_surf-Phi)
        
        # Thermal Dispersion 
        v2_disp = P / rho 
        
        # Relative Velocity
        v_rel = sqrt(v2_part + v2_disp)
        
        # Cross Section
        sigma = run_params.σ0 * (v_rel / run_params.w)^run_params.n
        dr_dlogm = exp(logm)/(4 * pi * r^2 * rho)
        
        return sigma * rho * dr_dlogm
    end

    #Accumulate Integral over Grid
    n_points = length(logm_face_grid)
    v_vals = zeros(length(logm_center_grid))
    τ_vals = zeros(length(logm_face_grid))
    
    v_vals[1] = sqrt(v_surf^2 + 2 * (Phi_surf - sol(logm_center_grid[1])[3]))
    τ_vals[1] = 0.0 
    
    Phi = Phi_surf
    for i in 1:(length(v_vals)-1)
         # Calculate velocity at the specific grid point for output
         dPhi = sol(logm_center_grid[i])[3] - Phi
         Phi += dPhi
         v2_next = v_vals[i]^2 - 2 * dPhi
         v_vals[i+1] = sqrt(v2_next)
    end
    for i in 1:n_points-1
        # Integrate from current logm to next smaller logm (moving inwards)
        dtau_signed, _ = quadgk(integrand, logm_face_grid[i], logm_face_grid[i+1],rtol=1e-3 )
        
        # Accumulate
        # dr_dlogm is negative (radius shrinks as mass shrinks), so dtau_signed is negative.
        # Optical depth is a positive accumulation of path, so we subtract to flip the sign.
        τ_vals[i+1] = τ_vals[i] - dtau_signed
        
       
    end
  #  loglog(m_face, reverse(τ_vals))
    
    return reverse(v_vals), reverse(τ_vals)
end

function calc_evaporation(run_params, sol, s, m_face, m_center, dm_face, ρ,r_center, v2s, v_vals)
    #each column is one geodesic
    #with its own initial (surface) velocity
    τ_array = zeros(length(m_face),length(v_vals))
    v_array = zeros(length(m_center), length(v_vals ))
    for (i, v_val) in enumerate(v_vals)
        vs, τs = calc_τ(run_params, sol, s, m_face, m_center, v_val)
        
        v_array[:, i] .= vs
        τ_array[:, i] .= τs
    end
    τ_centered_v = sqrt.(τ_array[:, 1:end-1] .* τ_array[:,2:end]) #we need to reduce the dimensionality along the v dimension, too

    dτ = diff(τ_centered_v, dims=1) #now the accumulation optical depth along each step on the geodesics
    τ_centered =sqrt.(τ_centered_v[2:end, :] .* τ_centered_v[1:end-1,:]) #and we need the value of \tau at the spatial grid centers
    v_centered = sqrt.(v_array[:, 2:end] .* v_array[:, 1:end-1]) #we need to integrate dv. 
    dvs = diff(v_array, dims=2)
    # \dot \u = \int dv 4 pi r^2 \rho/m v f(v) (0.5 v^2) e^-(\tau) d\tau
    #v2 = 2/3 u
    #TODO: there should be some geometric factor to only have outgoing v? Is it just 0.5 in the radial orbit approximation? 
    v2dot_vintegrated = -sum((@. 4. / 3.  * pi * r_center^2 * ρ/dm_face * v_centered^3 * maxwellian(v2s, v_centered) * dvs* exp(-τ_centered) * dτ),dims=2)
    #s = log(v2^(3/2)) - log(rho)
    #dsdt = 3/2 d log(v^2)/dt = 3/2 dv^2/dt /v2
    sdot = 3/2 * v2dot_vintegrated ./ v2s 
end=#
function calc_evaporation(run_params, Φ::Vector{Float64}, dm, r, dr, ρ, v2s, v_vals )
    v_shifted, τ_vals = calc_τ(run_params, dr, Φ, ρ, v2s, v_vals')
    dvs = diff(v_shifted, dims=2)
    v_centered = sqrt.(v_shifted[:,1:end-1] .* v_shifted[:,2:end])
    τ_centered = sqrt.(τ_vals[:,2:end] .* τ_vals[:,1:end-1])
  #  loglog(r, τ_vals,"--")
    #loglog(r, v_shifted, "--")
    dτ = -diff(vcat(τ_centered, zeros(1,length(v_vals)-1)), dims=1)
    # \dot \u = \int dv 4 pi r^2 \rho/m v f(v) (0.5 v^2) e^-(\tau) d\tau
    #v2 = 2/3 u
    #TODO: there should be some geometric factor to only have outgoing v? Is it just 0.5 in the radial orbit approximation? 
    
    v2dot_vintegrated = -sum((@. 4. / 3.  * pi * r^2 * ρ/dm * v_centered^3 * maxwellian(v2s, v_centered)) .* dvs.* exp.(-τ_centered) .* dτ, dims=2)
    #s = log(v2^(3/2)) - log(rho)
    #dsdt = 3/2 d log(v^2)/dt = 3/2 dv^2/dt /v2
    sdot = 3/2 * v2dot_vintegrated ./ v2s 
    #loglog(r, -v2dot_vintegrated./v2s, "--")
    return sdot
end
