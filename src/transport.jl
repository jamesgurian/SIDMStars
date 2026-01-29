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
#tau of surface velocity v (v should be scalar or row vector)
function calc_τ(run_params, dr, Φ_vals, ρ_center,v2_center, v)
    v_vals = sqrt.(v.^2 .- 2 * (Φ_vals .- Φ_vals[end]))

    integrand = ρ_center .* run_params.σ0 .* (sqrt.(v_vals.^2 .+ v2_center.^2) ./ run_params.w).^(run_params.n) .* dr
        
    # Integrate along columns
    τ_vals = reverse(cumsum(reverse(integrand, dims=1), dims=1), dims=1)    
    return v_vals, τ_vals
end


function calc_τ(sol, s, m_face, m_center, run_params, v_surf)
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
        
        # --- Retrieve Derivatives (Geometry) ---
        # sol(logm, Val{1}) gives [d(logr)/d(logm), d(logP)/d(logm), dPhi/d(logm)]
        du_dlogm = sol(logm, Val{1})
        dlogr_dlogm = du_dlogm[1]
        
        dr_dlogm = r * dlogr_dlogm
        
        rho = P^(3/5) * exp(-2*s(logr)/5)
        
        v2_part = v_surf^2 .+ 2*(Phi_surf-Phi)
        
        # Thermal Dispersion 
        v2_disp = P / rho 
        
        # Relative Velocity
        v_rel = sqrt(v2_part + v2_disp)
        
        # Cross Section
        sigma = run_params.σ0 * (v_rel / run_params.w)^run_params.n
        
        # Return dτ/d(logm) = σ * ρ * (dr/d(logm))
        return sigma * rho * dr_dlogm
    end

    #Accumulate Integral over Grid
    n_points = length(logm_face_grid)
    v_vals = zeros(n_points)
    τ_vals = zeros(length(logm_center_grid)+1)
    
    v_vals[1] = v_surf + 2* (Phi_surf - sol(logm_center_grid[1])[3])
    τ_vals[1] = 0.0 
    Phi = Phi_surf
    for i in 1:n_points-1
        # Integrate from current logm to next smaller logm (moving inwards)
        dtau_signed, err = quadgk(integrand, logm_face_grid[i], logm_face_grid[i+1],rtol=1e-3 )
        
        # Accumulate
        # dr_dlogm is negative (radius shrinks as mass shrinks), so dtau_signed is negative.
        # Optical depth is a positive accumulation of path, so we subtract to flip the sign.
        τ_vals[i+1] = τ_vals[i] - dtau_signed
        
        # Calculate velocity at the specific grid point for output
        dPhi = sol(logm_center_grid[i])[3] - Phi
        Phi += dPhi
        v2_next = v_vals[i]^2 - 2 * dPhi
        v_vals[i+1] = v2_next > 0 ? sqrt(v2_next) : 0.0
    end
    τ_vals = τ_vals[2:end]
    
    return reverse(v_vals), reverse(τ_vals)
end
function calc_evaporation(run_params, Φ, dm, r, dr, ρ, v2s, v_vals )
    v_shifted, τ_vals = calc_τ(run_params, dr, Φ, ρ, v2s, v_vals')
    dvs = diff(v_shifted, dims=2)
    v_centered = 0.5(v_shifted[:,1:end-1] .+ v_shifted[:,2:end])
    τ_centered = 0.5(τ_vals[:,2:end] .+ τ_vals[:,1:end-1])
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
