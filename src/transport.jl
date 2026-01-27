f_k5(v2_over_w,n) = 1/6 * sqrt(2 * v2_over_w)^n * gamma(n/2 + 4)

function Kinv_smfp(v2, dm_params) 
    b = 25 * sqrt(pi)/32 #smfp coefficient


    dm_params.σ0 * f_k5(sqrt(v2)/dm_params.w, dm_params.n)/(b * sqrt(v2)) 
end

function Kinv_lmfp(v2, ρ, dm_params,G)
    a = 4/sqrt(pi)
   
    4 * pi * G / ( a * dm_params.C * ρ * sqrt(v2)^3 * dm_params.σ0 * f_k5(sqrt(v2)/dm_params.w,dm_params.n))
end

function K_eff(v2,ρ, dm_params)
    1/(Kinv_smfp(v2, dm_params) + Kinv_lmfp(v2, ρ, dm_params, G))
end

#tau of surface velocity v (v should be scalar or row vector)
function calc_τ(run_params, dr, Φ_vals, ρ_center,v)
    v_vals = sqrt.(v.^2 .- 2 * (Φ_vals .- Φ_vals[end]))
    integrand = ρ_center .* run_params.σ0 .* (v_vals ./ run_params.w).^(run_params.n) .* dr
        
    # Integrate along columns
    #this allocates a bunch, so we should preallocate if this is ever expensive
    τ_vals = reverse(cumsum(reverse(integrand, dims=1), dims=1), dims=1)    
    return v_vals, τ_vals
end

function calc_evaporation(run_params, Φ, r, dr, ρ, v2s, v_vals )
    v_shifted, τ_vals = calc_τ(run_params, dr, Φ, ρ, v_vals')
    dvs = diff(v_shifted, dims=2)
    v_centered = 0.5(v_shifted[:,1:end-1] .+ v_shifted[:,2:end])
        # 4 pi r^2 \rho v f(v) (0.5 v^2) e^-(\tau) d\tau
    τ_centered = 0.5(τ_vals[:,2:end] .+ τ_vals[:,1:end-1])
    dτ = -diff(vcat(τ_centered, zeros(1,length(v_vals)-1)), dims=1)
    udot_vintegrated = sum((@. 2 * pi * r^2 * ρ * v_centered^3 * maxwellian(v2s, v_centered)) .* dvs.* exp.(-τ_centered) .* dτ, dims=2)
   #loglog(r, τ_vals)
    #figure()
    loglog(r,udot_vintegrated )
    return (udot_vintegrated )
#    println("size",sum((@. 2 * pi * r^2 * ρ * v_shifted^3 * maxwellian(v2s, v_shifted).*exp(-τ_vals) .* dτ))
end