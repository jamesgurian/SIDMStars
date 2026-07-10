K5(v2, run_params) = run_params.K5_prefactor * v2^(run_params.n / 2)
K7(v2, run_params) = run_params.K7_prefactor * v2^(run_params.n / 2)
K9(v2, run_params) = run_params.K9_prefactor * v2^(run_params.n / 2)

# Second-order Chapman-Enskog SMFP conductivity factor (Eq. k_secondorder):
#   K_eff^(2) = (28 K5^2 + 80 K5 K9 - 64 K7^2) / (77 K5 - 112 K7 + 80 K9)
function Keff2(v2, run_params)
    k5 = K5(v2, run_params); k7 = K7(v2, run_params); k9 = K9(v2, run_params)
    (28*k5^2 + 80*k5*k9 - 64*k7^2) / (77*k5 - 112*k7 + 80*k9)
end

# SMFP: κ_SMFP = (3/2)(b v / σ0) / K_eff^(2)   (Eq. kappa_smfp, b = 25√π/32),
# second Sonine approximation (Eq. k_secondorder).
function Kinv_smfp(v2, run_params)
    b = 25 * sqrt(pi)/32
    run_params.σ0 * Keff2(v2, run_params) / (b * sqrt(v2))
end


# LMFP: κ_LMFP uses K5 (Eq. kappa_lmfp), a = 4/√π
function Kinv_lmfp(v2, ρ, run_params,G)
    a = 4/sqrt(pi)
    4 * pi * G / ( a * run_params.C * ρ * sqrt(v2)^3 * run_params.σ0 * K5(v2,run_params))
end

# Total conductivity: SMFP and LMFP branches joined by the log-Gaussian
# resistivity bump measured from the DSMC heating trajectories (see common.jl),
#   1/κ = (Kinv_smfp + Kinv_lmfp) · [1 + A exp(−u²/2W²)],  u = ln(Kinv_smfp/Kinv_lmfp).
# Evaluated entirely in log space: both branch resistivities are power laws in
# (ρ, v²), so with the coefficients precomputed in RunParams
#   ln Kinv_smfp = lci_s + (n−1)/2 · ln v²
#   ln Kinv_lmfp = lci_l − ln ρ − (n+3)/2 · ln v²
# and u is their difference — no log or pow in the kernel, just three exp.
@inline function calc_κ_log(lρ, lv2, run_params)
    lki_s = run_params.lci_s + 0.5 * (run_params.n - 1) * lv2
    lki_l = run_params.lci_l - lρ - 0.5 * (run_params.n + 3) * lv2
    u = lki_s - lki_l
    return 1.0 / ((exp(lki_s) + exp(lki_l)) *
                  (1.0 + run_params.bump_A * exp(-u * u * run_params.inv2W2)))
end

# Value-space wrapper: ln ρ = (3/5)ln P − (2/5)s, ln v² = (2/5)(ln P + s).
# The hot loop in evolve.jl calls calc_κ_log directly with logs it already has.
@inline function calc_κ(P, s, run_params)
    lP = log(P)
    calc_κ_log(0.6 * lP - 0.4 * s, 0.4 * (lP + s), run_params)
end

struct EvaporationParams{T,U,V,W}
    logP::T 
    logr::U
    dPhi::V
    s::W
    run_params::RunParams
    G::Float64
    v2_i::Float64
end

#tau of surface velocity v 
#we integrate the optical depth, energy loss, and mass loss along particle trajectories
#this is a relatively heavy implementation, but it should be accurate even for low grid resolution 
function evaporation_integral!(du, u, p, logm)
    τ = u[1]
   #Edot = u[2]
    v2_part = p.v2_i + 2*p.dPhi(logm)
    P = exp(p.logP(logm))
    r = exp(p.logr(logm))
    m = exp(logm)
    ρ = P^(3/5) *exp(-2p.s(logm)/5)
    v2_local = P/ρ
    v2_rel = v2_part + 3*v2_local
    dr_dlogm = m/(4 * pi * r^2 * ρ)
    sigma = p.run_params.σ0 * (sqrt(v2_rel) / p.run_params.w)^p.run_params.n
    du[1] = -sigma * ρ * dr_dlogm
    shared_integrand = pi * r^2 * ρ *  du[1] *exp(-τ) 

    # \dot E(v) = \int  4 *(1/4) pi r^2 \rho v  (0.5 v^2) e^-(\tau) d\tau/dlogm dlogm/2 
    #flux is 1/4 n v 
    du[2] =    0.5 * v2_part^(3/2) * shared_integrand
    #\dot m = \int 0.5 pi r^2 \rho v e^(-\tau) d\tau/dlogm dlogm/2
    du[3] =  sqrt(v2_part) * shared_integrand 
end

condition_τ(u,t, integrator) =  integrator.u[1] > 10
            
cb_τ = DiscreteCallback(condition_τ, terminate!)

function calc_evaporation(run_params::RunParams{T, EvaporationOff},G, sol, s, m_face, m_center, v2s, v_vals) where T
    return zeros(length(m_center)), zeros(length(m_center))
end
function calc_evaporation(run_params::RunParams{T, EvaporationOn},G, sol, s, m_face, m_center, v2s, v_vals) where T
    #each column is one geodesic
    #with its own initial (surface) velocity
    udot_cell = zeros(length(m_face)-1,length(v_vals))
    mdot_cell = zeros(length(m_face)-1,length(v_vals))

    v_array = zeros(length(m_center),length(v_vals))

    Phi_surf =  sol(log(m_face[end]))[3] 
    dPhi(x) = Phi_surf - sol(x)[3] 
    log_m = log.(m_face)

    log_m_center = log.(m_center)
    
    u0 = [ 0., 0., 0.]
    m_span = (log_m[end], log_m[1])

    prob= ODEProblem{true, SciMLBase.FullSpecialize}(evaporation_integral!,u0, m_span)
    for i in eachindex(v_vals)
        pars = EvaporationParams(x->sol(x)[2], x->sol(x)[1], dPhi, s, run_params, G, v_vals[i]^2)
        prob= remake(prob, p=pars)
        sol_tau = solve(prob, callback = cb_τ, abstol = 1e-3, reltol = 1e-3, saveat=(log_m[1:end-1]))
        
      
        for j in 1:(length(m_face)-1)
            if log_m[j] < sol_tau.t[end] #if it's after we terminated (high tau)
                udot_cell[j,i] = 0.
                mdot_cell[j,i] = 0.
            else
                du_cell =(sol_tau(log_m[j+1])[2] - sol_tau(log_m[j])[2]) / (m_face[j+1] - m_face[j])
                udot_cell[j,i] = min(du_cell, 0.0) #guard against small positive du

                dm_cell = (sol_tau(log_m[j+1])[3] - sol_tau(log_m[j])[3])
                mdot_cell[j, i] = min(dm_cell, 0.0)
            end
        end
        v_array[:, i] .= sqrt.(v_vals[i]^2 .+ 2 .* dPhi.(log_m_center))
    end

 

    dv = diff(v_array, dims=2)
    v_integrand_u = maxwellian.(v2s, v_array) .* udot_cell
    v_integrand_u_centered = sqrt.(v_integrand_u[:, 1:end-1] .* v_integrand_u[:, 2:end])
    udot = -sum(v_integrand_u_centered .* dv, dims=2) #dv is negative, so we need a -

    v_integrand_m = maxwellian.(v2s, v_array) .* mdot_cell
    v_integrand_m_centered = sqrt.(v_integrand_m[:, 1:end-1] .* v_integrand_m[:, 2:end])
    mdot = -sum(v_integrand_m_centered .* dv, dims=2)
    #v2 = 2/3 u
    #s = (3/2)log(v2) - log(rho)
    sdot = udot ./ v2s .-  mdot ./ diff(m_face)
    #it's silly to diff then cumsum, when we have the integral directly but...
    mdot_cum = cumsum(mdot,dims=1)
    return sdot, mdot_cum

end
