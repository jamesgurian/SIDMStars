#f_k5(v_over_w,n) = 1/3 * 2^(n-1) * v_over_w^n * gamma(n/2 + 4)
f_k5(v_over_w,n) = (1/6) * (sqrt(2)*v_over_w)^n * gamma(n/2 + 4) 
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
    v2_rel = v2_part + v2_local
    dr_dlogm = m/(4 * pi * r^2 * ρ)
    sigma = p.run_params.σ0 * (sqrt(v2_rel) / p.run_params.w)^p.run_params.n
    du[1] = -sigma * ρ * dr_dlogm
    shared_integrand = pi * r^2 * ρ *  du[1] *exp(-τ) 

    # \dot E(v) = \int  4 pi r^2 \rho v  (0.5 v^2) e^-(\tau) d\tau/dlogm dlogm
    du[2] =    2 * v2_part^(3/2) * shared_integrand
    #\dot m = \int 4 pi r^2 \rho v e^(-\tau) d\tau/dlogm dlogm
    du[3] =  4 * sqrt(v2_part) * shared_integrand 
end

condition_τ(u,t, integrator) =  integrator.u[1] > 10
            
cb_τ = DiscreteCallback(condition_τ, terminate!)
function calc_evaporation(run_params,G, sol, s, m_face, m_center, v2s, v_vals)
    #each column is one geodesic
    #with its own initial (surface) velocity
    udot_cell = zeros(length(m_face)-1,length(v_vals))
    mdot_cell = zeros(length(m_face)-1,length(v_vals))

    v_array = zeros(length(m_center),length(v_vals))

    #v_array = zeros(length(m_center), length(v_vals ))
    Phi_surf =  sol(log(m_face[end]))[3] 
    dPhi(x) = Phi_surf - sol(x)[3] 
    log_m = log.(m_face)

    log_m_center = log.(m_center)
    
    u0 = [ 0., 0., 0.]
    m_span = (log_m[end], log_m[1])

    prob= ODEProblem(evaporation_integral!,u0, m_span)

    for i in eachindex(v_vals)

        pars = EvaporationParams(x->sol(x)[2], x->sol(x)[1], dPhi, s, run_params, G, v_vals[i]^2)
        prob= remake(prob, p=pars)
        sol_tau = solve(prob, callback = cb_τ, abstol = 1e-3, reltol = 1e-3)
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
