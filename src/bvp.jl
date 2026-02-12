struct StructureParams{T}
    s::T
    G::Float64
end

Base.eltype(::Type{<:StructureParams}) = Float64
#hse equations log transformed
function hse_system!(du,u,p,logm)
        m = exp(logm)
        r = exp(u[1])
        P = exp(u[2])
        #u[3] = Phi
        ρ = P^(3/5) *exp(-2p.s(logm)/5) 
        du[1] = m/(4 * pi * r^3 * ρ)
        du[2] = - p.G * m^2 /(4 * pi * r^4 * P)
        du[3] = p.G * m^2/(4 * pi * r^4 * ρ)
    end

#hse equations used by init (not log transformed, for no good reason)
function hse_system_rho!(du,u,p,m)
    r = u[1]
    P = u[2]
    
    ρ = p.ρ_func(r)
    
    du[1] = 1/(4 * pi * r^2 * ρ)
    du[2] = - p.G * m /(4 * pi * r^4)
end

function shoot_1p(mspan, log_P_guess_in, log_r_min, log_p_surf, params)

    # Define Problem Template
    logM_inner, logM_outer = mspan[1], mspan[2]
    u0_in_template = [log_r_min, log_P_guess_in, 0.] 
    prob_inner_tmpl = ODEProblem{true,SciMLBase.FullSpecialize}(hse_system!, u0_in_template, (logM_inner, logM_outer), params)
    
    # Callback: Terminate if pressure hits surface target early
    condition(u, t, integrator) = u[2] < log_p_surf
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition, affect!)
    
    # settings
    tol = 1e-10
    max_iter = 1000

    # Returns negative if Pc is too low, positive if Pc is too high
    function get_residual(current_logP)
        prob_in = remake(prob_inner_tmpl, u0=[log_r_min, current_logP, 0.])
        sol_in = solve(prob_in, callback=cb,reltol=1e-8)
        #if Integration completed (Reached M_outer)
        if sol_in.t[end] >= logM_outer
            return (sol_in(logM_outer)[2] - log_p_surf)/(abs(log_p_surf)) # Positive Residual (Too much pressure left)
            
        else
            # Negative Residual (Ran out of pressure too early)
           
            return -1.0 * (logM_outer - sol_in.t[end])
        end
    end

    #  Find a Bracket
    step_size = 0.1 # Step size for log(P)
    p_low, p_high = log_P_guess_in, log_P_guess_in
    res_low, res_high = get_residual(p_low), get_residual(p_high)
    
    # Expand until we have one positive and one negative residual
    expand_iter = 0
    while sign(res_low) == sign(res_high) && expand_iter < 20
        p_high += step_size/2
        p_low -= step_size/2
        res_high = get_residual(p_high)
        res_low = get_residual(p_low)


        expand_iter += 1
    end
    
    if sign(res_low) == sign(res_high)
        error("Could not bracket the root.")
    end
    
    if res_low > 0
        p_low, p_high = p_high, p_low
        res_low, res_high = res_high, res_low
    end
    
    iter = 0
    p_mid = (p_high + p_low)/2
    while abs(0.5(p_high - p_low)/(p_high+p_low)) > tol && iter < max_iter
        iter += 1
        p_mid = (p_low + p_high) / 2
        res_mid = get_residual(p_mid)
        
        if res_mid == 0
            break # Exact match found
        elseif res_mid < 0
            p_low = p_mid # Root is in upper half
        else
            p_high = p_mid # Root is in lower half
        end
    end
    if iter == max_iter
        error("Shooting did not converge")
    end
    

    final_prob = remake(prob_inner_tmpl, u0=[log_r_min, p_mid, 0.])
    sol = solve(final_prob,callback=cb,reltol=1e-8)
    return sol
end