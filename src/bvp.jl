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

function shoot_1p(mspan, log_P_guess_in, log_r_min, log_lastm, log_p_surf, params)

    # Define Problem Template
    logM_inner, logM_outer = mspan[1], mspan[2]
    u0_in_template = [log_r_min, log_P_guess_in, 0.] 
    prob_inner_tmpl = ODEProblem{true,SciMLBase.FullSpecialize}(hse_system!, u0_in_template, (logM_inner, logM_outer), params)
    
    # Callback: Terminate if pressure hits surface target early
    condition(u, t, integrator) = u[2] < log_p_surf/5
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition, affect!)
    
    # settings
    tol = 1e-2
    max_iter = 50

    # Returns negative if Pc is too low, positive if Pc is too high
    function get_residual(current_logP)
        prob_in = remake(prob_inner_tmpl, u0=[log_r_min, current_logP, 0.])
        sol_in = solve(prob_in, callback=cb,reltol=1e-12)
        #if Integration reached last grid point
        if sol_in.t[end] >= log_lastm
            return (sol_in(log_lastm)[2] - log_p_surf) 
            
        else
            
            return -1.0 
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
    res_mid = Inf
    while abs(res_mid) > tol && iter < max_iter
        iter += 1
        p_mid = (p_low + p_high) / 2
        res_mid = get_residual(p_mid)

        if res_mid == 0 || isapprox(p_low, p_mid)
            break # Exact match found, or we ran out of floating point precision
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
    sol = solve(final_prob,callback=cb,reltol=1e-12)
    return sol
end

function shoot_2p_fit(mspan, log_P_guess_in, log_R_guess_out, log_r_min, log_p_surf, logM_fit, params)
    # 1. Setup the Problem Domains
    logM_inner = mspan[1]
    logM_outer = mspan[2]
    # Outward domain: Center to Fitting Mass
    span_out = (logM_inner, logM_fit)
    span_in  = (logM_outer, logM_fit)

    # 2. Define the Residual Function for NonlinearSolve
    function get_residuals!(resid, u, p)
        curr_log_P_c    = u[1]
        curr_log_R_surf = u[2]

        # --- OUTWARD INTEGRATION ---
        u0_out = [log_r_min, curr_log_P_c, 0.0] 
        prob_out = ODEProblem{true,SciMLBase.FullSpecialize}(hse_system!, u0_out, span_out, params)
        sol_out = solve(prob_out, reltol=1e-12)

        # --- INWARD INTEGRATION ---
        # Note: adjust the 3rd element of u0_in if your surface BC for that variable isn't 0.0
        u0_in = [curr_log_R_surf, log_p_surf, 0.0] 
        prob_in = ODEProblem{true,SciMLBase.FullSpecialize}(hse_system!, u0_in, span_in, params)
        sol_in = solve(prob_in, reltol=1e-12)

        
            
        resid[1] = sol_out.u[end][1] - sol_in.u[end][1]
        resid[2] =  sol_out.u[end][2] - sol_in.u[end][2]
    end

    # 3. Setup and Solve the Nonlinear Problem
    u0_guess = [log_P_guess_in, log_R_guess_out]
    nl_prob = NonlinearProblem(get_residuals!, u0_guess, nothing)
    
    # TrustRegion is much safer for stiff ODE boundary values than NewtonRaphson
    nl_sol = solve(nl_prob, RobustMultiNewton(), abstol=1e-2, maxiters=500)

    if !SciMLBase.successful_retcode(nl_sol.retcode)
        @warn "Shooting did not fully converge. Retcode: $(nl_sol.retcode). Final Residual norm: $(norm(nl_sol.resid))"
    end

    # 4. Generate the Final Converged Profiles
    best_log_P_c, best_log_R_surf = nl_sol.u


    final_prob_out = ODEProblem(hse_system!, [log_r_min, best_log_P_c, 0.0], span_out, params)
    final_sol_out = solve(final_prob_out, reltol=1e-12)    
    final_prob_in = ODEProblem(hse_system!, [best_log_R_surf, log_p_surf, 0.0], span_in, params)
    final_sol_in = solve(final_prob_in, reltol=1e-12)


    # Package and return
    return CompositeODESolution(final_sol_out, final_sol_in)
end