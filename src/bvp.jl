
function hse_system!(du,u,p,logm)
        m = exp(logm)
        r = exp(u[1])
        P = exp(u[2])
        #u[3] = Phi
        ρ = P^(3/5) *exp(-2p.s(u[1])/5) 
        du[1] = m/(4 * pi * r^3 * ρ)
        du[2] = - p.G * m^2 /(4 * pi * r^4 * P)
        du[3] = p.G * m^2/(4 * pi * r^4 * ρ)
    end

function hse_system_rho!(du,u,p,m)
    r = u[1]
    P = u[2]
    
    ρ = p.ρ_func(r)
    
    du[1] = 1/(4 * pi * r^2 * ρ)
    du[2] = - p.G * m /(4 * pi * r^4)
end

function shoot_1p(mspan, log_P_guess_in, log_r_min, log_p_surf, params)
    # 1. Setup Grid
    logM_inner, logM_outer = mspan[1], mspan[2]

    # 2. Define Problem Templates
    u0_in_template = [log_r_min, log_P_guess_in, 0.] 
    prob_inner_tmpl = ODEProblem(hse_system!, u0_in_template, (logM_inner, logM_outer), params)
    
    # Callback: Terminate if Pressure hits surface target early
    condition(u, t, integrator) = u[2] < log_p_surf
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition, affect!)
    
    # settings
    tol = 1e-15
    max_iter = 1000

    # 3. Residual Function
    # Returns negative if Pc is too low, positive if Pc is too high
    function get_residual(current_logP)
        prob_in = remake(prob_inner_tmpl, u0=[log_r_min, current_logP, 0.])
        sol_in = solve(prob_in, callback=cb)
        # Case A: Integration completed (Reached M_outer)
        if sol_in.t[end] >= logM_outer || isapprox(sol_in.t[end], logM_outer, atol=1e-4)
            return (sol_in[end][2] - log_p_surf)/(abs(log_p_surf)) # Positive Residual (Too much pressure left)
            
        # Case B: Terminated early (Pressure dropped too fast)
        else
            # Negative Residual (Ran out of pressure too early)
           
            return -1.0 * (logM_outer - sol_in.t[end])
        end
    end

    # 4. Phase 1: Find a Bracket
    println("Searching for bracket...")
    step_size = 0.01 # Step size for log(P)
    p_low, p_high = log_P_guess_in, log_P_guess_in
    res_low, res_high = get_residual(p_low), get_residual(p_high)
    
    # Expand until we have one positive and one negative residual
    expand_iter = 0
    while sign(res_low) == sign(res_high) && expand_iter < 20
        if res_low < 0 
            # Both represent "too low", need to go higher
            p_high -= step_size
            res_high = get_residual(p_high)
            println(res_high)
        else 
            # Both represent "too high", need to go lower
            p_low += step_size
            res_low = get_residual(p_low)
        end
        expand_iter += 1
    end
    
    if sign(res_low) == sign(res_high)
        error("Could not bracket the root. Try a better initial guess.")
    end
    
    # Ensure p_low corresponds to negative residual for standard logic
    if res_low > 0
        p_low, p_high = p_high, p_low
        res_low, res_high = res_high, res_low
    end
    
    println("Bracket found: [$p_low, $p_high]")

    # 5. Phase 2: Bisection Method
    println("Starting Bisection...")
    iter = 0
    p_mid = (p_high + p_low)/2
    #println("0.5(p_high - p_low)/(p_high+p_low) ", 0.5(p_high - p_low)/(p_high+p_low) )
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
    
    println("Terminated after $iter iterations to Log(Pc) = $p_mid")

    # 6. Return Final Solution
    final_prob = remake(prob_inner_tmpl, u0=[log_r_min, p_mid, 0.])
    sol = solve(final_prob,callback=cb)
    return sol
end
function shoot(mspan, log_P_guess_in, log_R_guess_out, log_r_min, log_p_surf, params)
    # 1. Setup Grid
    logM_inner, logM_outer = mspan[1], mspan[2]
    logM_fit = 0.5 * (logM_inner + logM_outer) 
    
    # State Vector: x = [log_P_inner, log_R_outer]
    x = [log_P_guess_in, log_R_guess_out]
 #   println("Initial Guess: $x")
    println("Mfit", logM_fit)
    # 2. Define Problem Templates
    u0_in_template = [log_r_min, x[1]] 
    prob_inner_tmpl = ODEProblem(hse_system!, u0_in_template, (logM_inner, logM_fit), params)
    
    u0_out_template = [x[2], log_p_surf]
    prob_outer_tmpl = ODEProblem(hse_system!, u0_out_template, (logM_outer, logM_fit), params) 

    # Newton Parameters
    tol = 1e-4
    max_iter = 50
    epsilon = 1e-5

    condition_out(u, t, integrator) = u[2] > log_P_guess_in
    affect!(integrator) = terminate!(integrator)
    cb_out = DiscreteCallback(condition_out, affect!)
    
    # --- HELPER: Attempt a shot with specific parameters ---
    # Returns: (success::Bool, F::Vector, sol_in, sol_out)
    function try_shot(current_x)
        # Inner Shot
        prob_in = remake(prob_inner_tmpl, u0=[log_r_min, current_x[1]])
        #println("in")

        sol_in  = solve(prob_in, reltol=1e-3, abstol=1e-3,alg_hints=:stiff)
       # println("uend", sol_in.u[end])
        # Check Inner Failure
        

        # Outer Shot
        
        prob_out = remake(prob_outer_tmpl, u0=[current_x[2], log_p_surf])
        sol_out  = solve(prob_out, reltol=1e-4, abstol=1e-3,callback=cb_out)
        #println("uendout", sol_out.u[end])

        # Check Outer Failure
        

        # Calculate Residuals at fitting point
        end_in  = sol_in[end]
        end_out = sol_out[end]
        
        # F = Inner - Outer
        F_val = [(end_in[1] - end_out[1])/end_in[1], 
                 (end_in[2] - end_out[2])/end_out[2]]
        if sol_out.retcode != :Success || sol_in.retcode != :Success
            return (false, F_val, sol_in, sol_out)
        end
        return (true, F_val, sol_in, sol_out)
    end

    # --- INITIAL SHOT ---
    
    success, F, sol_in, sol_out = try_shot(x)
    while !success
        #println("uh oh", sol_out.u)
        #println("P_gu", log_P_guess_in)
        if sol_out.u[end][2] > log_P_guess_in
            x[2] *= 1.001
            success, F, sol_in, sol_out = try_shot(x)
        
        end

    end
    iter = 0
    println("Starting Two-Point Shooting ...")

    while iter < max_iter
        iter += 1
       # println("\n--- Iteration $iter, Error: $err ---")
        err = norm(F)
        # --- CONVERGENCE CHECK ---
        if err < tol
            println("Converged on iter ", iter)
            # Stitch solution for full domain return
            final_prob = remake(prob_inner_tmpl, 
                                u0=[log_r_min, x[1]], 
                                tspan=(logM_inner, logM_outer))
            return solve(final_prob, reltol=1e-4, abstol=1e-4, alg_hints=:stiff)
        end
        
        # --- JACOBIAN (Finite Difference) ---
        J = zeros(2,2)
        
        # Col 1: Perturb P_in
        s_p1, F_p1, _, _ = try_shot([x[1] + epsilon, x[2]])
       
        J[:,1] = (F_p1 - F) / epsilon
      

        # Col 2: Perturb R_out
        s_p2, F_p2, _, _ = try_shot([x[1], x[2] + epsilon])
        J[:,2] = (F_p2 - F) / epsilon

        # Calculate Full Newton Step Direction
        # We protect against Singular Matrix errors here
        delta_x = J\F

      #  println("Full Newton Step: $delta_x")

        # --- BACKTRACKING LINE SEARCH ---
        lambda = 0.5
        step_accepted = false
        min_lambda = 0.01

        while lambda >= min_lambda
            x_test = x .- (lambda .* delta_x)
            
            # Try the shot with reduced step
            s_test, F_test, _, _ = try_shot(x_test)
            
            # Criteria: 
            # 1. Integration must succeed (:Success)
            # 2. Residual norm should ideally decrease (Damping)
            if s_test && (norm(F_test) < err || lambda < 0.1)
           #     println("  Step accepted with lambda = $lambda. New Error: $(norm(F_test))")
                x = x_test
                F = F_test
                step_accepted = true
                break
            else
             #   println("  Failed with lambda = $lambda (Retcode Success: $s_test, Err: $(s_test ? norm(F_test) : "N/A"))")
                lambda *= 0.5
            end
        end

        if !step_accepted
          #  error("Line search failed: Step size became too small without converging or fixing integrator crash.")
        end
    end
    
    error("Shooting failed to converge after $max_iter iterations")
end