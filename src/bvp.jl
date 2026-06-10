struct StructureParams{T}
    s::T
    G::Float64
end

Base.eltype(::Type{<:StructureParams}) = Float64
#hse equations log transformed
function hse_system!(du, u, p, logm)
    m = exp(logm)
    r = exp(u[1])
    P = exp(u[2])
    
    # 1. Eliminate the slow fractional power and combine exponentials
    ρ = exp(0.6 * u[2] - 0.4 * p.s(logm))
    
    # 2. Precompute repeated terms
    four_pi_r3 = 4 * pi * r^3
    four_pi_r4 = four_pi_r3 * r
    Gm2 = p.G * m^2
    
    # 3. Update the preallocated array
    du[1] = m / (four_pi_r3 * ρ)
    du[2] = -Gm2 / (four_pi_r4 * P)
    du[3] = Gm2 / (four_pi_r4 * ρ)
    
    return nothing
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
    condition(u, t, integrator) = u[2] < log_p_surf/5
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition, affect!)
    
    # settings
    tol = 1e-2
    max_iter = 100

    # Returns negative if Pc is too low, positive if Pc is too high
    function get_residual(current_logP)
        prob_in = remake(prob_inner_tmpl, u0=[log_r_min, current_logP, 0.])
        sol_in = solve(prob_in, callback=cb,reltol=1e-8)
        #if Integration reached last grid point
        if sol_in.t[end] >= mspan[2]
            return (sol_in(logM_outer)[2] - log_p_surf) 
            
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
    sol = solve(final_prob,callback=cb,reltol=1e-8)
    return sol
end
# Fully parameterized, type-stable Exception
# Back to a strictly typed, lightweight struct
struct IntegrationFailure <: Exception
    t_out_end::Float64
    t_in_end::Float64
    u::Vector{Float64}
end
function shoot_2p_fit(mspan, log_P_guess_in, log_R_guess_out, log_r_min, log_p_surf, logM_fit, params, abstol, reltol)
    logM_inner = mspan[1]
    logM_outer = mspan[2]
    
    span_out = (logM_inner, logM_fit)
    span_in  = (logM_outer, logM_fit)
    
    prob_out = ODEProblem{true,SciMLBase.FullSpecialize}(hse_system!, (0., 0., 0.), span_out, params)
    prob_in  = ODEProblem{true,SciMLBase.FullSpecialize}(hse_system!, (0., 0., 0.), span_in, params)

    # Define the Residual Function for NonlinearSolve
    function get_residuals!(resid, u, p)
        curr_log_P_c    = u[1]
        curr_log_R_surf = u[2]

        # --- OUTWARD INTEGRATION ---
        u0_out = [log_r_min, curr_log_P_c, 0.0] 
        p_out = remake(prob_out, u0=u0_out, tspan=span_out)
        sol_out = solve(p_out, reltol=reltol, abstol=abstol, save_everystep=false)

        # --- INWARD INTEGRATION ---
        u0_in = [curr_log_R_surf, log_p_surf, 0.0] 
        p_in = remake(prob_in, u0=u0_in, tspan=span_in)
        sol_in = solve(p_in, reltol=reltol, abstol=abstol, save_everystep=false)

        if !SciMLBase.successful_retcode(sol_out.retcode) || !SciMLBase.successful_retcode(sol_in.retcode)
            # 2. Strip Dual tags safely at the throw site to maintain global type stability
            t_out_fail = ForwardDiff.value(sol_out.t[end])
            t_in_fail  = ForwardDiff.value(sol_in.t[end])
            u_fail     = ForwardDiff.value.(u)
            
            throw(IntegrationFailure(t_out_fail, t_in_fail, u_fail))
        end
            
        resid[1] = sol_out.u[end][1] - sol_in.u[end][1]
        resid[2] = sol_out.u[end][2] - sol_in.u[end][2]
    end

    nf = NonlinearFunction{true, SciMLBase.FullSpecialize}(get_residuals!)
    u0_guess = [log_P_guess_in, log_R_guess_out]
    local nl_sol = nothing
    
    max_retries = 50 
    retries = 0
    
    # 3. State flag to track our fallback logic
    moved_fitting_point = false

    while abstol >= 1e-13 && retries < max_retries
        nl_prob = NonlinearProblem(nf, u0_guess, nothing)
        
        try
            nl_sol = solve(nl_prob, RobustMultiNewton(), abstol=1e-2, maxiters=1000)
            
            # If we get here without throwing, integration succeeded for all steps.
            if SciMLBase.successful_retcode(nl_sol.retcode)
                break 
            else
                println("NLsolve failed to converge. Reducing ODE solve tolerance.")
                println("abstol = $abstol")
                abstol /= 10.
                reltol /= 10.
                u0_guess = nl_sol.u
                moved_fitting_point = false # Reset flag
            end
            
        catch e
            if e isa IntegrationFailure
                if !moved_fitting_point
                    println("Integrations did not meet, adjusting fitting point...")
                    println("Outward terminated at ", exp(e.t_out_end))
                    println("Inward terminated at ", exp(e.t_in_end))

                    if e.t_in_end != logM_fit && e.t_out_end == logM_fit
                        println("Inward integration failed")
                        logM_fit = e.t_in_end
                    elseif e.t_out_end != logM_fit && e.t_in_end == logM_fit
                        println("Outward integration failed")
                        logM_fit = e.t_out_end
                    else
                        println("Both integrations failed. Good luck!")
                        logM_fit = e.t_in_end 
                    end
                    
                    println("new M_fit = $(exp(logM_fit))")
                    
                    # Update integration spans for the next attempt
                    span_out = (logM_inner, logM_fit)
                    span_in  = (logM_outer, logM_fit)
                    
                    u0_guess = e.u 
                    moved_fitting_point = true # Mark that we tried moving it
                else
                    # It failed again, even after we just moved the fitting point
                    println("Integration failed again after moving fitting point. Reducing ODE solve tolerance.")
                    println("abstol = $abstol")
                    abstol /= 10.
                    reltol /= 10.
                    u0_guess = e.u 
                    moved_fitting_point = false # Reset flag so we are allowed to move it again on the next failure
                end
                retries += 1
            else
                rethrow(e) # Let things like DomainErrors or Interrupts bubble up
            end
        end
    end

    if isnothing(nl_sol) || !SciMLBase.successful_retcode(nl_sol.retcode)
        ret_msg = isnothing(nl_sol) ? "Max retries reached" : string(nl_sol.retcode)
        resid_msg = isnothing(nl_sol) ? "N/A" : string(norm(nl_sol.resid))
        @warn "Shooting did not converge. Retcode: $ret_msg. Final Residual norm: $resid_msg"
    end

    # Generate the Final Converged Profiles
    best_log_P_c, best_log_R_surf = isnothing(nl_sol) ? u0_guess : nl_sol.u

    p_out_final = remake(prob_out, u0=[log_r_min, best_log_P_c, 0.0], tspan=span_out)
    final_sol_out = solve(p_out_final, reltol=reltol, abstol=abstol)    
    
    p_in_final = remake(prob_in, u0=[best_log_R_surf, log_p_surf, 0.0], tspan=span_in)
    final_sol_in = solve(p_in_final, reltol=reltol, abstol=abstol)
    
    # Package and return
    return CompositeODESolution(final_sol_out, final_sol_in), abstol, reltol, logM_fit
end