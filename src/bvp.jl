
function hse_system!(du,u,p,logm)
        m = exp(logm)
        r = exp(u[1])
        P = exp(u[2])
        ρ = P^(3/5) *exp(-2p.s(m)/5) 
        du[1] = m/(4 * pi * r^3 * ρ)
        du[2] = - p.G * m^2 /(4 * pi * r^4 * P)
    end

function hse_system_rho!(du,u,p,m)
    r = u[1]
    P = u[2]
    
    ρ = p.ρ_func(r)

    du[1] = 1/(4 * pi * r^2 * ρ)
    du[2] = - p.G * m /(4 * pi * r^4)
end

function shoot(mspan, guess_inner, guess_outer, r_min, p_surf,params)
    fitting_point = 0.5*(mspan[2] - mspan[1])
    prob_inner = ODEProblem(hse_system!, guess_inner, (m_span[1],fitting_point),params)
    prob_outer = ODEProblem(hse_system!, guess_outer, (m_span[1],fitting_point),params)
    sol_inner = solve(prob_inner)
    sol_outer = solve(prob_outer)
    while #condition
        #adjust guess
    end


end
