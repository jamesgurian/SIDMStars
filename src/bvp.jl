
function hse_system!(du,u,p,m)
    r = u[1]
    P = u[2]

    ρ = P^(3/5) *exp(-2p.s(m)/5) 
    du[1] = 1/(4 * pi * r^2 * ρ)
    du[2] = - p.G * m /(4 * pi * r^4)
end



function hse_system_rho!(du,u,p,m)
    r = u[1]
    P = u[2]
    
    ρ = p.ρ_func(r)

    du[1] = 1/(4 * pi * r^2 * ρ)
    du[2] = - p.G * m /(4 * pi * r^4)
end