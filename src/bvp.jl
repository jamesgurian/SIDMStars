function structure_system!(du,u,p,t)
    r = u[1]
    P = u[2]
    s = u[3]
    L = u[4]
    m = t

    ρ = P^(3/5) *exp(-2s/5)
    du[1] = 1/(4 * pi * r^2 * ρ)
    du[2] = - G * m /(4 * pi * r^4)
    du[3] = - 5/2 * L/(p.κ(P/ρ) * P * 16 * pi^2 * r^4) - du[2]/P 
    dsdt = (s - p.s_prev(m))/p.dt
    du[4] = -P/ρ * dsdt

end


