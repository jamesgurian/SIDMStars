

function structure_system_T!(du,u,p,logm)
    #println(u)
    r = exp(u[1])
    P = exp(u[2])
    v2 = max(0,u[3])
    m = exp(logm)
       # println(m)

    L = u[4] * p.L_scale


    ρ = P/v2
    s = log(v2^(3/2)/ρ)
    du[1] = m/(4 * pi * r^3 * ρ)
    du[2] = - p.G * m^2 /(4 * P * pi * r^4)

    du[3] = -L*m/(16 * pi^2 * r^4 * ρ * p.κ(v2, ρ))
    println(m)
    println(p.s_prev(m))
    du[4] = -v2 * m *(s-p.s_prev(m))/(p.dt * p.L_scale)# p.s_prev(m)/p.L_scale

end

function bc_inner_general!(residual, u, p, r_inner, L_inner)
    #r = u[1]
    #L = u[4]
    residual[1] = u[1] - r_inner
    residual[2] = u[4] - L_inner

    #println("u, in", u)
   # println("residual", residual)

end


function bc_outer_general!(residual, u, p, P_surf, v2_surf)
    #P = u[2]
    #s = u[3]
    residual[1] = u[2] - P_surf
    residual[2] = u[3] - v2_surf
    println("u, out", u)
    println("residual", residual)

end