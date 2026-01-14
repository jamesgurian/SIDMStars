function structure_system!(du,u,p,m)
    r = u[1]
    P = (u[2]) #solve for log P
    s = u[3]
    L = u[4]
    println()
    println("m", m)
    println(u)
    #println("u_in", u)
    ρ = P^(3/5) *exp(-2s/5)
    println("rho", ρ)
    #println("ρ", ρ)
    if r < 1e-5
        du[1] = 0
        du[2] = 0
        du[3] = 0
        du[4] = 0 
    else
    du[1] = 1/(4 * pi * r^2 * ρ)
    du[2] = - p.G * m /(4 * P * pi * r^4)
    du[3] = - 5/2 * L/(p.κ(P/ρ,ρ) * 16 * pi^2 * r^4* ρ) - du[2]/P 
    #du[4] = -P/ρ * (u[3] - p.s_prev(m))/p.dt
    du[4] = - P/ρ * p.dsdt(m)
        println(du)
    end
end

function structure_system_T!(du,u,p,m)
    r = u[1]
    P = max(0, u[2]) #solve for log P
    v2 = max(0., u[3]) #negative will throw error in sqrt
    L = u[4]
   # println("u_in", u)
    ρ = P/v2
    s = log(v2^(3/2)/ρ)
    du[1] = 1/(4 * pi * r^2 * ρ)
    du[2] = - p.G * m /(4 * pi * r^4)
    #du[3] = - 5/2 * L/(p.κ(P/ρ,ρ) * P * 16 * pi^2 * r^4) - du[2]/P 
    #du[4] = -P/ρ * (u[3] - p.s_prev(m))/p.dt
    du[3] = -L/(16 * pi^2 * r^4 * ρ * p.κ(v2, ρ))
    du[4] = - v2 * p.dsdt(m)
    #println("du_out", du)
end

function bc_inner_general!(residual, u, p, r_inner, L_inner)
    #r = u[1]
    #L = u[4]
    residual[1] = u[1] - r_inner
    residual[2] = u[4] - L_inner

    println("u, in", u)
    println("residual", residual)

end


function bc_outer_general!(residual, u, p, P_surf, v2_surf)
    #P = u[2]
    #s = u[3]
    residual[1] = u[2] - P_surf
    residual[2] = u[3] - v2_surf
    println("u, out", u)
    println("residual", residual)

end