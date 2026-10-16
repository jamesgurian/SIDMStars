using Test
using SIDMStars
using SIDMStars: K5, K7, K9, Keff2, calc_κ_log

# InterpolatedCS built from a power-law shape σ̂(v²) = (v²/w²)^(n/2) must reproduce
# the closed-form PowerLawCS K5, K_eff^(2) and calc_κ_log to ≲1e−4, inside the
# table and in the extrapolation region.

G  = 4.30091727003628e-6
C  = 0.8505
σ0 = 1e-6
w  = 10.0

common_kwargs = (t_end = 1.0, r_max = 1.0, r_min = 1e-3, gridpoints = 10,
                  snap_interval_time = Inf, snap_interval_density = Inf, snap_loc = "")

# ln v² offsets (nats) relative to ln(v2_center): a mix of points well inside
# the default ±45-nat table half-width and points beyond it (extrapolation).
offsets = [-50.0, -44.0, -30.0, -10.0, -1.0, 0.0, 1.0, 10.0, 30.0, 44.0, 50.0]

# several decades of ln ρ
lρ_grid = range(-20.0, 20.0, length=9)

@testset "InterpolatedCS vs PowerLawCS (n=$n)" for n in (0.0, -2.0, -4.0)
    v2_center = w^2
    σ̂(v2) = (v2 / w^2)^(n / 2)
    cs_interp = InterpolatedCS(σ̂, σ0; v2_center = v2_center)

    rp_pl = RunParams(; G = G, C = C, σ0 = σ0, w = w, n = n, common_kwargs...)
    rp_it = RunParams(; G = G, C = C, σ0 = σ0, w = w, n = NaN, cs = cs_interp, common_kwargs...)

    @test rp_pl.cs isa PowerLawCS
    @test rp_it.cs isa InterpolatedCS
    @test isfinite(rp_pl.lci_s) && isfinite(rp_pl.lci_l)
    @test isnan(rp_it.lci_s)  # power-law field, unused by InterpolatedCS

    lv2_grid = [log(v2_center) + Δ for Δ in offsets]

    @testset "K5, Keff2 at lv2 offset $Δ" for (lv2, Δ) in zip(lv2_grid, offsets)
        v2 = exp(lv2)

        k5_pl = K5(v2, rp_pl)
        k5_it = K5(v2, rp_it)
        @test isapprox(k5_it, k5_pl; rtol = 1e-4)

        keff2_pl = Keff2(v2, rp_pl)
        keff2_it = Keff2(v2, rp_it)
        @test isapprox(keff2_it, keff2_pl; rtol = 1e-4)
    end

    @testset "calc_κ_log at lv2 offset $Δ, lρ=$lρ" for (lv2, Δ) in zip(lv2_grid, offsets), lρ in lρ_grid
        κinv_pl = calc_κ_log(lρ, lv2, rp_pl)
        κinv_it = calc_κ_log(lρ, lv2, rp_it)
        @test isapprox(κinv_it, κinv_pl; rtol = 1e-4)
    end
end

# LMFP K5/K7 mix (lmfp_K7_frac): the power-law (constant-factor) and interpolated (per-v² table)
# paths must agree, a = 0 must be bitwise the K5 closure, and in the deep LMFP κ must scale by
# K_lmfp/K5 = 1 + a·n/8 exactly (power law).
@testset "lmfp_K7_frac mix (n=$n)" for n in (0.0, -2.0, -4.0)
    a7 = 0.73
    σ̂(v2) = (v2 / w^2)^(n / 2)
    cs_interp = InterpolatedCS(σ̂, σ0; v2_center = w^2)
    rp_pl0 = RunParams(; G = G, C = C, σ0 = σ0, w = w, n = n, common_kwargs...)
    rp_plz = RunParams(; G = G, C = C, σ0 = σ0, w = w, n = n, lmfp_K7_frac = 0.0, common_kwargs...)
    rp_pl  = RunParams(; G = G, C = C, σ0 = σ0, w = w, n = n, lmfp_K7_frac = a7, common_kwargs...)
    rp_it  = RunParams(; G = G, C = C, σ0 = σ0, w = w, n = NaN, cs = cs_interp, lmfp_K7_frac = a7,
                         common_kwargs...)
    @test rp_plz.lci_l === rp_pl0.lci_l
    for Δ in (-30.0, -1.0, 0.0, 1.0, 30.0), lρ in (-20.0, 0.0, 20.0)
        lv2 = log(w^2) + Δ
        @test isapprox(K7(exp(lv2), rp_it), K7(exp(lv2), rp_pl); rtol = 1e-4)
        @test isapprox(calc_κ_log(lρ, lv2, rp_it), calc_κ_log(lρ, lv2, rp_pl); rtol = 1e-4)
    end
    # deep LMFP (small lρ, Kinv_lmfp ≫ Kinv_smfp): κ ∝ K_lmfp
    lv2 = log(w^2)
    @test isapprox(calc_κ_log(-60.0, lv2, rp_pl) / calc_κ_log(-60.0, lv2, rp_pl0), 1 + a7 * n / 8; rtol = 1e-6)
end

# Sanity check on the convenience constructor and the evolve() guard invariant
# (evolve.jl errors if cs isa PowerLawCS && !isfinite(lci_s); RegularizedCS
# runs must never trip that guard the other way — lb_s/lb_l must be finite).
@testset "RegularizedCS construction" begin
    reg = RegularizedCS(σ0, w)
    @test reg isa InterpolatedCS
    rp_reg = RunParams(; G = G, C = C, σ0 = σ0, w = w, n = NaN, cs = reg, common_kwargs...)
    @test isfinite(rp_reg.lb_s) && isfinite(rp_reg.lb_l)
    @test isfinite(calc_κ_log(0.0, log(w^2), rp_reg))
end

# --- End-to-end evolve() runs, demonstrating actual usage of the general
# CrossSection dispatch (not just the moment-level checks above). ---

r_max_run = 30.
prof = x -> nfw_profile(1.18, 2.73e7, x)
P_surf = calc_psurf(prof, r_max_run, G)

run_kwargs = (convection = ConvectionOff(), evaporation = EvaporationOff(),
              G = G, C = C, t_end = 0.05, r_max = r_max_run, r_min = 1e-3,
              gridpoints = 60, snap_interval_time = 0.01, snap_interval_density = 1.0)

# An InterpolatedCS built from a power-law shape must reproduce a PowerLawCS
# run to the InterpolatedCS table's precision, not just at the moment level
# (crosssection_test.jl above) but end-to-end through the shooting BVP /
# conduction solve in evolve().
@testset "InterpolatedCS (power-law shape) matches PowerLawCS in evolve()" begin
    n_run = -2.0
    σ̂(v2) = (v2 / w^2)^(n_run / 2)
    cs_interp = InterpolatedCS(σ̂, σ0; v2_center = w^2)

    snap_pl = tempname() * ".hdf5"
    snap_it = tempname() * ".hdf5"
    try
        rp_pl = RunParams(; σ0 = σ0, w = w, n = n_run, snap_loc = snap_pl, run_kwargs...)
        rp_it = RunParams(; σ0 = σ0, w = w, n = NaN, cs = cs_interp, snap_loc = snap_it, run_kwargs...)

        evolve(rp_pl, prof, eps_kh = .05, eps_diff = 0.2, mass_min_frac = 1e-6, p_outer = P_surf)
        evolve(rp_it, prof, eps_kh = .05, eps_diff = 0.2, mass_min_frac = 1e-6, p_outer = P_surf)

        snap1 = read_snapshot_file(snap_pl)
        snap2 = read_snapshot_file(snap_it)
        j1, j2 = size(snap1["rho"], 2), size(snap2["rho"], 2)
        @test j1 > 1 && j2 > 1

        @test isapprox(snap2["rho"][:, j2], snap1["rho"][:, j1]; rtol = 1e-3)
        @test isapprox(snap2["r"][:, j2], snap1["r"][:, j1]; rtol = 1e-3)
    finally
        rm(snap_pl, force = true)
        rm(snap_it, force = true)
    end
end

# RegularizedCS has no closed-form comparison (that's the point — it's the
# shape used by the companion KiSS-SIDM runs), so this is a smoke test: evolve()
# must complete and hand back a finite, positive, radially sorted profile.
@testset "RegularizedCS evolve() smoke test" begin
    reg = RegularizedCS(σ0, w)
    snap_reg = tempname() * ".hdf5"
    try
        rp_reg = RunParams(; σ0 = σ0, w = w, n = NaN, cs = reg, snap_loc = snap_reg, run_kwargs...)
        evolve(rp_reg, prof, eps_kh = .05, eps_diff = 0.2, mass_min_frac = 1e-6, p_outer = P_surf)

        snap = read_snapshot_file(snap_reg)
        j = size(snap["rho"], 2)
        @test j > 1
        @test all(isfinite, snap["rho"][:, j])
        @test all(>(0), snap["rho"][:, j])
        @test issorted(snap["r"][:, j])
    finally
        rm(snap_reg, force = true)
    end
end

println("crosssection_test.jl: all tests passed.")
