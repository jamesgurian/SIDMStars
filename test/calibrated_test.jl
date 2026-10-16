using Test
using SIDMStars

# calibrated_params fills in the adopted conductivity calibration; calib picks the
# LMFP velocity weighting (:smfp → Γ₅, :lmfp → Γ₅/Γ₇ mixture).

common = (t_end = 100.0, r_max = 40.0, snap_loc = "")
σ0 = 1e-6

@testset "calibrated_params" begin
    iso = calibrated_params(σ0; common...)
    @test iso.cs isa PowerLawCS && iso.n == 0.0 && iso.σ0 == σ0
    @test (iso.C, iso.bump_A, iso.bump_W) == (0.8505, 1.31, 1.68)
    # The library default is the calibrated C.
    @test RunParams(G = 1.0, σ0 = σ0, w = 1.0, n = 0.0; common..., r_min = 1e-3,
                    gridpoints = 10, snap_interval_time = Inf,
                    snap_interval_density = Inf).C == 0.8505
    @test iso.lmfp_K7_frac == 0.0
    @test calibrated_params(σ0; common..., calib = :lmfp).lmfp_K7_frac == 0.73
    @test iso.evaporation isa EvaporationOff && iso.convection isa ConvectionOn
    @test iso.snap_interval_time == 0.1 && iso.r_min ≈ 40e-6

    plcs = PowerLawCS(σ0 = σ0, w = 10.0, n = -4.0)
    @test calibrated_params(plcs; common...).lmfp_K7_frac == 0.0
    pl = calibrated_params(plcs; common..., calib = :lmfp)
    @test pl.lmfp_K7_frac == 0.73 && pl.w == 10.0 && pl.n == -4.0
    @test (pl.C, pl.bump_A, pl.bump_W) == (0.8505, 1.31, 1.68)

    regcs = RegularizedCS(σ0, 10.0)
    reg = calibrated_params(regcs; common..., calib = :lmfp)
    @test reg.lmfp_K7_frac == 0.73 && isnan(reg.n) && reg.σ0 == σ0
    @test calibrated_params(regcs; common..., calib = :smfp).lmfp_K7_frac == 0.0
    @test_throws ArgumentError calibrated_params(regcs; common..., calib = :imfp)

    # Overrides and pass-through keywords.
    ov = calibrated_params(regcs; common..., calib = :lmfp, lmfp_K7_frac = 0.0,
                           C = 0.9, collapse_density_factor = 100.0)
    @test ov.lmfp_K7_frac == 0.0 && ov.C == 0.9 && ov.collapse_density_factor == 100.0

    # Same RunParams as the long-hand construction.
    ref = RunParams(G = 4.30091727003628e-6, C = 0.8505, σ0 = σ0, w = 10.0, n = -4.0,
                    bump_A = 1.31, bump_W = 1.68, lmfp_K7_frac = 0.73,
                    t_end = 100.0, r_max = 40.0, r_min = 1e-6 * 40.0, gridpoints = 200,
                    snap_interval_time = 0.1, snap_interval_density = 0.1, snap_loc = "",
                    evaporation = EvaporationOff())
    @test all(getfield(pl, f) === getfield(ref, f) for f in fieldnames(typeof(ref)))
end
