# test_statistics.jl
# Unit tests for src/statistics.jl (compute_statistics)
#
# Strategy:
#   • Known-signal tests: construct simple synthetic TimeSeries where all
#     statistics can be computed analytically.
#   • NaN handling: verify that NaN time steps are excluded correctly.
#   • Multiple locations: verify each row in the output DataFrame.
#   • Error tests: mismatched times, mismatched location names.

using DataFrames

# ── helpers ───────────────────────────────────────────────────────────────────

"""Build a minimal TimeSeries from a matrix of values."""
function make_stats_ts(values::Matrix{Float32};
                       t0 = DateTime(2020, 1, 1),
                       dt_s = 3600.0,
                       names = ["S$i" for i in 1:size(values, 1)])
    N     = size(values, 2)
    times = [t0 + Millisecond(round(Int, dt_s * 1000 * (i-1))) for i in 1:N]
    lons  = zeros(Float64, length(names))
    lats  = zeros(Float64, length(names))
    return TimeSeries(values, times, names, lons, lats, "water level", "test")
end

# ── return type ───────────────────────────────────────────────────────────────

@testset "compute_statistics: returns DataFrame" begin
    obs   = make_stats_ts(Float32.([1.0 2.0 3.0]))
    model = make_stats_ts(Float32.([1.5 2.5 3.5]))
    df    = compute_statistics(obs, model)
    @test df isa DataFrame
    @test nrow(df) == 1
    @test names(df) == ["location_id", "location_name", "n_values",
                        "signal_rmse", "bias", "rmse", "mae", "max_error", "min_error"]
end

# ── location metadata ─────────────────────────────────────────────────────────

@testset "compute_statistics: location_id and location_name" begin
    obs   = make_stats_ts(Float32.([1.0 2.0; 3.0 4.0]), names=["A", "B"])
    model = make_stats_ts(Float32.([1.0 2.0; 3.0 4.0]), names=["A", "B"])
    df    = compute_statistics(obs, model)
    @test df.location_id   == [1, 2]
    @test df.location_name == ["A", "B"]
end

# ── n_values ──────────────────────────────────────────────────────────────────

@testset "compute_statistics: n_values counts valid steps" begin
    # 4 time steps, one NaN in obs → 3 valid
    obs_vals   = Float32.([1.0 NaN 3.0 4.0])
    model_vals = Float32.([1.0 2.0 3.0 4.0])
    obs   = make_stats_ts(obs_vals)
    model = make_stats_ts(model_vals)
    df    = compute_statistics(obs, model)
    @test df.n_values[1] == 3
end

# ── perfect model (zero error) ────────────────────────────────────────────────

@testset "compute_statistics: perfect model" begin
    vals  = Float32.([1.0 -2.0 3.0 -1.0])
    obs   = make_stats_ts(vals)
    model = make_stats_ts(vals)
    df    = compute_statistics(obs, model)
    @test df.bias[1]      ≈ 0.0   atol=1e-10
    @test df.rmse[1]      ≈ 0.0   atol=1e-10
    @test df.mae[1]       ≈ 0.0   atol=1e-10
    @test df.max_error[1] ≈ 0.0   atol=1e-10
    @test df.min_error[1] ≈ 0.0   atol=1e-10
end

# ── constant bias ─────────────────────────────────────────────────────────────

@testset "compute_statistics: constant positive bias" begin
    obs   = make_stats_ts(Float32.([1.0 2.0 3.0 4.0]))
    model = make_stats_ts(Float32.([2.0 3.0 4.0 5.0]))   # +1 everywhere
    df    = compute_statistics(obs, model)
    @test df.bias[1]      ≈  1.0   atol=1e-8
    @test df.rmse[1]      ≈  1.0   atol=1e-8
    @test df.mae[1]       ≈  1.0   atol=1e-8
    @test df.max_error[1] ≈  1.0   atol=1e-8
    @test df.min_error[1] ≈  1.0   atol=1e-8
end

@testset "compute_statistics: constant negative bias" begin
    obs   = make_stats_ts(Float32.([2.0 2.0 2.0]))
    model = make_stats_ts(Float32.([1.0 1.0 1.0]))   # −1 everywhere
    df    = compute_statistics(obs, model)
    @test df.bias[1]      ≈ -1.0   atol=1e-8
    @test df.min_error[1] ≈ -1.0   atol=1e-8
    @test df.max_error[1] ≈ -1.0   atol=1e-8
end

# ── signal_rmse ───────────────────────────────────────────────────────────────

@testset "compute_statistics: signal_rmse" begin
    # obs = [3, 4] → rms = sqrt((9+16)/2) = sqrt(12.5)
    obs   = make_stats_ts(Float32.([3.0 4.0]))
    model = make_stats_ts(Float32.([3.0 4.0]))
    df    = compute_statistics(obs, model)
    @test df.signal_rmse[1] ≈ sqrt(12.5)   atol=1e-8
end

# ── analytical RMSE ───────────────────────────────────────────────────────────

@testset "compute_statistics: analytical rmse / mae" begin
    # errors = [−1, +1, −1, +1] → rmse = 1, mae = 1
    obs   = make_stats_ts(Float32.([0.0 0.0 0.0 0.0]))
    model = make_stats_ts(Float32.([-1.0 1.0 -1.0 1.0]))
    df    = compute_statistics(obs, model)
    @test df.rmse[1]      ≈ 1.0   atol=1e-8
    @test df.mae[1]       ≈ 1.0   atol=1e-8
    @test df.bias[1]      ≈ 0.0   atol=1e-8
    @test df.max_error[1] ≈ 1.0   atol=1e-8
    @test df.min_error[1] ≈ -1.0  atol=1e-8
end

# ── NaN handling ──────────────────────────────────────────────────────────────

@testset "compute_statistics: NaN in observations are skipped" begin
    # obs has NaN at index 2; only indices 1, 3 used → bias = 0.5
    obs_vals   = Float32.([0.0 NaN 0.0])
    model_vals = Float32.([1.0 99.0 0.0])   # index 2 should be ignored
    obs   = make_stats_ts(obs_vals)
    model = make_stats_ts(model_vals)
    df    = compute_statistics(obs, model)
    # errors at valid indices: 1.0, 0.0 → bias = 0.5
    @test df.bias[1] ≈ 0.5   atol=1e-8
end

@testset "compute_statistics: NaN in model are skipped" begin
    obs_vals   = Float32.([0.0 0.0 0.0])
    model_vals = Float32.([1.0 NaN 1.0])   # index 2 NaN
    obs   = make_stats_ts(obs_vals)
    model = make_stats_ts(model_vals)
    df    = compute_statistics(obs, model)
    # valid errors: 1.0, 1.0 → bias = 1.0
    @test df.bias[1] ≈ 1.0   atol=1e-8
end

# ── multiple locations ────────────────────────────────────────────────────────

@testset "compute_statistics: multiple locations" begin
    # loc 1: bias +1; loc 2: bias −2
    obs_vals   = Float32.([0.0 0.0; 0.0 0.0])
    model_vals = Float32.([1.0 1.0; -2.0 -2.0])
    obs   = make_stats_ts(obs_vals, names=["A","B"])
    model = make_stats_ts(model_vals, names=["A","B"])
    df    = compute_statistics(obs, model)
    @test nrow(df) == 2
    @test df.bias[1] ≈  1.0   atol=1e-8
    @test df.bias[2] ≈ -2.0   atol=1e-8
end

# ── error: mismatched time axes ───────────────────────────────────────────────

@testset "compute_statistics: error on mismatched time length" begin
    obs   = make_stats_ts(Float32.([1.0 2.0 3.0]))          # 3 steps
    model = make_stats_ts(Float32.([1.0 2.0 3.0 4.0]))      # 4 steps
    @test_throws ErrorException compute_statistics(obs, model)
end

@testset "compute_statistics: error on shifted times" begin
    vals  = Float32.([1.0 2.0 3.0])
    obs   = make_stats_ts(vals, t0=DateTime(2020,1,1))
    model = make_stats_ts(vals, t0=DateTime(2020,1,2))   # one day later
    @test_throws ErrorException compute_statistics(obs, model)
end

# ── error: mismatched location names ─────────────────────────────────────────

@testset "compute_statistics: error on mismatched location names" begin
    obs   = make_stats_ts(Float32.([1.0 2.0]), names=["A"])
    model = make_stats_ts(Float32.([1.0 2.0]), names=["B"])
    @test_throws ErrorException compute_statistics(obs, model)
end

@testset "compute_statistics: error on mismatched location count" begin
    obs   = make_stats_ts(Float32.([1.0 2.0; 3.0 4.0]), names=["A","B"])
    model = make_stats_ts(Float32.([1.0 2.0]), names=["A"])
    @test_throws ErrorException compute_statistics(obs, model)
end
