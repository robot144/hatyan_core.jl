# test_plotting.jl
# Tests for src/plotting.jl
#
# Strategy
# ────────
# • Verify that each plot() method returns a Plots.Plot without error.
# • Save every plot to test/temp/ so they can be inspected after the run.
# • Cover single-location and multi-location inputs.
# • Cover optional keyword arguments (yunit, xscale, freq_unit,
#   max_constituents, location_index).
#
# test/temp/ is wiped at the start of every run so stale files never accumulate.

import Plots
Plots.gr()   # headless GR backend — no display required

# ── output directory: wipe and recreate on every run ──────────────────────────

const PLOT_DIR = joinpath(@__DIR__, "temp")
rm(PLOT_DIR; recursive=true, force=true)
mkpath(PLOT_DIR)

# ── shared fixtures ────────────────────────────────────────────────────────────

# 30-day TimeSeries from the real VLISSGN constituent file (single location).
# 30 days gives a frequency resolution of 1/30 cpd ≈ 0.033 cpd so tidal
# peaks (M2 ≈ 1.93 cpd, S2 = 2 cpd, K1 ≈ 1 cpd) are well resolved.
function _ts_vlissgn()
    tc    = read_donar_constituents(joinpath(TEST_DATA_DIR, "VLISSGN_ana.txt"))
    times = [DateTime(2019, 1, 1) + Minute(10*i) for i in 0:(30*144-1)]   # 4320 steps
    return prediction(tc, times, HatyanSettings(nodalfactors=true, fu_alltimes=true))
end

# Two-location TimeSeries built from the same data
function _ts_two_locations()
    tc = read_donar_constituents(joinpath(TEST_DATA_DIR, "VLISSGN_ana.txt"))
    amp2 = tc.amplitudes .* 0.7f0
    phi2 = mod.(tc.phases .+ 15.0f0, 360.0f0)
    tc2  = TidalConstituents(amp2, phi2, tc.constituent_names,
                              ["LOC2"], [5.0], [53.0], tc.quantity, tc.source)
    tc_both = TidalConstituents(
        vcat(tc.amplitudes, tc2.amplitudes),
        vcat(tc.phases,     tc2.phases),
        tc.constituent_names,
        ["VLISSGN", "LOC2"],
        vcat(tc.longitudes, tc2.longitudes),
        vcat(tc.latitudes,  tc2.latitudes),
        tc.quantity, tc.source,
    )
    times = [DateTime(2019, 1, 1) + Minute(10*i) for i in 0:(30*144-1)]
    return prediction(tc_both, times, HatyanSettings(nodalfactors=true, fu_alltimes=true))
end

# Two-location TidalConstituents
function _tc_two_locations()
    tc1  = read_donar_constituents(joinpath(TEST_DATA_DIR, "VLISSGN_ana.txt"))
    amp2 = tc1.amplitudes .* 0.8f0
    phi2 = mod.(tc1.phases .+ 20.0f0, 360.0f0)
    tc2  = TidalConstituents(amp2, phi2, tc1.constituent_names,
                              ["LOC2"], [5.0], [53.0], tc1.quantity, tc1.source)
    return TidalConstituents(
        vcat(tc1.amplitudes, tc2.amplitudes),
        vcat(tc1.phases,     tc2.phases),
        tc1.constituent_names,
        ["VLISSGN", "LOC2"],
        vcat(tc1.longitudes, tc2.longitudes),
        vcat(tc1.latitudes,  tc2.latitudes),
        tc1.quantity, tc1.source,
    )
end

# ── helper ─────────────────────────────────────────────────────────────────────

"""Save `p` to `test/temp/<name>.png` and assert the file is non-empty."""
function save_and_check(p::Plots.Plot, name::String)
    path = joinpath(PLOT_DIR, name * ".png")
    Plots.savefig(p, path)
    @test isfile(path)
    @test filesize(path) > 1_000
end

# ── TimeSeries plots ───────────────────────────────────────────────────────────

@testset "plot(TimeSeries): single location" begin
    ts = _ts_vlissgn()
    p  = Plots.plot(ts)
    @test p isa Plots.Plot
    save_and_check(p, "ts_single")
end

@testset "plot(TimeSeries): location_index selects one location" begin
    ts = _ts_two_locations()
    p1 = Plots.plot(ts; location_index=1)
    @test p1 isa Plots.Plot
    save_and_check(p1, "ts_loc1")

    p2 = Plots.plot(ts; location_index=2)
    @test p2 isa Plots.Plot
    save_and_check(p2, "ts_loc2")
end

@testset "plot(TimeSeries): location_index=nothing overlays all locations" begin
    ts = _ts_two_locations()
    p  = Plots.plot(ts; location_index=nothing)
    @test p isa Plots.Plot
    save_and_check(p, "ts_all_locs")
end

@testset "plot(TimeSeries): out-of-range location_index errors" begin
    ts = _ts_vlissgn()   # 1 location
    @test_throws ErrorException Plots.plot(ts; location_index=2)
end

@testset "plot(TimeSeries): yunit keyword" begin
    ts = _ts_vlissgn()
    p  = Plots.plot(ts; yunit="m")
    @test p isa Plots.Plot
    save_and_check(p, "ts_yunit")
end

@testset "plot(TimeSeries): size and dpi kwargs forwarded" begin
    ts = _ts_vlissgn()
    p  = Plots.plot(ts; size=(800, 400), dpi=150)
    save_and_check(p, "ts_custom_size")
end

# ── TidalConstituents plots ────────────────────────────────────────────────────

@testset "plot(TidalConstituents): single location" begin
    tc = read_donar_constituents(joinpath(TEST_DATA_DIR, "VLISSGN_ana.txt"))
    p  = Plots.plot(tc)
    @test p isa Plots.Plot
    save_and_check(p, "tc_single")
end

@testset "plot(TidalConstituents): default shows top 20 by amplitude" begin
    tc = read_donar_constituents(joinpath(TEST_DATA_DIR, "VLISSGN_ana.txt"))
    # M2 is the dominant constituent at Vlissingen — it must appear in the default view
    @test "M2" in get_constituent_names(tc)
    p = Plots.plot(tc)
    @test p isa Plots.Plot
    save_and_check(p, "tc_default_top20")
end

@testset "plot(TidalConstituents): max_constituents selects top N by amplitude" begin
    tc = read_donar_constituents(joinpath(TEST_DATA_DIR, "VLISSGN_ana.txt"))
    p  = Plots.plot(tc; max_constituents=5)
    @test p isa Plots.Plot
    save_and_check(p, "tc_top5")
end

@testset "plot(TidalConstituents): location_index selects one location" begin
    tc_both = _tc_two_locations()
    p1 = Plots.plot(tc_both; location_index=1, max_constituents=20)
    @test p1 isa Plots.Plot
    save_and_check(p1, "tc_loc1")

    p2 = Plots.plot(tc_both; location_index=2, max_constituents=20)
    @test p2 isa Plots.Plot
    save_and_check(p2, "tc_loc2")
end

@testset "plot(TidalConstituents): location_index=nothing shows all locations" begin
    tc_both = _tc_two_locations()
    p = Plots.plot(tc_both; location_index=nothing, max_constituents=20)
    @test p isa Plots.Plot
    save_and_check(p, "tc_all_locs")
end

@testset "plot(TidalConstituents): out-of-range location_index errors" begin
    tc = read_donar_constituents(joinpath(TEST_DATA_DIR, "VLISSGN_ana.txt"))
    @test_throws ErrorException Plots.plot(tc; location_index=2)
end

# ── FourierSeries plots ────────────────────────────────────────────────────────

@testset "plot(FourierSeries): single location" begin
    fs = fft_series(_ts_vlissgn())
    p  = Plots.plot(fs)
    @test p isa Plots.Plot
    save_and_check(p, "fs_single")
end

@testset "plot(FourierSeries): freq_max limits x-axis" begin
    fs = fft_series(_ts_vlissgn())
    p  = Plots.plot(fs; freq_max=5.0)   # zoom into sub-diurnal band
    save_and_check(p, "fs_freq_max")
end

@testset "plot(FourierSeries): log x-scale" begin
    fs = fft_series(_ts_vlissgn())
    p  = Plots.plot(fs; xscale=:log10)
    save_and_check(p, "fs_logx")
end

@testset "plot(FourierSeries): cycles-per-hour frequency unit" begin
    fs = fft_series(_ts_vlissgn())
    p  = Plots.plot(fs; freq_unit="cph")
    save_and_check(p, "fs_cph")
end

@testset "plot(FourierSeries): cycles-per-day frequency unit" begin
    fs = fft_series(_ts_vlissgn())
    p  = Plots.plot(fs; freq_unit="cpd")
    save_and_check(p, "fs_cpd")
end

@testset "plot(FourierSeries): location_index selects one location" begin
    fs = fft_series(_ts_two_locations())
    p1 = Plots.plot(fs; location_index=1)
    @test p1 isa Plots.Plot
    save_and_check(p1, "fs_loc1")

    p2 = Plots.plot(fs; location_index=2)
    @test p2 isa Plots.Plot
    save_and_check(p2, "fs_loc2")
end

@testset "plot(FourierSeries): location_index=nothing shows all locations" begin
    fs = fft_series(_ts_two_locations())
    p  = Plots.plot(fs; location_index=nothing)
    save_and_check(p, "fs_all_locs")
end

@testset "plot(FourierSeries): out-of-range location_index errors" begin
    fs = fft_series(_ts_vlissgn())   # 1 location
    @test_throws ErrorException Plots.plot(fs; location_index=2)
end
