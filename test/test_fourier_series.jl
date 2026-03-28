# test_fourier_series.jl
# Unit tests for src/fourier_series.jl
#
# Tests cover: construction, getters, copy constructor, location selection,
# frequency selection, merge_by_locations, Base.show, and error handling.

# ── shared test fixtures ───────────────────────────────────────────────────────

# Two-location, three-frequency FourierSeries used across most tests.
function make_test_fs()
    amplitudes = Float32[1.0 2.0 3.0;
                         4.0 5.0 6.0]   # [2 locations × 3 frequencies]
    phases     = Float32[10.0 20.0 30.0;
                         40.0 50.0 60.0]
    frequencies = [0.001, 0.002, 0.003]  # Hz
    return FourierSeries(
        amplitudes, phases, frequencies,
        ["StationA", "StationB"],
        [4.0, 5.0], [52.0, 53.0],
        "water level", "test_source",
    )
end

# ── construction and getters ──────────────────────────────────────────────────

@testset "FourierSeries: construction and getters" begin
    fs = make_test_fs()

    @test fs isa FourierSeries
    @test fs isa AbstractFourierSeries

    @test get_amplitudes(fs)  == Float32[1 2 3; 4 5 6]
    @test get_phases(fs)      == Float32[10 20 30; 40 50 60]
    @test get_frequencies(fs) == [0.001, 0.002, 0.003]
    @test get_names(fs)       == ["StationA", "StationB"]
    @test get_longitudes(fs)  == [4.0, 5.0]
    @test get_latitudes(fs)   == [52.0, 53.0]
    @test get_quantity(fs)    == "water level"
    @test get_source(fs)      == "test_source"

    @test size(get_amplitudes(fs)) == (2, 3)
    @test size(get_phases(fs))     == (2, 3)
    @test length(get_frequencies(fs)) == 3
end

# ── copy constructor ──────────────────────────────────────────────────────────

@testset "FourierSeries: copy constructor" begin
    fs   = make_test_fs()
    fs2  = FourierSeries(fs)

    @test fs2 isa FourierSeries
    @test get_amplitudes(fs2)  == get_amplitudes(fs)
    @test get_phases(fs2)      == get_phases(fs)
    @test get_frequencies(fs2) == get_frequencies(fs)
    @test get_names(fs2)       == get_names(fs)
    @test get_source(fs2)      == get_source(fs)
end

# ── select by location ────────────────────────────────────────────────────────

@testset "FourierSeries: select_location_by_id" begin
    fs  = make_test_fs()
    fs1 = select_location_by_id(fs, 1)

    @test get_names(fs1)       == ["StationA"]
    @test get_amplitudes(fs1)  == Float32[1 2 3]
    @test get_phases(fs1)      == Float32[10 20 30]
    @test get_frequencies(fs1) == get_frequencies(fs)   # frequency axis unchanged
    @test size(get_amplitudes(fs1)) == (1, 3)
end

@testset "FourierSeries: select_locations_by_ids" begin
    fs   = make_test_fs()
    fss  = select_locations_by_ids(fs, [2, 1])

    @test get_names(fss)      == ["StationB", "StationA"]
    @test get_amplitudes(fss) == Float32[4 5 6; 1 2 3]
    @test get_phases(fss)     == Float32[40 50 60; 10 20 30]
    @test size(get_amplitudes(fss)) == (2, 3)
end

@testset "FourierSeries: select_location_by_name" begin
    fs  = make_test_fs()
    fsb = select_location_by_name(fs, "StationB")

    @test get_names(fsb)      == ["StationB"]
    @test get_amplitudes(fsb) == Float32[4 5 6]
    @test get_phases(fsb)     == Float32[40 50 60]
end

@testset "FourierSeries: select_locations_by_names" begin
    fs  = make_test_fs()
    fss = select_locations_by_names(fs, ["StationB", "StationA"])

    @test get_names(fss)      == ["StationB", "StationA"]
    @test get_amplitudes(fss) == Float32[4 5 6; 1 2 3]
end

@testset "FourierSeries: select_locations_by_names — missing location error" begin
    fs = make_test_fs()
    @test_throws ErrorException select_locations_by_names(fs, ["StationX"])
end

# ── select by frequency ───────────────────────────────────────────────────────

@testset "FourierSeries: select_frequencies_by_indices" begin
    fs   = make_test_fs()
    fss  = select_frequencies_by_indices(fs, [1, 3])

    @test get_frequencies(fss) == [0.001, 0.003]
    @test get_amplitudes(fss)  == Float32[1 3; 4 6]
    @test get_phases(fss)      == Float32[10 30; 40 60]
    @test size(get_amplitudes(fss)) == (2, 2)
    @test get_names(fss)       == get_names(fs)        # locations unchanged
end

@testset "FourierSeries: select_frequencies_by_indices — single bin" begin
    fs  = make_test_fs()
    fss = select_frequencies_by_indices(fs, [2])

    @test get_frequencies(fss) == [0.002]
    @test get_amplitudes(fss)  == Float32[2; 5;;]   # 2×1 matrix
    @test size(get_amplitudes(fss)) == (2, 1)
end

# ── merge_by_locations ────────────────────────────────────────────────────────

@testset "FourierSeries: merge_by_locations — round-trip" begin
    fs   = make_test_fs()
    fs1  = select_location_by_id(fs, 1)
    fs2  = select_location_by_id(fs, 2)
    merged = merge_by_locations([fs1, fs2])

    @test merged isa FourierSeries
    @test get_names(merged)       == get_names(fs)
    @test get_amplitudes(merged)  == get_amplitudes(fs)
    @test get_phases(merged)      == get_phases(fs)
    @test get_frequencies(merged) == get_frequencies(fs)
    @test get_quantity(merged)    == get_quantity(fs)
    @test get_source(merged)      == get_source(fs)
end

@testset "FourierSeries: merge_by_locations — empty vector error" begin
    @test_throws ErrorException merge_by_locations(FourierSeries[])
end

@testset "FourierSeries: merge_by_locations — mismatched frequency axis error" begin
    fs1 = make_test_fs() |> x -> select_location_by_id(x, 1)
    fs2 = FourierSeries(
        Float32[4 5 6], Float32[40 50 60], [0.001, 0.002, 0.999],  # different frequencies
        ["StationB"], [5.0], [53.0], "water level", "test_source",
    )
    @test_throws ErrorException merge_by_locations([fs1, fs2])
end

@testset "FourierSeries: merge_by_locations — mismatched quantity error" begin
    fs1 = make_test_fs() |> x -> select_location_by_id(x, 1)
    fs2 = FourierSeries(
        Float32[4 5 6], Float32[40 50 60], [0.001, 0.002, 0.003],
        ["StationB"], [5.0], [53.0], "salinity", "test_source",   # different quantity
    )
    @test_throws ErrorException merge_by_locations([fs1, fs2])
end

@testset "FourierSeries: merge_by_locations — mismatched source error" begin
    fs1 = make_test_fs() |> x -> select_location_by_id(x, 1)
    fs2 = FourierSeries(
        Float32[4 5 6], Float32[40 50 60], [0.001, 0.002, 0.003],
        ["StationB"], [5.0], [53.0], "water level", "other_source",  # different source
    )
    @test_throws ErrorException merge_by_locations([fs1, fs2])
end

# ── Base.show ─────────────────────────────────────────────────────────────────

@testset "FourierSeries: Base.show" begin
    fs  = make_test_fs()
    buf = IOBuffer()

    show(buf, fs)
    out = String(take!(buf))
    @test occursin("FourierSeries", out)
    @test occursin("water level",   out)
    @test occursin("test_source",   out)
    @test occursin("StationA",      out)
    @test occursin("3",             out)   # N frequencies

    show(buf, MIME"text/plain"(), fs)
    out2 = String(take!(buf))
    @test occursin("FourierSeries", out2)
    @test occursin("water level",   out2)
end

# ── single-location edge case ─────────────────────────────────────────────────

@testset "FourierSeries: single location, single frequency" begin
    fs = FourierSeries(
        Float32[7.0;;], Float32[45.0;;], [0.01],
        ["Solo"], [3.0], [51.0], "velocity", "src",
    )
    @test size(get_amplitudes(fs))  == (1, 1)
    @test get_frequencies(fs)       == [0.01]
    @test get_amplitudes(fs)[1, 1]  ≈ 7.0f0
    @test get_phases(fs)[1, 1]      ≈ 45.0f0
end
