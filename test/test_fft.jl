# test_fft.jl
# Unit tests for src/fft.jl  (fft_series / ifft_series)
#
# Strategy:
#   • Known-signal tests: construct a TimeSeries from pure sinusoids whose
#     amplitude and phase are analytically known, then verify the spectrum.
#   • Round-trip tests: fft_series → ifft_series must recover the original
#     values to within Float32 precision.
#   • Edge cases: even/odd N, DC-only, single location vs multiple locations.
#   • Error tests: non-uniform times, mismatched bin count.

# ── helpers ───────────────────────────────────────────────────────────────────

"""Build a minimal TimeSeries from a matrix of values and a time step (seconds)."""
function make_ts(values::Matrix{Float32}, dt_s::Float64;
                 t0 = DateTime(2020, 1, 1),
                 names = ["S$i" for i in 1:size(values,1)])
    N     = size(values, 2)
    times = [t0 + Millisecond(round(Int, dt_s * 1000 * (i-1))) for i in 1:N]
    lons  = zeros(Float64, length(names))
    lats  = zeros(Float64, length(names))
    return TimeSeries(values, times, names, lons, lats, "water level", "test")
end

"""Return the 1-based Julia index in `freqs` closest to `target_hz`."""
function freq_index(freqs, target_hz)
    return argmin(abs.(freqs .- target_hz))
end

"""
Signal that is exactly the k-th FFT basis function (0-based bin k).
This avoids floating-point disagreement between independently-computed Hz values.
"""
function pure_cosine(N::Int, k::Int, A::Float32, phi_deg::Float64)
    return A .* cos.(2π .* k ./ N .* (0:N-1) .+ deg2rad(phi_deg))
end

# ── output type ───────────────────────────────────────────────────────────────

@testset "fft_series: returns FourierSeries" begin
    ts = make_ts(ones(Float32, 1, 8), 1.0)
    fs = fft_series(ts)
    @test fs isa FourierSeries
    @test fs isa AbstractFourierSeries
end

# ── frequency axis ────────────────────────────────────────────────────────────

@testset "fft_series: frequency axis (even N)" begin
    N = 8; dt = 1.0   # 1 Hz sampling, 8 points
    ts = make_ts(zeros(Float32, 1, N), dt)
    fs = fft_series(ts)

    freqs = get_frequencies(fs)
    @test length(freqs) == N ÷ 2 + 1   # 5 bins
    @test freqs[1] ≈ 0.0               # DC
    @test freqs[2] ≈ 1.0 / (N * dt)   # fundamental
    @test freqs[end] ≈ 0.5 / dt        # Nyquist = 0.5 Hz
end

@testset "fft_series: frequency axis (odd N)" begin
    N = 9; dt = 1.0
    ts = make_ts(zeros(Float32, 1, N), dt)
    fs = fft_series(ts)

    freqs = get_frequencies(fs)
    @test length(freqs) == N ÷ 2 + 1   # 5 bins
    @test freqs[1] ≈ 0.0
    @test freqs[2] ≈ 1.0 / (N * dt)
end

# ── DC component ──────────────────────────────────────────────────────────────

@testset "fft_series: DC (constant signal)" begin
    N = 16; A = 2.5f0
    ts = make_ts(fill(A, 1, N), 1.0)
    fs = fft_series(ts)

    @test get_amplitudes(fs)[1, 1] ≈ A   atol=1e-5   # DC bin == mean
    @test all(abs.(get_amplitudes(fs)[1, 2:end]) .< 1e-5)  # no other energy
end

# ── pure cosine ───────────────────────────────────────────────────────────────

@testset "fft_series: pure cosine (even N, integer cycles)" begin
    N  = 64; dt = 1.0
    A  = 3.0f0; φ = 30.0   # degrees
    k0 = 4   # 0-based FFT bin (4 complete cycles in window)

    sig = pure_cosine(N, k0, A, φ)
    ts  = make_ts(Float32.(reshape(sig, 1, N)), dt)
    fs  = fft_series(ts)

    @test get_amplitudes(fs)[1, k0+1] ≈ A          atol=1e-4
    @test get_phases(fs)[1, k0+1]     ≈ Float32(φ) atol=1e-3
end

@testset "fft_series: pure cosine (odd N)" begin
    N  = 63; dt = 1.0
    A  = 1.5f0; φ = -45.0
    k0 = 2   # 0-based FFT bin

    # cos(2π*k0/N*n + φ) is exactly the k0-th DFT basis function
    sig = pure_cosine(N, k0, A, φ)
    ts  = make_ts(Float32.(reshape(sig, 1, N)), dt)
    fs  = fft_series(ts)

    @test get_amplitudes(fs)[1, k0+1] ≈ A          atol=1e-4
    @test get_phases(fs)[1, k0+1]     ≈ Float32(φ) atol=1e-3
end

# ── two sinusoids ─────────────────────────────────────────────────────────────

@testset "fft_series: two superimposed cosines" begin
    N = 128; dt = 600.0   # 10-minute tidal steps
    k1 = 1; A1 = 2.0f0; φ1 = 0.0
    k2 = 3; A2 = 0.5f0; φ2 = 90.0

    sig = pure_cosine(N, k1, A1, φ1) .+ pure_cosine(N, k2, A2, φ2)
    ts  = make_ts(Float32.(reshape(sig, 1, N)), dt)
    fs  = fft_series(ts)

    @test get_amplitudes(fs)[1, k1+1] ≈ A1   atol=1e-4
    @test get_amplitudes(fs)[1, k2+1] ≈ A2   atol=1e-4
end

# ── metadata propagation ──────────────────────────────────────────────────────

@testset "fft_series: metadata" begin
    ts = make_ts(ones(Float32, 1, 8), 1.0)
    fs = fft_series(ts)

    @test get_names(fs)      == get_names(ts)
    @test get_longitudes(fs) == get_longitudes(ts)
    @test get_latitudes(fs)  == get_latitudes(ts)
    @test get_quantity(fs)   == get_quantity(ts)
    @test occursin("fft",    get_source(fs))
    @test occursin("test",   get_source(fs))
end

# ── multiple locations ────────────────────────────────────────────────────────

@testset "fft_series: multiple locations" begin
    N = 32; dt = 1.0
    k0 = 2   # 0-based FFT bin
    A  = [1.0f0, 3.0f0]

    vals = Float32.(vcat(pure_cosine(N, k0, A[1], 0.0)',
                         pure_cosine(N, k0, A[2], 0.0)'))
    ts = make_ts(vals, dt, names=["L1","L2"])
    fs = fft_series(ts)

    @test size(get_amplitudes(fs)) == (2, N÷2+1)
    @test get_amplitudes(fs)[1, k0+1] ≈ A[1]   atol=1e-4
    @test get_amplitudes(fs)[2, k0+1] ≈ A[2]   atol=1e-4
end

# ── round-trip: fft → ifft ────────────────────────────────────────────────────

@testset "ifft_series: round-trip (even N)" begin
    N  = 64; dt = 600.0
    f0 = 3.0 / (N * dt)
    t  = (0:N-1) .* dt
    sig = Float32.(2.0 .* cos.(2π .* f0 .* t) .+ 0.5)
    ts  = make_ts(reshape(sig, 1, N), dt)

    fs    = fft_series(ts)
    ts_rt = ifft_series(fs, get_times(ts))

    @test ts_rt isa TimeSeries
    @test maximum(abs.(get_values(ts_rt) .- get_values(ts))) < 1e-4
end

@testset "ifft_series: round-trip (odd N)" begin
    N  = 63; dt = 1.0
    f0 = 2.0 / (N * dt)
    t  = (0:N-1) .* dt
    sig = Float32.(1.5 .* cos.(2π .* f0 .* t .+ deg2rad(45.0)))
    ts  = make_ts(reshape(sig, 1, N), dt)

    ts_rt = ifft_series(fft_series(ts), get_times(ts))
    @test maximum(abs.(get_values(ts_rt) .- get_values(ts))) < 1e-4
end

@testset "ifft_series: round-trip (multiple locations)" begin
    N = 32; dt = 1.0
    f0 = 1.0 / (N * dt)
    t  = (0:N-1) .* dt
    v1 = Float32.(2.0 .* cos.(2π .* f0 .* t))
    v2 = Float32.(0.3 .* cos.(2π .* 3f0 .* t .+ deg2rad(120.0)))
    vals = Float32.(hcat(v1, v2)')
    ts   = make_ts(vals, dt, names=["A","B"])

    ts_rt = ifft_series(fft_series(ts), get_times(ts))
    @test maximum(abs.(get_values(ts_rt) .- get_values(ts))) < 1e-4
end

# ── ifft_series metadata ──────────────────────────────────────────────────────

@testset "ifft_series: metadata" begin
    ts = make_ts(ones(Float32, 1, 8), 1.0)
    fs = fft_series(ts)
    ts_rt = ifft_series(fs, get_times(ts))

    @test get_names(ts_rt)      == get_names(ts)
    @test get_longitudes(ts_rt) == get_longitudes(ts)
    @test get_latitudes(ts_rt)  == get_latitudes(ts)
    @test get_quantity(ts_rt)   == get_quantity(ts)
    @test occursin("ifft",      get_source(ts_rt))
    @test occursin("fft",       get_source(ts_rt))
    @test occursin("test",      get_source(ts_rt))
end

# ── error handling ────────────────────────────────────────────────────────────

@testset "fft_series: non-uniform times error" begin
    N  = 8
    t0 = DateTime(2020)
    # irregular gap at step 3
    times = [t0 + Second(i) for i in 0:N-1]
    times[4] = times[4] + Second(2)
    ts = TimeSeries(ones(Float32, 1, N), times, ["S"], [0.0], [0.0],
                    "water level", "test")
    @test_throws ErrorException fft_series(ts)
end

@testset "ifft_series: bin count mismatch error" begin
    ts = make_ts(ones(Float32, 1, 8), 1.0)
    fs = fft_series(ts)   # 8-point FFT → 5 bins

    # Build a times vector of length 10 → expects 6 bins
    times_10 = [DateTime(2020) + Second(i) for i in 0:9]
    @test_throws ErrorException ifft_series(fs, times_10)
end
