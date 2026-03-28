# fft.jl
#
# Forward and inverse FFT conversion between TimeSeries and FourierSeries.
#
# Conventions
# ───────────
# • Uses FFTW.rfft / FFTW.irfft (real-valued, one-sided spectrum).
# • Amplitude-preserving normalisation: a pure tone x(t) = A·cos(2πft + φ)
#   maps to amplitude A (metres) and phase φ (degrees) at frequency f.
# • Frequency axis: Hz, using FFTW.rfftfreq(N, 1/dt) — rfftfreq takes the
#   sampling rate (Hz), not the sample spacing.
# • DC bin (k=0) and Nyquist bin (k=N÷2, even N only) are NOT doubled in the
#   amplitude; all other bins represent the one-sided spectrum with factor 2.
#
# Round-trip guarantee
# ────────────────────
# For a TimeSeries ts with N uniformly-spaced steps:
#   ifft_series(fft_series(ts), get_times(ts)) ≈ ts   (float32 precision)

using FFTW

# ── helpers ────────────────────────────────────────────────────────────────────

"""
    _check_uniform_times(times)

Throw an error when `times` are not uniformly spaced (to within 1 ms tolerance).
"""
function _check_uniform_times(times::Vector{DateTime})
    N = length(times)
    N < 2 && error("TimeSeries must have at least 2 time steps for FFT.")
    expected_dt = times[2] - times[1]
    for i in 2:N-1
        actual = times[i+1] - times[i]
        if abs((actual - expected_dt).value) > 1   # 1 ms tolerance
            error("FFT requires uniformly-spaced times; step at index $i differs " *
                  "from the first step ($(expected_dt) vs $(actual)).")
        end
    end
end

"""
    _dt_seconds(times) -> Float64

Return the uniform time step of `times` in seconds.
"""
function _dt_seconds(times::Vector{DateTime})::Float64
    return (times[2] - times[1]).value / 1000.0   # Dates.Millisecond → seconds
end

# ── forward FFT ───────────────────────────────────────────────────────────────

"""
    fft_series(ts) -> FourierSeries

Compute the one-sided amplitude/phase spectrum of each location in `ts`.

The amplitude at each frequency bin equals the peak amplitude of the
corresponding sinusoidal component (metres for water-level data).
Phases are in degrees.  The frequency axis is in Hz.

The input `ts` must have uniformly-spaced time steps.
"""
function fft_series(ts::AbstractTimeSeries)::FourierSeries
    times  = get_times(ts)
    vals   = get_values(ts)          # Float32 [locations × times]
    N      = length(times)

    _check_uniform_times(times)
    dt = _dt_seconds(times)

    n_locs  = size(vals, 1)
    n_freqs = N ÷ 2 + 1

    amplitudes = Matrix{Float32}(undef, n_locs, n_freqs)
    phases     = Matrix{Float32}(undef, n_locs, n_freqs)

    for i in 1:n_locs
        X = rfft(Float64.(vals[i, :]))    # length n_freqs complex vector

        # Amplitude: 2/N for interior bins, 1/N for DC and (even-N) Nyquist
        a = abs.(X) .* (2.0 / N)
        a[1] *= 0.5                       # DC
        if N % 2 == 0
            a[end] *= 0.5                 # Nyquist
        end

        amplitudes[i, :] = Float32.(a)
        phases[i, :]     = Float32.(rad2deg.(angle.(X)))
    end

    frequencies = collect(rfftfreq(N, inv(dt)))   # Hz, length n_freqs

    return FourierSeries(
        amplitudes, phases, frequencies,
        get_names(ts), get_longitudes(ts), get_latitudes(ts),
        get_quantity(ts), get_source(ts) * " | fft",
    )
end

# ── inverse FFT ───────────────────────────────────────────────────────────────

"""
    ifft_series(fs, times) -> TimeSeries

Reconstruct a TimeSeries from a FourierSeries.

`times` defines the output time grid and must have length `N` satisfying
`N ÷ 2 + 1 == length(get_frequencies(fs))`.  For a round-trip, pass the
original `get_times(ts)`.

The source string of the result ends with `" | ifft"`.
"""
function ifft_series(fs::AbstractFourierSeries, times::Vector{DateTime})::TimeSeries
    amps   = get_amplitudes(fs)      # Float32 [locations × frequencies]
    phases_deg = get_phases(fs)      # Float32 [locations × frequencies]

    N       = length(times)
    n_freqs = N ÷ 2 + 1

    if size(amps, 2) != n_freqs
        error("FourierSeries has $(size(amps, 2)) frequency bins but " *
              "length(times) = $N requires $n_freqs bins.")
    end

    n_locs = size(amps, 1)
    vals   = Matrix{Float32}(undef, n_locs, N)

    for i in 1:n_locs
        a   = Float64.(amps[i, :])
        phi = deg2rad.(Float64.(phases_deg[i, :]))

        # Rebuild the unnormalised complex spectrum that rfft would have produced
        X       = @. (N / 2.0) * a * exp(im * phi)
        X[1]   *= 2.0                        # DC: undo the 1/N factor
        if N % 2 == 0
            X[end] *= 2.0                    # Nyquist: undo the 1/N factor
        end

        vals[i, :] = Float32.(irfft(X, N))
    end

    return TimeSeries(
        vals, times,
        get_names(fs), get_longitudes(fs), get_latitudes(fs),
        get_quantity(fs), get_source(fs) * " | ifft",
    )
end
