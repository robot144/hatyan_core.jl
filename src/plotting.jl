# plotting.jl
#
# Plots.jl methods for the three main data types.
#
# Each function returns a `Plots.Plot` and accepts any keyword arguments that
# Plots.jl itself accepts (title, size, dpi, …).  Additional domain-specific
# keyword arguments are documented per function.
#
# Conventions
# ───────────
# • One series per location; location names become the legend labels.
# • Multi-location TidalConstituents / FourierSeries use layout sub-panels so
#   every location gets its own axes and is easy to read.
# • TimeSeries locations are overlaid on a single axes (standard for comparison).

import Plots
import Plots: mm

# ── TimeSeries ─────────────────────────────────────────────────────────────────

"""
    plot(ts::AbstractTimeSeries; location_index=1, kwargs...) -> Plots.Plot

Line plot of one location in `ts` versus time.

# Keyword arguments (in addition to all standard Plots.jl kwargs)
- `location_index`: which location to plot (default `1`).  Pass `nothing` to
  overlay all locations on the same axes.
- `yunit`: string appended to the y-axis label (default `""`).
"""
function Plots.plot(ts::AbstractTimeSeries;
                    location_index::Union{Integer, Nothing} = 1,
                    yunit::String = "",
                    kwargs...)
    times  = get_times(ts)
    vals   = get_values(ts)
    names  = get_names(ts)
    qty    = get_quantity(ts)

    indices = isnothing(location_index) ? eachindex(names) : (location_index:location_index)

    if !isnothing(location_index) && location_index ∉ eachindex(names)
        error("location_index $location_index is out of range " *
              "($(length(names)) location(s) available).")
    end

    ylabel = isempty(yunit) ? qty : "$qty ($yunit)"

    # Plot the first series in the initial plot() call so that axis attributes
    # (xlabel in particular) are guaranteed to be applied with data present.
    first_i = first(indices)
    p = Plots.plot(times, vals[first_i, :];
        label         = names[first_i],
        xlabel        = "Time",
        ylabel        = ylabel,
        title         = get_source(ts),
        legend        = :outertopright,
        bottom_margin = 5mm,
        left_margin   = 5mm,
        kwargs...,
    )
    for i in Iterators.drop(indices, 1)
        Plots.plot!(p, times, vals[i, :]; label = names[i])
    end
    return p
end

# ── TidalConstituents ──────────────────────────────────────────────────────────

"""
    plot(tc::AbstractTidalConstituents; location_index=1, kwargs...) -> Plots.Plot

Two-row panel per location: amplitude bar chart (top) and phase scatter (bottom).
Constituent names appear on the shared x-axis; x-tick labels are rotated 90°.

# Keyword arguments (in addition to all standard Plots.jl kwargs)
- `location_index`: which location to plot (default `1`).  Pass `nothing` to
  show all locations side-by-side in a `(2 × N)` grid.
- `max_constituents`: show only the top N constituents ranked by mean amplitude
  across the plotted locations (default `20`).  Set to `typemax(Int)` to show all.
"""
function Plots.plot(tc::AbstractTidalConstituents;
                    location_index::Union{Integer, Nothing} = 1,
                    max_constituents::Integer = 20,
                    kwargs...)
    cnames = get_constituent_names(tc)
    amps   = get_amplitudes(tc)
    phases = get_phases(tc)
    lnames = get_names(tc)
    qty    = get_quantity(tc)

    if !isnothing(location_index) && location_index ∉ eachindex(lnames)
        error("location_index $location_index is out of range " *
              "($(length(lnames)) location(s) available).")
    end

    indices = isnothing(location_index) ? eachindex(lnames) : (location_index:location_index)
    n_plot  = length(indices)

    # Select the top max_constituents by mean amplitude across the plotted locations,
    # then sort the selection by constituent frequency (ascending) for display.
    n_show      = min(length(cnames), max_constituents)
    mean_amp    = vec(sum(amps[collect(indices), :], dims=1)) ./ length(indices)
    top_indices = sortperm(mean_amp, rev=true)[1:n_show]
    dummy_date  = [DateTime(2000, 1, 1)]
    top_freqs   = get_schureman_freqs(cnames[top_indices], dummy_date)
    order       = top_indices[sortperm(vec(top_freqs))]
    cnames     = cnames[order]
    amps       = amps[:,   order]
    phases     = phases[:, order]
    x          = 1:n_show
    xtick_args = (collect(x), cnames)

    p = Plots.plot(;
        layout        = (2, n_plot),
        bottom_margin = 8mm,
        left_margin   = 5mm,
        kwargs...,
    )

    for (col, j) in enumerate(indices)
        amp_idx   = col           # top row: column col
        phase_idx = n_plot + col  # bottom row: column col

        Plots.bar!(p[amp_idx], x, amps[j, :];
            label     = lnames[j],
            ylabel    = "$qty amplitude (m)",
            title     = lnames[j],
            xticks    = xtick_args,
            xrotation = 90,
        )

        Plots.scatter!(p[phase_idx], x, phases[j, :];
            label     = lnames[j],
            ylabel    = "Phase (°)",
            ylims     = (0, 360),
            xticks    = xtick_args,
            xrotation = 90,
            markersize = 3,
            markerstrokewidth = 0,
        )
    end

    return p
end

# ── FourierSeries ──────────────────────────────────────────────────────────────

"""
    plot(fs::AbstractFourierSeries; location_index=1, kwargs...) -> Plots.Plot

Amplitude spectrum: frequency on the x-axis, amplitude on the y-axis.

# Keyword arguments (in addition to all standard Plots.jl kwargs)
- `location_index`: which location to plot (default `1`).  Pass `nothing` to
  show all locations, each in its own sub-panel.
- `xscale`: x-axis scale passed to Plots.jl, e.g. `:log10` (default `:identity`).
- `freq_unit`: `"cpd"` (cycles per day, default) or `"Hz"` or `"cph"`
  (cycles per hour).  Frequencies are converted before plotting.
- `freq_max`: upper limit of the x-axis in the chosen `freq_unit` (default
  `20.0` cpd — covers all major tidal bands up to the 8th diurnal harmonic).
  Set to `Inf` to show the full spectrum up to the Nyquist.
- `yunit`: string appended to the y-axis label (default `""`).
"""
function Plots.plot(fs::AbstractFourierSeries;
                    location_index::Union{Integer, Nothing} = 1,
                    xscale::Symbol    = :identity,
                    freq_unit::String = "cpd",
                    freq_max::Real    = 20.0,
                    yunit::String     = "",
                    kwargs...)
    freqs  = get_frequencies(fs)    # always in Hz
    amps   = get_amplitudes(fs)
    names  = get_names(fs)
    qty    = get_quantity(fs)

    if !isnothing(location_index) && location_index ∉ eachindex(names)
        error("location_index $location_index is out of range " *
              "($(length(names)) location(s) available).")
    end

    indices = isnothing(location_index) ? eachindex(names) : (location_index:location_index)
    n_plot  = length(indices)

    # Convert frequency axis
    scale_factor, xlabel = if freq_unit == "cph"
        3600.0, "Frequency (cycles/hour)"
    elseif freq_unit == "cpd"
        86400.0, "Frequency (cycles/day)"
    else
        1.0, "Frequency (Hz)"
    end
    freqs_plot = freqs .* scale_factor

    ylabel = isempty(yunit) ? "$qty amplitude" : "$qty amplitude ($yunit)"

    # Drop the DC bin (f=0) when a log x-scale is requested to avoid log(0) warnings.
    if xscale === :log10 || xscale === :ln || xscale === :log2
        mask       = freqs_plot .> 0
        freqs_plot = freqs_plot[mask]
        amps       = amps[:, mask]
    end

    p = Plots.plot(;
        layout        = (n_plot, 1),
        bottom_margin = 5mm,
        left_margin   = 5mm,
        kwargs...,
    )

    xlims = (0.0, Float64(freq_max))

    for (panel, i) in enumerate(indices)
        Plots.plot!(p[panel], freqs_plot, amps[i, :];
            label     = get(kwargs, :label, names[i]),
            xlabel    = xlabel,
            ylabel    = ylabel,
            title     = names[i],
            xscale    = xscale,
            xlims     = xlims,
        )
    end

    return p
end
