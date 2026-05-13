# plotting.jl
#
# Plots.jl methods for tidal-specific data types (TidalConstituents, FourierSeries).
# The generic AbstractTimeSeries plot method lives in MultiTimeSeries.jl.
#
# Conventions
# ───────────
# • One series per plot; location_index selects which location to show.
# • TidalConstituents uses a (2 × 1) layout: amplitude bar chart (top) and
#   phase scatter (bottom).

import Plots
import Plots: mm

# ── TidalConstituents ──────────────────────────────────────────────────────────

"""
    plot(tc::AbstractTidalConstituents; location_index=1, kwargs...) -> Plots.Plot

Two-panel plot: amplitude bar chart (top) and phase scatter (bottom).
Constituent names appear on the shared x-axis; x-tick labels are rotated 90°.

# Keyword arguments (in addition to all standard Plots.jl kwargs)
- `location_index`: which location to plot (default `1`).
- `max_constituents`: show only the top N constituents ranked by amplitude
  (default `20`).  Set to `typemax(Int)` to show all.
"""
function Plots.plot(tc::AbstractTidalConstituents;
                    location_index::Integer = 1,
                    label = nothing,
                    max_constituents::Integer = 20,
                    kwargs...)
    cnames = get_constituent_names(tc)
    amps   = get_amplitudes(tc)
    phases = get_phases(tc)
    lnames = get_names(tc)
    qty    = get_quantity(tc)

    if location_index ∉ eachindex(lnames)
        error("location_index $location_index is out of range " *
              "($(length(lnames)) location(s) available).")
    end

    # Select the top max_constituents by amplitude, then sort by frequency.
    n_show      = min(length(cnames), max_constituents)
    top_indices = sortperm(vec(amps[location_index, :]), rev=true)[1:n_show]
    dummy_date  = [DateTime(2000, 1, 1)]
    top_freqs   = get_schureman_freqs(cnames[top_indices], dummy_date)
    order       = top_indices[sortperm(vec(top_freqs))]
    cnames      = cnames[order]
    amps_loc    = amps[location_index, order]
    phases_loc  = phases[location_index, order]
    x           = 1:n_show
    xtick_args  = (collect(x), cnames)
    lbl         = isnothing(label) ? lnames[location_index] : label

    p = Plots.plot(;
        layout        = (2, 1),
        bottom_margin = 8mm,
        left_margin   = 5mm,
        kwargs...,
    )

    Plots.bar!(p[1], x, amps_loc;
        label     = lbl,
        ylabel    = "$qty amplitude (m)",
        title     = lnames[location_index],
        xticks    = xtick_args,
        xrotation = 90,
    )

    Plots.scatter!(p[2], x, phases_loc;
        label     = lbl,
        ylabel    = "Phase (°)",
        ylims     = (0, 360),
        xticks    = xtick_args,
        xrotation = 90,
        markersize = 3,
        markerstrokewidth = 0,
    )

    return p
end

"""
    plot(tc_ref, tc_comp; location_index=1, label_ref=nothing, label_comp=nothing,
         max_constituents=20, kwargs...) -> Plots.Plot

Two-panel comparison plot of two `TidalConstituents` objects (e.g. observed vs
predicted).  Constituents are selected and sorted the same way as the single-TC
`plot`: top `max_constituents` by `tc_ref` amplitude, then re-sorted by frequency.

- Panel 1 (top): grouped amplitude bar chart — `tc_ref` and `tc_comp` side by side.
- Panel 2 (bottom): phase scatter — `tc_ref` circles, `tc_comp` diamonds.

`tc_ref` and `tc_comp` must share the same constituent list.

# Keyword arguments
- `location_index`: which location to plot (default `1`).
- `label_ref`: legend label for `tc_ref` (default: source string of `tc_ref`).
- `label_comp`: legend label for `tc_comp` (default: source string of `tc_comp`).
- `max_constituents`: number of constituents to show (default `20`).
"""
function Plots.plot(tc_ref::AbstractTidalConstituents,
                    tc_comp::AbstractTidalConstituents;
                    location_index::Integer = 1,
                    label_ref               = nothing,
                    label_comp              = nothing,
                    max_constituents::Integer = 20,
                    kwargs...)

    cnames  = get_constituent_names(tc_ref)
    lnames  = get_names(tc_ref)
    qty     = get_quantity(tc_ref)

    if location_index ∉ eachindex(lnames)
        error("location_index $location_index is out of range " *
              "($(length(lnames)) location(s) available).")
    end
    length(cnames) == length(get_constituent_names(tc_comp)) ||
        error("tc_ref and tc_comp have different numbers of constituents.")

    amps_ref  = get_amplitudes(tc_ref)
    amps_comp = get_amplitudes(tc_comp)
    phi_ref   = get_phases(tc_ref)
    phi_comp  = get_phases(tc_comp)

    # Select top N from tc_ref by amplitude, then re-sort by frequency.
    n_show      = min(length(cnames), max_constituents)
    top_indices = sortperm(vec(amps_ref[location_index, :]), rev=true)[1:n_show]
    dummy_date  = [DateTime(2000, 1, 1)]
    top_freqs   = get_schureman_freqs(cnames[top_indices], dummy_date)
    order       = top_indices[sortperm(vec(top_freqs))]

    cnames_show = cnames[order]
    x           = 1:n_show
    xtick_args  = (collect(x), cnames_show)
    lbl_ref     = isnothing(label_ref)  ? get_source(tc_ref)  : label_ref
    lbl_comp    = isnothing(label_comp) ? get_source(tc_comp) : label_comp

    # Grouped bar chart: (n_show × 2) matrix → Plots groups automatically.
    amp_matrix = hcat(amps_ref[location_index, order], amps_comp[location_index, order])

    p = Plots.plot(;
        layout        = (2, 1),
        bottom_margin = 8mm,
        left_margin   = 5mm,
        kwargs...,
    )

    Plots.bar!(p[1], x, amp_matrix;
        label     = [lbl_ref lbl_comp],
        ylabel    = "$qty amplitude (m)",
        title     = lnames[location_index],
        xticks    = xtick_args,
        xrotation = 90,
        alpha     = 0.75,
    )

    Plots.scatter!(p[2], x, phi_ref[location_index, order];
        label             = lbl_ref,
        ylabel            = "Phase (°)",
        ylims             = (0, 360),
        xticks            = xtick_args,
        xrotation         = 90,
        markersize        = 4,
        markerstrokewidth = 0,
    )
    Plots.scatter!(p[2], x, phi_comp[location_index, order];
        label             = lbl_comp,
        markersize        = 4,
        markerstrokewidth = 0,
        marker            = :diamond,
    )

    return p
end

# ── FourierSeries ──────────────────────────────────────────────────────────────

"""
    plot(fs::AbstractFourierSeries; location_index=1, kwargs...) -> Plots.Plot

Amplitude spectrum: frequency on the x-axis, amplitude on the y-axis.

# Keyword arguments (in addition to all standard Plots.jl kwargs)
- `location_index`: which location to plot (default `1`).
- `xscale`: x-axis scale passed to Plots.jl, e.g. `:log10` (default `:identity`).
- `freq_unit`: `"cpd"` (cycles per day, default) or `"Hz"` or `"cph"`
  (cycles per hour).  Frequencies are converted before plotting.
- `freq_max`: upper limit of the x-axis in the chosen `freq_unit` (default
  `20.0` cpd — covers all major tidal bands up to the 8th diurnal harmonic).
  Set to `Inf` to show the full spectrum up to the Nyquist.
- `yunit`: string appended to the y-axis label (default `""`).
"""
function Plots.plot(fs::AbstractFourierSeries;
                    location_index::Integer = 1,
                    label = nothing,
                    xscale::Symbol    = :identity,
                    freq_unit::String = "cpd",
                    freq_max::Real    = 20.0,
                    yunit::String     = "",
                    kwargs...)
    freqs  = get_frequencies(fs)    # always in Hz
    amps   = get_amplitudes(fs)
    names  = get_names(fs)
    qty    = get_quantity(fs)

    if location_index ∉ eachindex(names)
        error("location_index $location_index is out of range " *
              "($(length(names)) location(s) available).")
    end

    # Convert frequency axis
    scale_factor, xlabel = if freq_unit == "cph"
        3600.0, "Frequency (cycles/hour)"
    elseif freq_unit == "cpd"
        86400.0, "Frequency (cycles/day)"
    else
        1.0, "Frequency (Hz)"
    end
    freqs_plot = freqs .* scale_factor
    amps_loc   = amps[location_index, :]

    ylabel = isempty(yunit) ? "$qty amplitude" : "$qty amplitude ($yunit)"
    lbl    = isnothing(label) ? names[location_index] : label

    # Drop the DC bin (f=0) when a log x-scale is requested to avoid log(0) warnings.
    if xscale === :log10 || xscale === :ln || xscale === :log2
        mask       = freqs_plot .> 0
        freqs_plot = freqs_plot[mask]
        amps_loc   = amps_loc[mask]
    end

    return Plots.plot(freqs_plot, amps_loc;
        label         = lbl,
        xlabel        = xlabel,
        ylabel        = ylabel,
        title         = names[location_index],
        xscale        = xscale,
        xlims         = (0.0, Float64(freq_max)),
        bottom_margin = 5mm,
        left_margin   = 5mm,
        kwargs...,
    )
end
