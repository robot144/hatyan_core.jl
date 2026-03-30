# statistics.jl
#
# Per-station validation statistics comparing a model TimeSeries against
# measurements.
#
# Output columns
# ──────────────
# location_id   : 1-based location index
# location_name : station name (from measurements)
# signal_rmse   : √(mean(obs²)) — RMS amplitude of the observed signal
# bias          : mean(model − obs)
# rmse          : √(mean((model − obs)²))
# mae           : mean(|model − obs|)
# max_error     : max(model − obs)
# min_error     : min(model − obs)
#
# NaN handling
# ────────────
# Time steps where either observations or model contain NaN are excluded
# from all statistics at that location.
#
# Time alignment
# ──────────────
# Both TimeSeries must share the same time axis.  An error is thrown if the
# lengths differ or if any corresponding timestamps do not match.  Use
# select_timespan / merge_by_times to align the series before calling this
# function.

using Statistics: mean
using DataFrames

"""
    compute_statistics(meas, model) -> DataFrame

Compute per-station validation statistics comparing `model` against `meas`.

Both `meas` and `model` must share the same time axis and have the same
locations in matching order.  Time steps where either series contains `NaN`
are excluded from all metrics at that location.

# Output columns

| Column        | Description                                            |
|:--------------|:-------------------------------------------------------|
| `location_id`   | 1-based index in the input TimeSeries                |
| `location_name` | Station name from `meas`                             |
| `n_values`      | Number of valid (non-NaN) time steps used            |
| `signal_rmse`   | √(mean(obs²)) — RMS amplitude of the observed signal |
| `bias`          | mean(model − obs)                                    |
| `rmse`          | √(mean((model − obs)²))                              |
| `mae`           | mean(\\|model − obs\\|)                              |
| `max_error`     | max(model − obs)                                     |
| `min_error`     | min(model − obs)                                     |

# Example
```julia
stats = compute_statistics(obs1990, pred1990)
# → DataFrame with one row per station
```
"""
function compute_statistics(
    meas::AbstractTimeSeries,
    model::AbstractTimeSeries,
)::DataFrame
    _check_time_alignment(meas, model)
    _check_location_alignment(meas, model)

    obs_vals   = get_values(meas)    # Float32 [locations × times]
    model_vals = get_values(model)   # Float32 [locations × times]
    loc_names  = get_names(meas)

    n_locs = size(obs_vals, 1)

    location_ids   = Vector{Int}(undef, n_locs)
    location_names = Vector{String}(undef, n_locs)
    n_values       = Vector{Int}(undef, n_locs)
    signal_rmse    = Vector{Float64}(undef, n_locs)
    bias           = Vector{Float64}(undef, n_locs)
    rmse           = Vector{Float64}(undef, n_locs)
    mae            = Vector{Float64}(undef, n_locs)
    max_error      = Vector{Float64}(undef, n_locs)
    min_error      = Vector{Float64}(undef, n_locs)

    for i in 1:n_locs
        obs = Float64.(obs_vals[i, :])
        mod = Float64.(model_vals[i, :])

        # Drop time steps where either series is NaN
        valid = .!isnan.(obs) .& .!isnan.(mod)
        obs   = obs[valid]
        mod   = mod[valid]

        if isempty(obs)
            error("Location $(loc_names[i]): no valid (non-NaN) overlapping " *
                  "time steps found.")
        end

        err = mod .- obs

        location_ids[i]   = i
        location_names[i] = loc_names[i]
        n_values[i]       = length(obs)
        signal_rmse[i]    = sqrt(mean(obs .^ 2))
        bias[i]           = mean(err)
        rmse[i]           = sqrt(mean(err .^ 2))
        mae[i]            = mean(abs.(err))
        max_error[i]      = maximum(err)
        min_error[i]      = minimum(err)
    end

    return DataFrame(
        location_id   = location_ids,
        location_name = location_names,
        n_values      = n_values,
        signal_rmse   = signal_rmse,
        bias          = bias,
        rmse          = rmse,
        mae           = mae,
        max_error     = max_error,
        min_error     = min_error,
    )
end

# ── internal helpers ──────────────────────────────────────────────────────────

function _check_time_alignment(a::AbstractTimeSeries, b::AbstractTimeSeries)
    t_a = get_times(a)
    t_b = get_times(b)
    if length(t_a) != length(t_b)
        error("Time axes have different lengths: $(length(t_a)) vs $(length(t_b)). " *
              "Align the TimeSeries before calling compute_statistics.")
    end
    mismatch = findfirst(t_a .!= t_b)
    if mismatch !== nothing
        error("Time axes do not match at index $mismatch " *
              "($(t_a[mismatch]) vs $(t_b[mismatch])). " *
              "Align the TimeSeries before calling compute_statistics.")
    end
end

function _check_location_alignment(a::AbstractTimeSeries, b::AbstractTimeSeries)
    n_a = size(get_values(a), 1)
    n_b = size(get_values(b), 1)
    if n_a != n_b
        error("Location counts differ: measurements has $n_a, model has $n_b.")
    end
    names_a = get_names(a)
    names_b = get_names(b)
    if names_a != names_b
        error("Location names do not match.\n" *
              "  measurements : $names_a\n" *
              "  model        : $names_b")
    end
end
