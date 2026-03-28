# fourier_series.jl
#
# Data type for Fourier spectra across multiple locations.
# Analogous to constituents.jl, with the following mapping:
#
#   TidalConstituents                    FourierSeries
#   ────────────────────────────────     ────────────────────────────────────────
#   amplitudes[locations × consts]    →  amplitudes[locations × frequencies]
#   phases[locations × consts]        →  phases[locations × frequencies]
#   constituent_names::Vector{String} →  frequencies::Vector{Float64}  [Hz]
#   names, lons, lats, qty, src       →  names, lons, lats, qty, src  (identical)
#
# All locations share the same frequency axis: either all have a given bin or none do.

# Extend MultiTimeSeries generic functions so that FourierSeries shares the same
# interface as TimeSeries / TidalConstituents (single-dispatch, no ambiguity).
import MultiTimeSeries: get_names, get_longitudes, get_latitudes, get_quantity, get_source
import MultiTimeSeries: select_locations_by_ids, select_location_by_id
import MultiTimeSeries: select_locations_by_names, select_location_by_name

struct FourierSeries <: AbstractFourierSeries
    amplitudes::Matrix{Float32}   # one-sided amplitude spectrum [locations × frequencies]
    phases::Matrix{Float32}       # phase spectrum in degrees   [locations × frequencies]
    frequencies::Vector{Float64}  # frequency axis in Hz
    names::Vector{String}         # station names
    longitudes::Vector{Float64}   # longitudes of stations
    latitudes::Vector{Float64}    # latitudes of stations
    quantity::String              # physical quantity (e.g. "water level")
    source::String                # source of the data or analysis
end

"""
    FourierSeries(fs::AbstractFourierSeries)

Copy constructor from any AbstractFourierSeries object.
"""
function FourierSeries(fs::AbstractFourierSeries)
    return FourierSeries(
        get_amplitudes(fs), get_phases(fs), get_frequencies(fs),
        get_names(fs), get_longitudes(fs), get_latitudes(fs),
        get_quantity(fs), get_source(fs),
    )
end

#
# Getters
#
get_amplitudes(fs::FourierSeries)  = fs.amplitudes
get_phases(fs::FourierSeries)      = fs.phases
get_frequencies(fs::FourierSeries) = fs.frequencies
get_names(fs::FourierSeries)       = fs.names
get_longitudes(fs::FourierSeries)  = fs.longitudes
get_latitudes(fs::FourierSeries)   = fs.latitudes
get_quantity(fs::FourierSeries)    = fs.quantity
get_source(fs::FourierSeries)      = fs.source

#
# Selection by location — analogous to select_locations_by_ids / select_location_by_name in constituents.jl
#
function select_locations_by_ids(fs::AbstractFourierSeries, location_indices::Vector{T} where T<:Integer)
    return FourierSeries(
        get_amplitudes(fs)[location_indices, :],
        get_phases(fs)[location_indices, :],
        get_frequencies(fs),
        get_names(fs)[location_indices],
        get_longitudes(fs)[location_indices],
        get_latitudes(fs)[location_indices],
        get_quantity(fs),
        get_source(fs),
    )
end

function select_location_by_id(fs::AbstractFourierSeries, location_index::Integer)
    return select_locations_by_ids(fs, [location_index])
end

function select_locations_by_names(fs::AbstractFourierSeries, location_names::Vector{String})
    all_names = get_names(fs)
    indices = [findfirst(==(n), all_names) for n in location_names]
    missing_names = location_names[findall(isnothing, indices)]
    if !isempty(missing_names)
        error("Locations not found in FourierSeries: " * join(missing_names, ", "))
    end
    return select_locations_by_ids(fs, Int.(indices))
end

function select_location_by_name(fs::AbstractFourierSeries, location_name::String)
    return select_locations_by_names(fs, [location_name])
end

#
# Selection by frequency — analogous to select_constituents_by_names in constituents.jl
#
"""
    select_frequencies_by_indices(fs, indices) -> FourierSeries

Return a new FourierSeries containing only the frequency bins at `indices`.
All locations are kept; the shared frequency axis is reduced to `indices`.
"""
function select_frequencies_by_indices(fs::AbstractFourierSeries, indices::Vector{<:Integer})
    return FourierSeries(
        get_amplitudes(fs)[:, indices],
        get_phases(fs)[:, indices],
        get_frequencies(fs)[indices],
        get_names(fs),
        get_longitudes(fs),
        get_latitudes(fs),
        get_quantity(fs),
        get_source(fs),
    )
end

#
# Merge by locations — analogous to merge_by_locations in constituents.jl
# Each element of the vector must share the same frequency axis.
#
"""
    merge_by_locations(fourier_vector) -> FourierSeries

Merge a vector of FourierSeries (one location each) into a single object.
All elements must share the same frequency axis, quantity, and source.
"""
function merge_by_locations(fourier_vector::Vector{<:AbstractFourierSeries})
    if isempty(fourier_vector)
        error("Cannot merge an empty vector of FourierSeries.")
    end
    ref         = fourier_vector[1]
    quantity    = get_quantity(ref)
    source      = get_source(ref)
    frequencies = get_frequencies(ref)
    n_freqs     = length(frequencies)
    n_locations = length(fourier_vector)
    names       = Vector{String}()
    longitudes  = Vector{Float64}()
    latitudes   = Vector{Float64}()
    amplitudes  = zeros(Float32, n_locations, n_freqs)
    phases      = zeros(Float32, n_locations, n_freqs)
    for (i, fs) in enumerate(fourier_vector)
        if get_frequencies(fs) != frequencies
            error("Cannot merge FourierSeries with different frequency axes.")
        end
        if get_quantity(fs) != quantity
            error("Cannot merge FourierSeries with different quantities: $(get_quantity(fs)) != $quantity")
        end
        if get_source(fs) != source
            error("Cannot merge FourierSeries with different sources: $(get_source(fs)) != $source")
        end
        push!(names,      get_names(fs)[1])
        push!(longitudes, get_longitudes(fs)[1])
        push!(latitudes,  get_latitudes(fs)[1])
        amplitudes[i, :] = get_amplitudes(fs)[1, :]
        phases[i, :]     = get_phases(fs)[1, :]
    end
    return FourierSeries(amplitudes, phases, frequencies,
                         names, longitudes, latitudes, quantity, source)
end

#
# Show — analogous to Base.show in constituents.jl
#
function Base.show(io::IO, fs::AbstractFourierSeries)
    freqs = get_frequencies(fs)
    println(io, "AbstractFourierSeries:")
    println(io, "   Quantity:         ", get_quantity(fs))
    println(io, "   Source:           ", get_source(fs))
    println(io, "   Locations:        ", join(get_names(fs), ", "))
    println(io, "   N locations:      ", length(get_names(fs)))
    println(io, "   N frequencies:    ", length(freqs))
    println(io, "   Freq range [Hz]:  ", isempty(freqs) ? "—" : "$(freqs[1]) … $(freqs[end])")
    println(io, "   Data shape:       ", size(get_amplitudes(fs)), " [locations × frequencies]")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", fs::AbstractFourierSeries)
    println(io, "AbstractFourierSeries: $(get_quantity(fs)) from $(get_source(fs)), " *
                "$(length(get_names(fs))) locations, $(length(get_frequencies(fs))) frequencies.")
    return nothing
end

function Base.show(io::IO, fs::FourierSeries)
    freqs = get_frequencies(fs)
    println(io, "FourierSeries:")
    println(io, "   Quantity:         ", get_quantity(fs))
    println(io, "   Source:           ", get_source(fs))
    println(io, "   Locations:        ", join(get_names(fs), ", "))
    println(io, "   N locations:      ", length(get_names(fs)))
    println(io, "   N frequencies:    ", length(freqs))
    println(io, "   Freq range [Hz]:  ", isempty(freqs) ? "—" : "$(freqs[1]) … $(freqs[end])")
    println(io, "   Data shape:       ", size(get_amplitudes(fs)), " [locations × frequencies]")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", fs::FourierSeries)
    println(io, "FourierSeries: $(get_quantity(fs)) from $(get_source(fs)), " *
                "$(length(get_names(fs))) locations, $(length(get_frequencies(fs))) frequencies.")
    return nothing
end
