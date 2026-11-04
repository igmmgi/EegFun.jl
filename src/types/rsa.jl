"""
    NoiseCeiling <: EegFunData

Leave-one-participant-out RSA noise ceiling estimates across time (Nili et al., 2014).

# Fields
- `lower_bound::Vector{Float64}`: Conservative lower bound at each time point
- `upper_bound::Vector{Float64}`: Liberal upper bound at each time point
- `n_participants::Int`: Number of participants included
"""
struct NoiseCeiling <: EegFunData
    lower_bound::Vector{Float64}
    upper_bound::Vector{Float64}
    n_participants::Int
end

Base.copy(nc::NoiseCeiling)::NoiseCeiling = NoiseCeiling(copy(nc.lower_bound), copy(nc.upper_bound), nc.n_participants)

function _bound_range_str(v::Vector{Float64}; include_mean::Bool = false)
    isempty(v) && return "N/A"
    lo, hi = round.(extrema(v); digits = 3)
    return include_mean ? "$(round(mean(v); digits = 3)) (range: $lo-$hi)" : "$lo to $hi"
end

function Base.show(io::IO, nc::NoiseCeiling)
    println(io, "NoiseCeiling:")
    println(io, "  Participants: $(nc.n_participants)")
    println(io, "  Time points: $(length(nc.lower_bound))")
    println(io, "  Lower bound: $(_bound_range_str(nc.lower_bound))")
    println(io, "  Upper bound: $(_bound_range_str(nc.upper_bound))")
end

"""
    RsaData <: EegFunData

Time-resolved Representational Similarity Analysis (RSA) results.

# Fields
- `file::String`: Source filename
- `condition_names::Vector{String}`: Analyzed condition labels
- `times::Vector{Float64}`: Time points in seconds
- `rdm::Array{Float64,3}`: Representational dissimilarity matrices `[time × condition × condition]`
- `dissimilarity_measure::Symbol`: Dissimilarity metric (`:correlation`, `:euclidean`, `:mahalanobis`, `:spearman`)
- `channels::Vector{Symbol}`: Included electrode labels
- `layout::Layout`: Electrode layout
- `sample_rate::Int64`: Sample rate in Hz
- `model_correlations::Union{Array{Float64,2},Nothing}`: Model RDM correlations `[time × model]`
- `model_names::Union{Vector{String},Nothing}`: Model names
- `p_values::Union{Array{Float64,2},Nothing}`: Model correlation p-values `[time × model]`
- `noise_ceiling::Union{NoiseCeiling,Nothing}`: Group noise ceiling estimates
- `analysis_info::AnalysisInfo`: Preprocessing metadata
"""
mutable struct RsaData <: EegFunData
    file::String
    condition_names::Vector{String}
    times::Vector{Float64}
    rdm::Array{Float64,3}
    dissimilarity_measure::Symbol
    channels::Vector{Symbol}
    layout::Layout
    sample_rate::Int64
    model_correlations::Union{Array{Float64,2},Nothing}
    model_names::Union{Vector{String},Nothing}
    p_values::Union{Array{Float64,2},Nothing}
    noise_ceiling::Union{NoiseCeiling,Nothing}
    analysis_info::AnalysisInfo
end

function RsaData(
    file::String,
    condition_names::Vector{String},
    times::Vector{Float64},
    rdm::Array{Float64,3},
    dissimilarity_measure::Symbol,
    channels::Vector{Symbol},
    layout::Layout,
    sample_rate::Int64;
    model_correlations::Union{Array{Float64,2},Nothing} = nothing,
    model_names::Union{Vector{String},Nothing} = nothing,
    p_values::Union{Array{Float64,2},Nothing} = nothing,
    noise_ceiling::Union{NoiseCeiling,Nothing} = nothing,
    analysis_info::AnalysisInfo = AnalysisInfo(),
)
    return RsaData(
        file,
        condition_names,
        times,
        rdm,
        dissimilarity_measure,
        channels,
        copy(layout),
        sample_rate,
        model_correlations,
        model_names,
        p_values,
        noise_ceiling,
        copy(analysis_info),
    )
end

function Base.copy(rsa::RsaData)::RsaData
    return RsaData(
        rsa.file,
        copy(rsa.condition_names),
        copy(rsa.times),
        copy(rsa.rdm),
        rsa.dissimilarity_measure,
        copy(rsa.channels),
        copy(rsa.layout),
        rsa.sample_rate,
        isnothing(rsa.model_correlations) ? nothing : copy(rsa.model_correlations),
        isnothing(rsa.model_names) ? nothing : copy(rsa.model_names),
        isnothing(rsa.p_values) ? nothing : copy(rsa.p_values),
        isnothing(rsa.noise_ceiling) ? nothing : copy(rsa.noise_ceiling),
        copy(rsa.analysis_info),
    )
end

function Base.show(io::IO, rsa::RsaData)
    time_str = isempty(rsa.times) ? "N/A" : "$(first(rsa.times)) to $(last(rsa.times)) s"
    println(io, "RsaData:")
    println(io, "  File: $(rsa.file)")
    println(io, "  Conditions: $(join(rsa.condition_names, ", ")) (n=$(length(rsa.condition_names)))")
    println(io, "  Dissimilarity measure: $(rsa.dissimilarity_measure)")
    println(io, "  Time points: $(length(rsa.times)) ($time_str)")
    println(io, "  Channels: $(length(rsa.channels))")

    if !isnothing(rsa.model_correlations) && !isnothing(rsa.model_names) && size(rsa.model_correlations, 1) > 0
        println(io, "  Model comparisons: $(length(rsa.model_names))")
        for (idx, model_name) in enumerate(rsa.model_names)
            idx > size(rsa.model_correlations, 2) && break
            col = @view rsa.model_correlations[:, idx]
            max_idx = argmax(col)
            max_r = round(col[max_idx]; digits = 3)
            if max_idx <= length(rsa.times)
                println(io, "    - $model_name: max r=$max_r at $(round(rsa.times[max_idx]; digits = 3)) s")
            else
                println(io, "    - $model_name: max r=$max_r")
            end
        end
    end

    if !isnothing(rsa.noise_ceiling)
        nc = rsa.noise_ceiling
        println(io, "  Noise ceiling ($(nc.n_participants) participants):")
        println(io, "    - Lower: $(_bound_range_str(nc.lower_bound; include_mean = true))")
        println(io, "    - Upper: $(_bound_range_str(nc.upper_bound; include_mean = true))")
    end
end
