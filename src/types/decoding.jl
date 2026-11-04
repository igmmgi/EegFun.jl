"""
    DecodingParameters

Cross-validation and multi-class coding parameters for MVPA/decoding.

# Fields
- `chance_level::Float64`: Theoretical chance accuracy (e.g., `0.5` for binary, `1 / n_classes` for multi-class)
- `n_iterations::Int64`: Number of cross-validation iterations
- `n_folds::Int64`: Number of cross-validation folds
- `class_coding::Symbol`: Multi-class coding scheme (`:one_vs_one`, `:one_vs_all`, `:binary`)
- `n_classes::Int64`: Number of decoded classes
"""
struct DecodingParameters
    chance_level::Float64
    n_iterations::Int64
    n_folds::Int64
    class_coding::Symbol
    n_classes::Int64
end

"""
    DecodedData <: EegFunData

Time-resolved MVPA/decoding results for a single participant.

# Fields
- `file::String`: Source filename
- `condition_names::Vector{String}`: Decoded condition names
- `times::Vector{Float64}`: Decoded time points in seconds
- `average_score::Vector{Float64}`: Mean classification accuracy at each time point
- `channels::Vector{Symbol}`: Channels included as features
- `parameters::DecodingParameters`: Cross-validation parameters
- `stderror::Union{Vector{Float64},Nothing}`: Standard error of accuracy across iterations/folds
- `confusion_matrix::Union{Array{Float64,3},Nothing}`: Confusion matrices `[time × true_class × predicted_class]`
- `raw_predictions::Union{Array{Float64,4},Nothing}`: Raw predictions `[iteration × fold × time × class]`
"""
mutable struct DecodedData <: EegFunData
    file::String
    condition_names::Vector{String}
    times::Vector{Float64}
    average_score::Vector{Float64}
    channels::Vector{Symbol}
    parameters::DecodingParameters
    stderror::Union{Vector{Float64},Nothing}
    confusion_matrix::Union{Array{Float64,3},Nothing}
    raw_predictions::Union{Array{Float64,4},Nothing}

    function DecodedData(
        file::String,
        condition_names::Vector{String},
        times::Vector{Float64},
        average_score::Vector{Float64},
        channels::Vector{Symbol},
        parameters::DecodingParameters;
        stderror::Union{Vector{Float64},Nothing} = nothing,
        confusion_matrix::Union{Array{Float64,3},Nothing} = nothing,
        raw_predictions::Union{Array{Float64,4},Nothing} = nothing,
    )
        return new(file, condition_names, times, average_score, channels, parameters, stderror, confusion_matrix, raw_predictions)
    end
end

function Base.copy(decoded::DecodedData)::DecodedData
    return DecodedData(
        decoded.file,
        copy(decoded.condition_names),
        copy(decoded.times),
        copy(decoded.average_score),
        copy(decoded.channels),
        decoded.parameters;
        stderror = isnothing(decoded.stderror) ? nothing : copy(decoded.stderror),
        confusion_matrix = isnothing(decoded.confusion_matrix) ? nothing : copy(decoded.confusion_matrix),
        raw_predictions = isnothing(decoded.raw_predictions) ? nothing : copy(decoded.raw_predictions),
    )
end

function Base.show(io::IO, decoded::DecodedData)
    time_str = isempty(decoded.times) ? "N/A" : "$(first(decoded.times)) to $(last(decoded.times)) s"
    println(io, "DecodedData:")
    println(io, "  File: $(decoded.file)")
    println(io, "  Conditions: $(join(decoded.condition_names, ", "))")
    println(io, "  Classes: $(decoded.parameters.n_classes)")
    println(io, "  Time range: $time_str")
    println(io, "  Iterations: $(decoded.parameters.n_iterations)")
    println(io, "  Cross-validation folds: $(decoded.parameters.n_folds)")
    println(io, "  Channels: $(length(decoded.channels))")
    if !isnothing(decoded.stderror) && !isempty(decoded.average_score) && !isempty(decoded.times)
        max_idx = argmax(decoded.average_score)
        max_acc = round(decoded.average_score[max_idx]; digits = 3)
        if max_idx <= length(decoded.times)
            println(io, "  Max accuracy: $max_acc at $(round(decoded.times[max_idx]; digits = 3)) s")
        else
            println(io, "  Max accuracy: $max_acc")
        end
    end
end
