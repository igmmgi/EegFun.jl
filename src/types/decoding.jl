"""
Types for MVPA/decoding analysis results.

This module defines types for storing multivariate pattern analysis (MVPA)
and decoding results, including classification accuracy over time and
confusion matrices.
"""

"""
    DecodedData

Stores MVPA/decoding analysis results for a single participant.

This type represents the results of multivariate pattern classification
analysis, including classification accuracy over time, confusion matrices,
and metadata about the analysis parameters.

# Fields
- `file::String`: Source filename
- `condition_names::Vector{String}`: Names of conditions/bins being decoded
- `times::Vector{Float64}`: Time points in seconds where decoding was performed
- `average_score::Vector{Float64}`: Average classification accuracy at each time point
- `stderror::Union{Vector{Float64}, Nothing}`: Standard error of accuracy (if computed)
- `chance_level::Float64`: Expected chance-level performance (e.g., 0.5 for binary, 1/n_classes for n-class)
- `n_classes::Int64`: Number of classes/conditions being decoded
- `n_iterations::Int64`: Number of iterations/permutations performed
- `n_folds::Int64`: Number of cross-validation folds
- `channels::Vector{Symbol}`: Channel names used in the analysis
- `layout::Layout`: Layout object containing electrode positioning information
- `sample_rate::Int64`: Sample rate of the original data in Hz
- `method::Symbol`: Classification method used (`:svm`, `:lda`, `:crossnobis`, etc.)
- `class_coding::Symbol`: Multi-class coding scheme (`:one_vs_one`, `:one_vs_all`, `:binary`)
- `confusion_matrix::Union{Array{Float64, 3}, Nothing}`: Confusion matrices [time × true_class × predicted_class] (if computed)
- `raw_predictions::Union{Array{Float64, 4}, Nothing}`: Raw predictions [iteration × fold × time × class] (if saved)
- `analysis_info::AnalysisInfo`: Analysis information and preprocessing metadata
"""
mutable struct DecodedData <: SingleDataFrameEeg
    file::String
    condition_names::Vector{String}
    times::Vector{Float64}
    average_score::Vector{Float64}
    stderror::Union{Vector{Float64}, Nothing}
    chance_level::Float64
    n_classes::Int64
    n_iterations::Int64
    n_folds::Int64
    channels::Vector{Symbol}
    layout::Layout
    sample_rate::Int64
    method::Symbol
    class_coding::Symbol
    confusion_matrix::Union{Array{Float64, 3}, Nothing}
    raw_predictions::Union{Array{Float64, 4}, Nothing}
    analysis_info::AnalysisInfo
end

"""
    DecodedData constructor with default values

Create a DecodedData object with required fields.
"""
function DecodedData(
    file::String,
    condition_names::Vector{String},
    times::Vector{Float64},
    average_score::Vector{Float64},
    channels::Vector{Symbol},
    layout::Layout,
    sample_rate::Int64,
    method::Symbol;
    stderror::Union{Vector{Float64}, Nothing} = nothing,
    chance_level::Float64 = 1.0 / length(condition_names),
    n_iterations::Int64 = 100,
    n_folds::Int64 = 3,
    class_coding::Symbol = :one_vs_one,
    confusion_matrix::Union{Array{Float64, 3}, Nothing} = nothing,
    raw_predictions::Union{Array{Float64, 4}, Nothing} = nothing,
    analysis_info::AnalysisInfo = AnalysisInfo(),
)
    n_classes = length(condition_names)
    return DecodedData(
        file,
        condition_names,
        times,
        average_score,
        stderror,
        chance_level,
        n_classes,
        n_iterations,
        n_folds,
        channels,
        layout,
        sample_rate,
        method,
        class_coding,
        confusion_matrix,
        raw_predictions,
        analysis_info,
    )
end

"""
    Base.show(io::IO, decoded::DecodedData)

Display DecodedData in a readable format.
"""
function Base.show(io::IO, decoded::DecodedData)
    println(io, "DecodedData:")
    println(io, "  File: $(decoded.file)")
    println(io, "  Conditions: $(join(decoded.condition_names, ", "))")
    println(io, "  Method: $(decoded.method)")
    println(io, "  Classes: $(decoded.n_classes)")
    println(io, "  Time range: $(decoded.times[1]) to $(decoded.times[end]) s")
    println(io, "  Iterations: $(decoded.n_iterations)")
    println(io, "  Cross-validation folds: $(decoded.n_folds)")
    println(io, "  Channels: $(length(decoded.channels))")
    if !isnothing(decoded.stderror)
        max_acc = maximum(decoded.average_score)
        max_idx = argmax(decoded.average_score)
        println(io, "  Max accuracy: $(round(max_acc, digits=3)) at $(round(decoded.times[max_idx], digits=3)) s")
    end
end

