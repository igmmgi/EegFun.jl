"""
Types for Representational Similarity Analysis (RSA) results.

This module defines types for storing RSA results, including Representational
Dissimilarity Matrices (RDMs) and correlations with model RDMs.
"""

"""
    RsaData

Stores RSA analysis results for a single participant.

This type represents the results of Representational Similarity Analysis,
including Representational Dissimilarity Matrices (RDMs) computed at each
time point, correlations with model RDMs, and metadata about the analysis.

# Fields
- `file::String`: Source filename
- `condition_names::Vector{String}`: Names of conditions/bins analyzed
- `times::Vector{Float64}`: Time points in seconds where RDMs were computed
- `rdm::Array{Float64, 3}`: Representational Dissimilarity Matrices [time × condition × condition]
- `dissimilarity_measure::Symbol`: Measure used (:correlation, :euclidean, :mahalanobis, etc.)
- `channels::Vector{Symbol}`: Channel names used in the analysis
- `layout::Layout`: Layout object containing electrode positioning information
- `sample_rate::Int64`: Sample rate of the original data in Hz
- `model_correlations::Union{Array{Float64, 2}, Nothing}`: Correlations with model RDMs [time × model]
- `model_names::Union{Vector{String}, Nothing}`: Names of models compared (if any)
- `p_values::Union{Array{Float64, 2}, Nothing}`: P-values for model correlations [time × model]
- `analysis_info::AnalysisInfo`: Analysis information and preprocessing metadata
"""
mutable struct RsaData <: SingleDataFrameEeg
    file::String
    condition_names::Vector{String}
    times::Vector{Float64}
    rdm::Array{Float64, 3}  # [time × condition × condition]
    dissimilarity_measure::Symbol
    channels::Vector{Symbol}
    layout::Layout
    sample_rate::Int64
    model_correlations::Union{Array{Float64, 2}, Nothing}  # [time × model]
    model_names::Union{Vector{String}, Nothing}
    p_values::Union{Array{Float64, 2}, Nothing}  # [time × model]
    analysis_info::AnalysisInfo
end

"""
    RsaData constructor with default values

Create an RsaData object with required fields.
"""
function RsaData(
    file::String,
    condition_names::Vector{String},
    times::Vector{Float64},
    rdm::Array{Float64, 3},
    dissimilarity_measure::Symbol,
    channels::Vector{Symbol},
    layout::Layout,
    sample_rate::Int64;
    model_correlations::Union{Array{Float64, 2}, Nothing} = nothing,
    model_names::Union{Vector{String}, Nothing} = nothing,
    p_values::Union{Array{Float64, 2}, Nothing} = nothing,
    analysis_info::AnalysisInfo = AnalysisInfo(),
)
    return RsaData(
        file,
        condition_names,
        times,
        rdm,
        dissimilarity_measure,
        channels,
        layout,
        sample_rate,
        model_correlations,
        model_names,
        p_values,
        analysis_info,
    )
end

"""
    Base.show(io::IO, rsa::RsaData)

Display RsaData in a readable format.
"""
function Base.show(io::IO, rsa::RsaData)
    n_conditions = length(rsa.condition_names)
    n_times = length(rsa.times)
    
    println(io, "RsaData:")
    println(io, "  File: $(rsa.file)")
    println(io, "  Conditions: $(join(rsa.condition_names, ", ")) (n=$(n_conditions))")
    println(io, "  Dissimilarity measure: $(rsa.dissimilarity_measure)")
    println(io, "  Time points: $(n_times) ($(rsa.times[1]) to $(rsa.times[end]) s)")
    println(io, "  Channels: $(length(rsa.channels))")
    
    if !isnothing(rsa.model_correlations) && !isnothing(rsa.model_names)
        println(io, "  Model comparisons: $(length(rsa.model_names))")
        for (idx, model_name) in enumerate(rsa.model_names)
            max_corr = maximum(rsa.model_correlations[:, idx])
            max_idx = argmax(rsa.model_correlations[:, idx])
            println(io, "    - $model_name: max r=$(round(max_corr, digits=3)) at $(round(rsa.times[max_idx], digits=3)) s")
        end
    end
end

