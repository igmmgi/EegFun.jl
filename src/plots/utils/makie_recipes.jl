import Makie

_tf_1d_recipe_error(dat) =
    @minimal_error("Direct 1D Makie recipes (e.g. `lines`) are not supported for $(typeof(dat)); use `EegFun.plot_tf` instead.")

Makie.convert_arguments(::Makie.PointBased, dat::Union{TimeFreqData,TimeFreqEpochData}) = _tf_1d_recipe_error(dat)
Makie.convert_arguments(::Makie.PointBased, dat::Union{TimeFreqData,TimeFreqEpochData}, ::Symbol) = _tf_1d_recipe_error(dat)
Makie.convert_arguments(::Type{<:Makie.Series}, dat::Union{TimeFreqData,TimeFreqEpochData}) = _tf_1d_recipe_error(dat)
Makie.convert_arguments(::Type{<:Makie.Series}, dat::TimeFreqEpochData, ::Symbol) = _tf_1d_recipe_error(dat)

@inline function _copy_col_to_row!(mat::Matrix{Float64}, row::Int, col::AbstractVector)
    @inbounds for j in eachindex(col)
        mat[row, j] = Float64(col[j])
    end
    return nothing
end

function Makie.convert_arguments(P::Makie.PointBased, dat::SpectrumData, channel::Symbol)
    hasproperty(dat.data, channel) || @minimal_error("Channel $channel not found in SpectrumData.")
    return Makie.convert_arguments(P, dat.data[!, :freq], dat.data[!, channel])
end

function Makie.convert_arguments(P::Makie.PointBased, dat::SingleDataFrameEeg, channel::Symbol)
    hasproperty(dat.data, channel) || @minimal_error("Channel $channel not found in data.")
    return Makie.convert_arguments(P, dat.data[!, :time], dat.data[!, channel])
end

function Makie.convert_arguments(P::Type{<:Makie.Series}, dat::MultiDataFrameEeg, channel::Symbol)
    (isempty(dat.data) || !(channel in channel_labels(dat))) && @minimal_error("Channel $channel not found in data.")
    times_vec = dat.data[1][!, :time]
    n_ep = length(dat.data)
    n_samp = length(times_vec)
    y_matrix = Matrix{Float64}(undef, n_ep, n_samp)
    for ep = 1:n_ep
        _copy_col_to_row!(y_matrix, ep, dat.data[ep][!, channel])
    end
    return Makie.convert_arguments(P, times_vec, y_matrix)
end

function Makie.convert_arguments(P::Type{<:Makie.Series}, dat::SpectrumData)
    freqs_vec = dat.data[!, :freq]
    chans = channel_labels(dat)
    y_matrix = Matrix{Float64}(undef, length(chans), length(freqs_vec))
    for (i, ch) in enumerate(chans)
        _copy_col_to_row!(y_matrix, i, dat.data[!, ch])
    end
    return Makie.convert_arguments(P, freqs_vec, y_matrix)
end

function Makie.convert_arguments(P::Type{<:Makie.Series}, dat::SingleDataFrameEeg)
    return Makie.convert_arguments(P, dat.data[!, :time], transpose(Matrix(dat)))
end

function Makie.convert_arguments(P::Makie.PointBased, layout::Layout)
    if !hasproperty(layout.data, :x2) || !hasproperty(layout.data, :y2)
        has_valid_coordinates(layout) || @minimal_error("Layout has no valid polar or Cartesian coordinates.")
        layout = copy(layout)
        _ensure_coordinates_2d!(layout)
    end
    return Makie.convert_arguments(P, layout.data[!, :x2], layout.data[!, :y2])
end
