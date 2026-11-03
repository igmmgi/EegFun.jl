import Makie

# 0. Guard 2D Time-Frequency types against 1D line/series recipes
function Makie.convert_arguments(::Makie.PointBased, dat::Union{TimeFreqData,TimeFreqEpochData})
    @minimal_error("Direct 1D Makie recipes (e.g. `lines`) are not supported for $(typeof(dat)); use `EegFun.plot_tf` instead.")
end
function Makie.convert_arguments(::Makie.PointBased, dat::Union{TimeFreqData,TimeFreqEpochData}, ::Symbol)
    @minimal_error("Direct 1D Makie recipes (e.g. `lines`) are not supported for $(typeof(dat)); use `EegFun.plot_tf` instead.")
end
function Makie.convert_arguments(::Type{<:Makie.Series}, dat::Union{TimeFreqData,TimeFreqEpochData})
    @minimal_error("Direct 1D Makie recipes (e.g. `series`) are not supported for $(typeof(dat)); use `EegFun.plot_tf` instead.")
end
function Makie.convert_arguments(::Type{<:Makie.Series}, dat::TimeFreqEpochData, ::Symbol)
    @minimal_error("Direct 1D Makie recipes (e.g. `series`) are not supported for $(typeof(dat)); use `EegFun.plot_tf` instead.")
end

# 1a. SpectrumData Single Channel (uses `:freq` instead of `:time`)
function Makie.convert_arguments(P::Makie.PointBased, dat::SpectrumData, channel::Symbol)
    if !hasproperty(dat.data, channel)
        @minimal_error("Channel $channel not found in SpectrumData.")
    end
    freqs_vec = dat.data.freq
    y_vec = dat.data[!, channel]
    return Makie.convert_arguments(P, freqs_vec, y_vec)
end

# 1b. Single Channel (Continuous / ERP)
# Allows: lines(dat, :Fz)
function Makie.convert_arguments(P::Makie.PointBased, dat::SingleDataFrameEeg, channel::Symbol)
    if !hasproperty(dat.data, channel)
        @minimal_error("Channel $channel not found in data.")
    end
    times_vec = dat.data.time
    y_vec = dat.data[!, channel]

    # Delegate back to Makie's standard conversion for 1D arrays
    return Makie.convert_arguments(P, times_vec, y_vec)
end

# 2. Multi Channel (EpochData)
# Allows: series(dat, :Fz) -> plots all epochs as separate overlapping lines
function Makie.convert_arguments(P::Type{<:Makie.Series}, dat::MultiDataFrameEeg, channel::Symbol)
    # Find channel index safely
    chans = EegFun.channel_labels(dat)
    c_idx = findfirst(==(channel), chans)
    if isnothing(c_idx)
        @minimal_error("Channel $channel not found in data.")
    end

    times_vec = dat.data[1].time

    # Array(dat) returns (epochs × samples × channels)
    # Extract the specific channel: (epochs × samples)
    # Makie's series expects (N_series × N_points)
    y_matrix = Array(dat)[:, :, c_idx]

    return Makie.convert_arguments(P, times_vec, y_matrix)
end

# 3a. SpectrumData Butterfly Plot (all channels overlaid vs frequency)
function Makie.convert_arguments(P::Type{<:Makie.Series}, dat::SpectrumData)
    freqs_vec = dat.data.freq
    chans = channel_labels(dat)
    y_matrix = Matrix{Float64}(undef, length(chans), length(freqs_vec))
    for (i, ch) in enumerate(chans)
        y_matrix[i, :] = dat.data[!, ch]
    end
    return Makie.convert_arguments(P, freqs_vec, y_matrix)
end

# 3b. Butterfly Plot (No channel specified, all channels overlaid)
# Works natively for time-domain SingleDataFrameEeg types (Continuous, ERP)
function Makie.convert_arguments(P::Type{<:Makie.Series}, dat::SingleDataFrameEeg)
    times_vec = dat.data.time
    # Matrix(dat) returns a (samples × channels) matrix
    # Makie's series expects (N_series × N_points), so we transpose to (channels × samples)
    y_matrix = transpose(Matrix(dat))
    return Makie.convert_arguments(P, times_vec, y_matrix)
end

# 4. Layout spatial plotting
# Allows: scatter(dat.layout) or lines(dat.layout)
function Makie.convert_arguments(P::Makie.PointBased, layout::EegFun.Layout)
    if !hasproperty(layout.data, :x2) || !hasproperty(layout.data, :y2)
        if !has_valid_coordinates(layout)
            @minimal_error("Layout has no valid polar or Cartesian coordinates.")
        end
        layout = copy(layout)
        _ensure_coordinates_2d!(layout)
    end
    x_vec = layout.data.x2
    y_vec = layout.data.y2
    return Makie.convert_arguments(P, x_vec, y_vec)
end
