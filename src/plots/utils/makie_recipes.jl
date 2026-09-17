import Makie

# 1. Single Channel (Continuous / ERP)
# Allows: lines(dat, :Fz)
function Makie.convert_arguments(P::Makie.PointBased, dat::SingleDataFrameEeg, channel::Symbol)
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
        error("Channel $channel not found in data.")
    end
    
    times_vec = dat.data[1].time
    
    # Array(dat) returns (epochs × samples × channels)
    # Extract the specific channel: (epochs × samples)
    # Makie's series expects (N_series × N_points)
    y_matrix = Array(dat)[:, :, c_idx]
    
    return Makie.convert_arguments(P, times_vec, y_matrix)
end

# 3. Butterfly Plot (No channel specified, all channels overlaid)
# Works natively for all SingleDataFrameEeg types (Continuous, ERP, Spectrum)
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
        error("Layout is missing x2/y2 Cartesian coordinates. Run polar_to_cartesian_xy!(layout) first.")
    end
    x_vec = layout.data.x2
    y_vec = layout.data.y2
    return Makie.convert_arguments(P, x_vec, y_vec)
end
