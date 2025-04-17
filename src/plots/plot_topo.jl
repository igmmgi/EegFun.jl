##########################################
# 2D topographic plot
##########################################
function plot_topoplot!(
    fig, 
    ax,
    dat,
    layout;
    xlim = nothing,
    ylim = nothing,
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict(),
    topo_kwargs = Dict(),
    colorbar_kwargs = Dict(),
)

    if (:x2 ∉ propertynames(layout) || :y2 ∉ propertynames(layout))
        polar_to_cartesian_xy!(layout)
    end

    head_default_kwargs = Dict(:color => :black, :linewidth => 2)
    head_kwargs = merge(head_default_kwargs, head_kwargs)

    point_default_kwargs = Dict(:plot_points => true, :marker => :circle, :markersize => 12, :color => :black)
    point_kwargs = merge(point_default_kwargs, point_kwargs)

    label_default_kwargs =
        Dict(:plot_labels => true, :fontsize => 20, :color => :black, :color => :black, :xoffset => 0, :yoffset => 0)
    label_kwargs = merge(label_default_kwargs, label_kwargs)
    xoffset = pop!(label_kwargs, :xoffset)
    yoffset = pop!(label_kwargs, :yoffset)

    topo_default_kwargs = Dict(:colormap => :jet, :gridscale => 300)
    topo_kwargs = merge(topo_default_kwargs, topo_kwargs)
    gridscale = pop!(topo_kwargs, :gridscale)

    colorbar_default_kwargs = Dict(:plot_colorbar => true, :width => 30)
    colorbar_kwargs = merge(colorbar_default_kwargs, colorbar_kwargs)
    plot_colorbar = pop!(colorbar_kwargs, :plot_colorbar)

    if isnothing(xlim)
        xlim = [dat.time[1], dat.time[end]]
        xlim_idx = 1:nrow(dat)
    end

    # convert xlim to index
    xlim_idx = find_idx_range(dat.time, xlim[1], xlim[2])

    # interpolate data
    data = data_interpolation_topo(
        mean.(eachcol(dat[xlim_idx, layout.label])),
        permutedims(Matrix(layout[!, [:x2, :y2]])),
        gridscale,
    )

    if isnothing(ylim)
        ylim = minimum(data[.!isnan.(data)]), maximum(data[.!isnan.(data)])
    end

    radius = 88 # mm
    co = contourf!(
        range(-radius * 2, radius * 2, length = gridscale),
        range(-radius * 2, radius * 2, length = gridscale),
        data,
        levels = range(ylim[1], ylim[2], div(gridscale, 2));
        topo_kwargs...,
    )

    if plot_colorbar
        Colorbar(fig[1, 2], co; colorbar_kwargs...)
    end

    # head shape
    plot_layout_2d!(fig, ax, layout, head_kwargs = head_kwargs, point_kwargs = point_kwargs, label_kwargs = label_kwargs)

    display(GLMakie.Screen(), fig)
    return fig, ax

end

function plot_topoplot(
    dat,
    layout;
    kwargs...
)
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_topoplot!(fig, ax, dat, layout; kwargs...)
    display(fig)
    return fig, ax
end

"""
    plot_topoplot(dat::ContinuousData; kwargs...)

Create a topographic plot from continuous EEG data.

# Arguments
- `dat`: ContinuousData object
- `kwargs...`: Additional keyword arguments passed to plot_topoplot!

# Returns
- Figure and Axis objects
"""
function plot_topoplot(dat::ContinuousData; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_topoplot!(fig, ax, dat.data, dat.layout; kwargs...)
    display(fig)
    return fig, ax
end

"""
    plot_topoplot!(fig, ax, dat::ContinuousData; kwargs...)

Add a topographic plot to existing figure/axis from continuous EEG data.

# Arguments
- `fig`: Figure object
- `ax`: Axis object
- `dat`: ContinuousData object
- `kwargs...`: Additional keyword arguments
"""
function plot_topoplot!(fig, ax, dat::ContinuousData; kwargs...)
    plot_topoplot!(fig, ax, dat.data, dat.layout; kwargs...)
end

"""
    plot_topoplot(dat::EpochData, epoch::Int; kwargs...)

Create a topographic plot from epoched EEG data.

# Arguments
- `dat`: EpochData object
- `kwargs...`: Additional keyword arguments passed to plot_topoplot!

# Returns
- Figure and Axis objects
"""
function plot_topoplot(dat::EpochData, epoch::Int; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])
    println("heheheheheheheheheheh")
    plot_topoplot!(fig, ax, dat.data[epoch], dat.layout; kwargs...)
    display(fig)
    return fig, ax
end

"""
    plot_topoplot!(fig, ax, dat::EpochData, epoch::Int; kwargs...)

Add a topographic plot to existing figure/axis from epoched EEG data.

# Arguments
- `fig`: Figure object
- `ax`: Axis object
- `dat`: EpochData object
- `kwargs...`: Additional keyword arguments
"""
function plot_topoplot!(fig, ax, dat::EpochData, epoch::Int; kwargs...)
    plot_topoplot!(fig, ax, dat.data[epoch], dat.layout; kwargs...)
end

"""
    plot_topoplot(dat::ErpData; kwargs...)

Create a topographic plot from ERP data.

# Arguments
- `dat`: ErpData object
- `kwargs...`: Additional keyword arguments passed to plot_topoplot!

# Returns
- Figure and Axis objects
"""
function plot_topoplot(dat::ErpData; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_topoplot!(fig, ax, dat.data, dat.layout; kwargs...)
    display(fig)
    return fig, ax
end

"""
    plot_topoplot!(fig, ax, dat::ErpData; kwargs...)

Add a topographic plot to existing figure/axis from ERP data.

# Arguments
- `fig`: Figure object
- `ax`: Axis object
- `dat`: ErpData object
- `kwargs...`: Additional keyword arguments
"""
function plot_topoplot!(fig, ax, dat::ErpData; kwargs...)
    plot_topoplot!(fig, ax, dat.data, dat.layout; kwargs...)
end


# Basic Tests
# TODO: Implement proper tests
# layout = read_layout("./layouts/biosemi72.csv");
# subject = 3
# dat = read_bdf("../Flank_C_$(subject).bdf");
# dat = create_eeg_dataframe(dat, layout);
# epochs = []
# for (idx, epoch) in enumerate([1, 4, 5, 3])
#      push!(epochs, extract_epochs(dat, idx, epoch, -2, 4))
# end
# 
# # Continuous Data
# plot_topoplot(dat)
# plot_topoplot(dat.data, dat.layout)
# 
# # Epoch Data
# epoch = extract_epochs(dat, 1, 1, -2, 4)
# plot_topoplot(epoch, 1) # 1st epoch
# plot_topoplot(epoch, 2) # 2nd epoch
# 
# # ERP Data
# erp = average_epochs(epochs)
# plot_topoplot(erp)
# 
# # Try some keyword arguments
# plot_topoplot(erp, xlim = (-0.1, 0.2), ylim = (-10, 10), 
#     head_kwargs = Dict(:linewidth => 5),
#     point_kwargs = Dict(:markersize => 20),
#     label_kwargs = Dict(:fontsize => 12),
#     topo_kwargs = Dict(:colormap => :viridis),
#     colorbar_kwargs = Dict(:width => 20),
# )


"""
    circle_mask!(dat::Matrix{<:AbstractFloat}, grid_scale::Int)

Apply a circular mask to a topographic data matrix, setting values outside
the circle to NaN.

# Arguments
- `dat::Matrix{Float64}`: Data matrix to mask (modified in-place)
- `grid_scale::Int`: Size of the grid (must be square matrix)

# Throws
- `ArgumentError`: If grid_scale is not positive
- `DimensionMismatch`: If matrix is not square
"""
function circle_mask!(dat::Matrix{<:AbstractFloat}, grid_scale::Int)
    if grid_scale <= 0
        throw(ArgumentError("grid_scale must be positive"))
    end
    if size(dat) != (grid_scale, grid_scale)
        throw(DimensionMismatch("Data matrix must be $(grid_scale)×$(grid_scale)"))
    end

    center = grid_scale / 2
    @inbounds for col = 1:grid_scale
        for row = 1:grid_scale
            x_dist = center - col
            y_dist = center - row
            if sqrt(x_dist^2 + y_dist^2) > center
                dat[col, row] = NaN
            end
        end
    end
    return dat
end



"""
    data_interpolation_topo(dat::Vector{Float64}, points::Matrix{Float64}, grid_scale::Int)

Interpolate EEG data using scattered interpolation 

# Arguments
- `dat::Vector{<:AbstractFloat}`: EEG values at electrode positions
- `points::Matrix{<:AbstractFloat}`: 2×N matrix of electrode coordinates
- `grid_scale::Int`: Size of the output grid
"""
function data_interpolation_topo(dat::Vector{<:AbstractFloat}, points::Matrix{<:AbstractFloat}, grid_scale::Int)
    # Check input data
    if any(isnan, dat) || any(isinf, dat)
        throw(ArgumentError("Input data contains NaN or Inf values"))
    end
    if any(isnan, points) || any(isinf, points)
        throw(ArgumentError("Input points contain NaN or Inf values"))
    end

    # Create grid based on the range of input points
    x_range = range(minimum(points[1, :]), maximum(points[1, :]), length = grid_scale)
    y_range = range(minimum(points[2, :]), maximum(points[2, :]), length = grid_scale)

    # Create regular grid
    grid_points = zeros(2, grid_scale^2)
    @inbounds for (idx, (i, j)) in enumerate(Iterators.product(x_range, y_range))
        grid_points[1, idx] = i
        grid_points[2, idx] = j
    end

    # Perform interpolation
    try
        itp = ScatteredInterpolation.interpolate(Multiquadratic(), points, dat)
        result = reshape(ScatteredInterpolation.evaluate(itp, grid_points), grid_scale, grid_scale)
        circle_mask!(result, grid_scale)
        return result
    catch e
        throw(ErrorException("Interpolation failed: $(e)"))
    end

end