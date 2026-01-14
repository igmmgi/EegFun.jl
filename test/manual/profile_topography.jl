using eegfun
using BenchmarkTools
using Profile
using ProfileView
using DataFrames
using Statistics

# Get some basic data
data_file = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "recoded", "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
eegfun.polar_to_cartesian_xy!(layout_file)
dat = eegfun.read_bdf(data_file);
dat = eegfun.create_eeg_dataframe(dat, layout_file);
eegfun.rereference!(dat, :avg)
eegfun.filter_data!(dat, "hp", 1)

# Prepare data for direct function call
dat_subset = eegfun.subset(dat, sample_selection = x -> x.time .>= 7.984 .&& x.time .<= 8.168)
channel_data = vec(mean(Matrix(dat_subset.data[:, Not(:time)]), dims = 1))
layout = dat_subset.layout
gridscale = 67

println("=== Profiling spherical spline interpolation ===")
println("\nBaseline benchmark:")
@btime eegfun._data_interpolation_topo_spherical_spline($channel_data, $layout, $gridscale)

println("\n\nDetailed profiling (run for 10 seconds):")
# Warm up
for i = 1:5
    eegfun._data_interpolation_topo_spherical_spline(channel_data, layout, gridscale)
end

# Profile
Profile.clear()
@profile for i = 1:100
    eegfun._data_interpolation_topo_spherical_spline(channel_data, layout, gridscale)
end

# Show results
Profile.print(format = :flat, sortedby = :count, mincount = 10)

println("\n\n=== Breaking down individual components ===")

# Test 1: Coordinate extraction and normalization
println("\n1. Coordinate extraction and normalization:")
@btime begin
    eegfun._ensure_coordinates_3d!($layout)
    n_channels = length($channel_data)
    x3_col = $layout.data.x3::Vector{Float64}
    y3_col = $layout.data.y3::Vector{Float64}
    z3_col = $layout.data.z3::Vector{Float64}

    coords = Matrix{Float64}(undef, n_channels, 3)
    @inbounds for i = 1:n_channels
        coords[i, 1] = x3_col[i]
        coords[i, 2] = y3_col[i]
        coords[i, 3] = z3_col[i]
    end
    coords
end

# Test 2: G matrix calculation
println("\n2. G matrix calculation:")
n_channels = length(channel_data)
coords = Matrix{Float64}(undef, n_channels, 3)
@inbounds for i = 1:n_channels
    coords[i, 1] = layout.data.x3[i]
    coords[i, 2] = layout.data.y3[i]
    coords[i, 3] = layout.data.z3[i]
end
coords_unit = Matrix{Float64}(undef, n_channels, 3)
@inbounds for i = 1:n_channels
    norm_factor = sqrt(coords[i, 1]^2 + coords[i, 2]^2 + coords[i, 3]^2)
    coords_unit[i, 1] = coords[i, 1] / norm_factor
    coords_unit[i, 2] = coords[i, 2] / norm_factor
    coords_unit[i, 3] = coords[i, 3] / norm_factor
end

@btime begin
    cosang = $coords_unit * $coords_unit'
    G = eegfun._calc_g_matrix(cosang)
end

# Test 3: System solve
println("\n3. System solve:")
cosang = coords_unit * coords_unit'
G = eegfun._calc_g_matrix(cosang)
@inbounds for i = 1:n_channels
    G[i, i] += 1e-5
end
G_extended = Matrix{Float64}(undef, n_channels + 1, n_channels + 1)
@inbounds for i = 1:n_channels
    for j = 1:n_channels
        G_extended[i, j] = G[i, j]
    end
    G_extended[i, n_channels+1] = 1.0
    G_extended[n_channels+1, i] = 1.0
end
G_extended[n_channels+1, n_channels+1] = 0.0
data_vector = Vector{Float64}(undef, n_channels + 1)
@inbounds for i = 1:n_channels
    data_vector[i] = channel_data[i]
end
data_vector[n_channels+1] = 0.0

@btime $G_extended \ $data_vector

# Test 4: Grid point generation and projection
println("\n4. Grid point generation and stereographic projection:")
electrode_radius = mean(sqrt.(sum(coords .^ 2, dims = 2)))
x_range = range(-1.0, 1.0, length = gridscale)
y_range = range(-1.0, 1.0, length = gridscale)

@btime begin
    grid_x = vec(repeat($x_range', length($y_range), 1))
    grid_y = repeat($y_range, length($x_range))
    r_2d = sqrt.(grid_x .^ 2 .+ grid_y .^ 2)
    valid_mask = r_2d .<= 1.0
    valid_indices = findall(valid_mask)

    valid_x = grid_x[valid_indices]
    valid_y = grid_y[valid_indices]
    valid_r = r_2d[valid_indices]

    r_norm = valid_r ./ 1.0
    z3 = $electrode_radius .* (1.0 .- r_norm .^ 2) ./ (1.0 .+ r_norm .^ 2)
    x3 = valid_y .* (1.0 .+ z3 ./ $electrode_radius)
    y3 = valid_x .* (1.0 .+ z3 ./ $electrode_radius)

    grid_points_3d = hcat(x3, y3, z3)

    n_valid = length(valid_indices)
    grid_points_unit = Matrix{Float64}(undef, n_valid, 3)
    @inbounds for i = 1:n_valid
        norm_val = sqrt(grid_points_3d[i, 1]^2 + grid_points_3d[i, 2]^2 + grid_points_3d[i, 3]^2)
        if norm_val > 0
            grid_points_unit[i, 1] = grid_points_3d[i, 1] / norm_val
            grid_points_unit[i, 2] = grid_points_3d[i, 2] / norm_val
            grid_points_unit[i, 3] = grid_points_3d[i, 3] / norm_val
        end
    end
    grid_points_unit
end

# Test 5: Interpolation calculation
println("\n5. Final interpolation calculation:")
weights = G_extended \ data_vector
grid_x = vec(repeat(x_range', length(y_range), 1))
grid_y = repeat(y_range, length(x_range))
r_2d = sqrt.(grid_x .^ 2 .+ grid_y .^ 2)
valid_mask = r_2d .<= 1.0
valid_indices = findall(valid_mask)
valid_x = grid_x[valid_indices]
valid_y = grid_y[valid_indices]
valid_r = r_2d[valid_indices]
r_norm = valid_r ./ 1.0
z3 = electrode_radius .* (1.0 .- r_norm .^ 2) ./ (1.0 .+ r_norm .^ 2)
x3 = valid_y .* (1.0 .+ z3 ./ electrode_radius)
y3 = valid_x .* (1.0 .+ z3 ./ electrode_radius)
grid_points_3d = hcat(x3, y3, z3)
n_valid = length(valid_indices)
grid_points_unit = Matrix{Float64}(undef, n_valid, 3)
@inbounds for i = 1:n_valid
    norm_val = sqrt(grid_points_3d[i, 1]^2 + grid_points_3d[i, 2]^2 + grid_points_3d[i, 3]^2)
    grid_points_unit[i, 1] = grid_points_3d[i, 1] / norm_val
    grid_points_unit[i, 2] = grid_points_3d[i, 2] / norm_val
    grid_points_unit[i, 3] = grid_points_3d[i, 3] / norm_val
end

@btime begin
    # New optimized interpolation path
    cosang_grid = $grid_points_unit * $coords_unit'
    factors = eegfun._get_g_factors(15)
    @inbounds for i in eachindex(cosang_grid)
        cosang_grid[i] = eegfun._legendre_val(clamp(cosang_grid[i], -1.0, 1.0), factors)
    end
    w_channels = @view weights[1:n_channels]
    interpolated_valid = cosang_grid * w_channels
    interpolated_valid .+= weights[n_channels+1]
end

println("\n\n=== Summary ===")
println("Total function time: ~435ms")
println("Check which components dominate the runtime above.")
