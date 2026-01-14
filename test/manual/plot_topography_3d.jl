using eegfun
using DataFrames

# Create some test data
println("Creating test data...")

# Create a simple layout
layout = eegfun.read_layout("./data/layouts/biosemi/biosemi64.csv")
eegfun.polar_to_cartesian_xyz!(layout)

# Create synthetic ERP data with a clear pattern
n_samples = 100
times = range(-0.2, 0.8, length = n_samples)
n_channels = nrow(layout.data)

# Create data with a clear spatial pattern (frontal positivity)
# Build as NamedTuple to preserve order, with :time first (required for metadata detection)
data_pairs = [:time => collect(times)]

for (idx, ch) in enumerate(layout.data.label)
    # Create a simple ERP-like pattern with spatial variation
    # Frontal channels (negative z) get positive values
    z_pos = layout.data.z3[idx]
    amplitude = -z_pos * 10.0  # Frontal (negative z) = positive amplitude
    
    # Add some temporal structure
    signal = amplitude * exp.(-((times .- 0.2).^2) / (2 * 0.1^2))
    push!(data_pairs, ch => signal)
end

# Create DataFrame from NamedTuple to preserve column order (:time first)
df = DataFrame(NamedTuple(data_pairs))
erp = eegfun.ErpData("test_data", 1, "Test Condition", df, layout, 250, eegfun.AnalysisInfo(), 10)

println("Plotting 3D topography...")
# Plot at a specific time point
fig, ax = eegfun.plot_topography_3d(
    erp;
    sample_selection = eegfun.samples((0.15, 0.25)),
    sphere_resolution = 200,
    show_electrodes = true,
    display_plot = true,
)

println("Done! Check the plot window.")
