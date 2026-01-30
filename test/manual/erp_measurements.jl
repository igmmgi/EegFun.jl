"""
Tutorial: ERP Measurement Options

This script provides an introduction to the ERP measurement capabilities 
in EegFun for extracting quantitative features from ERP data.

1. Amplitude measurements (mean, peak)
2. Latency measurements (peak, fractional)
3. Area/integral measurements
4. Peak-to-peak measurements
"""

using EegFun

input_dir = "./data/files/erps"
file_pattern = "erps_good"

# ----------------------------------------------------------------------------
# Amplitude Measurements
# ----------------------------------------------------------------------------

dat = EegFun.load_data("./data/files/erps/example1_erps_good.jld2")

EegFun.plot_erp(
    dat,
    condition_selection = EegFun.conditions([1]),
    channel_selection = EegFun.channels([:Pz]),
    # baseline_interval = (0, 0),  
    baseline_interval = -0.2:0,
)


# Mean amplitude in a time window
println("\n=== Mean Amplitude ===")
mean_amp = EegFun.erp_measurements(
    file_pattern,
    "max_peak_amplitude",
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1]),
    channel_selection = EegFun.channels([:Pz, :Cz, :Fz]),
    analysis_interval = (0.6, 0.8),  # P3 window: 300-500ms
    baseline_interval = (-0.2, 0.0),
)
println(mean_amp)

# Maximum peak amplitude (e.g., for P3)
println("\n=== Maximum Peak Amplitude ===")
max_peak_amp = EegFun.erp_measurements(
    file_pattern,
    "max_peak_amplitude",
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1]),
    channel_selection = EegFun.channels([:Pz, :Cz]),
    analysis_interval = EegFun.samples((0.25, 0.6)),  # P3 search window
    baseline_interval = EegFun.samples((-0.2, 0.0)),
    local_window = 3,  # Robust peak detection
)
println(max_peak_amp)

# Minimum peak amplitude (e.g., for N1 or N2)
println("\n=== Minimum Peak Amplitude ===")
min_peak_amp = EegFun.erp_measurements(
    file_pattern,
    "min_peak_amplitude",
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1]),
    channel_selection = EegFun.channels([:Fz, :Cz]),
    analysis_interval = EegFun.samples((0.15, 0.25)),  # N2 search window
    baseline_interval = EegFun.samples((-0.2, 0.0)),
    local_window = 3,
)
println(min_peak_amp)

# ----------------------------------------------------------------------------
# Latency Measurements
# ----------------------------------------------------------------------------

# Maximum peak latency
println("\n=== Maximum Peak Latency ===")
max_peak_lat = EegFun.erp_measurements(
    file_pattern,
    "max_peak_latency",
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1]),
    channel_selection = EegFun.channels([:Pz]),
    analysis_interval = EegFun.samples((0.25, 0.6)),
    baseline_interval = EegFun.samples((-0.2, 0.0)),
    local_window = 3,
)
println(max_peak_lat)

# Minimum peak latency
println("\n=== Minimum Peak Latency ===")
min_peak_lat = EegFun.erp_measurements(
    file_pattern,
    "min_peak_latency",
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1]),
    channel_selection = EegFun.channels([:Fz]),
    analysis_interval = EegFun.samples((0.15, 0.25)),
    baseline_interval = EegFun.samples((-0.2, 0.0)),
    local_window = 3,
)
println(min_peak_lat)

# Onset latency (50% fractional peak)
println("\n=== Onset Latency (50% fractional peak) ===")
onset_lat = EegFun.erp_measurements(
    file_pattern,
    "fractional_peak_latency",
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1]),
    channel_selection = EegFun.channels([:Pz]),
    analysis_interval = EegFun.samples((0.25, 0.6)),
    baseline_interval = EegFun.samples((-0.2, 0.0)),
    local_window = 3,
    fractional_peak_fraction = 0.5,  # 50% of peak amplitude
    fractional_peak_direction = :onset,
)
println(onset_lat)

# Offset latency (50% fractional peak)
println("\n=== Offset Latency (50% fractional peak) ===")
offset_lat = EegFun.erp_measurements(
    file_pattern,
    "fractional_peak_latency",
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1]),
    channel_selection = EegFun.channels([:Pz]),
    analysis_interval = EegFun.samples((0.25, 0.6)),
    baseline_interval = EegFun.samples((-0.2, 0.0)),
    local_window = 3,
    fractional_peak_fraction = 0.5,
    fractional_peak_direction = :offset,
)
println(offset_lat)

# Fractional area latency (50% of area)
println("\n=== Fractional Area Latency (50%) ===")
fract_area_lat = EegFun.erp_measurements(
    file_pattern,
    "fractional_area_latency",
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1]),
    channel_selection = EegFun.channels([:Pz]),
    analysis_interval = EegFun.samples((0.3, 0.5)),
    baseline_interval = EegFun.samples((-0.2, 0.0)),
    fractional_area_fraction = 0.5,  # 50% of total area
)
println(fract_area_lat)

# ----------------------------------------------------------------------------
# Area/Integral Measurements
# ----------------------------------------------------------------------------

# Rectified area (absolute value integration)
println("\n=== Rectified Area ===")
rect_area = EegFun.erp_measurements(
    file_pattern,
    "rectified_area",
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1]),
    channel_selection = EegFun.channels([:Pz, :Cz]),
    analysis_interval = EegFun.samples((0.3, 0.5)),
    baseline_interval = EegFun.samples((-0.2, 0.0)),
)
println(rect_area)

# Integral (signed area under curve)
println("\n=== Integral (Signed Area) ===")
integral_area = EegFun.erp_measurements(
    file_pattern,
    "integral",
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1]),
    channel_selection = EegFun.channels([:Pz, :Cz]),
    analysis_interval = EegFun.samples((0.3, 0.5)),
    baseline_interval = EegFun.samples((-0.2, 0.0)),
)
println(integral_area)

# Positive area only
println("\n=== Positive Area ===")
pos_area = EegFun.erp_measurements(
    file_pattern,
    "positive_area",
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1]),
    channel_selection = EegFun.channels([:Pz]),
    analysis_interval = EegFun.samples((0.3, 0.5)),
    baseline_interval = EegFun.samples((-0.2, 0.0)),
)
println(pos_area)

# Negative area only
println("\n=== Negative Area ===")
neg_area = EegFun.erp_measurements(
    file_pattern,
    "negative_area",
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1]),
    channel_selection = EegFun.channels([:Fz]),
    analysis_interval = EegFun.samples((0.15, 0.25)),
    baseline_interval = EegFun.samples((-0.2, 0.0)),
)
println(neg_area)

# ----------------------------------------------------------------------------
# Peak-to-Peak Measurements
# ----------------------------------------------------------------------------

# Peak-to-peak amplitude
println("\n=== Peak-to-Peak Amplitude ===")
p2p_amp = EegFun.erp_measurements(
    file_pattern,
    "peak_to_peak_amplitude",
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1]),
    channel_selection = EegFun.channels([:Pz, :Cz]),
    analysis_interval = EegFun.samples((0.0, 0.6)),
    baseline_interval = EegFun.samples((-0.2, 0.0)),
    local_window = 3,
)
println(p2p_amp)

# Peak-to-peak latency
println("\n=== Peak-to-Peak Latency ===")
p2p_lat = EegFun.erp_measurements(
    file_pattern,
    "peak_to_peak_latency",
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1]),
    channel_selection = EegFun.channels([:Pz, :Cz]),
    analysis_interval = EegFun.samples((0.0, 0.6)),
    baseline_interval = EegFun.samples((-0.2, 0.0)),
    local_window = 3,
)
println(p2p_lat)

println("\n=== All measurements complete! ===")
