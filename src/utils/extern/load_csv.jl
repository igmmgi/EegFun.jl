using eegfun
using CSV
using DataFrames
using Statistics

"""
    load_csv(data_dir::String; 
                      file::String="fieldtrip_data",
                      condition::Int=1,
                      condition_name::String="fieldtrip",
                      fsample::Union{Int,Nothing}=nothing) -> EpochData

Load FieldTrip epochs data saved as CSV and convert to EpochData.

The CSV file should be created using the MATLAB function `matlab_to_csv.m`.

# Arguments
- `data_dir::String`: Directory containing `epochs_data.csv`
  - The CSV has columns: [trial, time, channel1, channel2, ...]

# Keyword Arguments
- `file::String="fieldtrip_data"`: Source filename for EpochData
- `condition::Int=1`: Condition number
- `condition_name::String="fieldtrip"`: Condition name
- `fsample::Union{Int,Nothing}=nothing`: Sample rate in Hz (if not provided, will be estimated from time column)

# Returns
- `EpochData`: EpochData object with the loaded data

# Example
```julia
epochs = load_csv("/path/to/csv/files", fsample=256)
```
"""
function load_csv(
    data_dir::String;
    file::String = "csv_data",
    condition::Int = 1,
    condition_name::String = "csv_data",
    fsample::Union{Int,Nothing} = nothing,
)

    data_file = joinpath(data_dir, "epochs_data.csv")
    if !isfile(data_file)
        error("Data file not found: $data_file")
    end
    all_data = CSV.read(data_file, DataFrame)

    # Column names: trial, time, channel1, channel2, ...
    actual_cols = names(all_data)
    
    # First two columns should be 'trial' and 'time'
    # CSV.read returns column names as strings, so compare with strings
    println("hahahah")
    if length(actual_cols) < 2 || actual_cols[1] != "trial" || actual_cols[2] != "time"
        error("Expected first two columns to be 'trial' and 'time', got: $(actual_cols[1:min(2, length(actual_cols))])")
    end
    
    channel_labels = [Symbol(col) for col in actual_cols[3:end]]
    n_trials = length(unique(all_data.trial))
    
    # Estimate sample rate from time column if not provided
    if isnothing(fsample)
        # Get time differences for first trial
        trial1_data = Base.filter(row -> row.trial == 1, all_data)
        if length(trial1_data.time) > 1
            time_diffs = diff(sort(trial1_data.time))
            # Most common time difference should be 1/fsample
            median_dt = median(time_diffs)
            fsample = Int(round(1.0 / median_dt))
            @info "Estimated sample rate from time column: $fsample Hz"
        else
            error("Cannot estimate sample rate. Please provide fsample parameter.")
        end
    end

    # Create layout
    layout_df = DataFrame(:label => channel_labels, :inc => zeros(length(channel_labels)), :azi => zeros(length(channel_labels)))
    layout = eegfun.Layout(layout_df, nothing, nothing)

    # Create AnalysisInfo
    analysis_info = eegfun.AnalysisInfo()

    # Split data into trials
    trial_dfs = Vector{DataFrame}(undef, n_trials)
    for trial_idx = 1:n_trials
        # Filter data for this trial
        trial_data = Base.filter(row -> row.trial == trial_idx, all_data)

        # Create DataFrame with time and channel columns
        trial_df = DataFrame(:time => trial_data.time, :epoch => fill(trial_idx, nrow(trial_data)))

        # Add channel columns
        for ch in channel_labels
            trial_df[!, ch] = trial_data[!, ch]
        end

        trial_dfs[trial_idx] = trial_df
    end

    # Create EpochData
    return eegfun.EpochData(file, condition, condition_name, trial_dfs, layout, fsample, analysis_info)
end

