"""
FieldTrip data importer for EegFun.jl

Loads FieldTrip .mat files and converts to EegFun data structures.
"""


"""
    load_fieldtrip(filepath::String, layout::Layout) → ContinuousData/EpochData/ErpData

Load FieldTrip .mat file and convert to EegFun data structure.

# Arguments
- `filepath::String`: Path to FieldTrip .mat file
- `layout::Layout`: Channel layout (FieldTrip doesn't store layout in data files)

# Returns
- `ContinuousData`: For single-trial continuous data
- `EpochData`: For multi-trial epoched data
- `ErpData`: For averaged ERP data

# Example
```julia
using EegFun
# We need to specify a layout to pass to our EegFun data types
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi64.csv")
data = EegFun.load_fieldtrip("./resources/data/fieldtrip/continuous.mat", layout)
```
"""
function load_fieldtrip(filepath::String, layout::Layout)
    @info "Loading FieldTrip .mat file: $filepath"

    mat_data = matread(filepath)

    # We need to try and figure our what kind of data is in the Matlab mat file
    ft = nothing
    is_erp = false
    for (key, value) in mat_data
        if value isa Dict
            # Check for avg (ERP) first, then trial (continuous/epoched)
            if haskey(value, "avg") && haskey(value, "time")
                ft = value
                is_erp = true
                @info "Found FieldTrip ERP data in variable: $key"
                break
            elseif haskey(value, "trial") && haskey(value, "fsample")
                ft = value
                @info "Found FieldTrip data in variable: $key"
                break
            end
        end
    end

    if ft === nothing
        @minimal_error("Could not find FieldTrip data structure in file. Expected a dict with 'trial'+'fsample' or 'avg'+'time' fields.")
    end

    # Handle ERP data differently (has 'avg' instead of 'trial', infer sample rate from time)
    if is_erp
        # Infer sample rate from time vector
        time_vec = vec(ft["time"])
        if length(time_vec) > 1
            dt = time_vec[2] - time_vec[1]
            sample_rate = Int(round(1.0 / dt))
        else
            error("Invalid ERP data: time vector has only $(length(time_vec)) point(s). Cannot infer sample rate.")
        end
        @info "Detected averaged ERP data, inferred sample rate: $sample_rate Hz"
        return _fieldtrip_to_erpdata(ft, filepath, layout, sample_rate)
    end

    # Handle trial-based data (continuous (1 trial) or epoched)
    n_trials = length(ft["trial"])
    sample_rate = Int(ft["fsample"])
    @info "Found $n_trials trial(s), sample rate: $sample_rate Hz"
    # Detect data type and convert
    if n_trials == 1 # Single trial: continuous data
        @info "Detected continuous data"
        return _fieldtrip_to_continuousdata(ft, filepath, layout, sample_rate)
    else # Multiple trials: epoched data
        @info "Detected epoched data with $n_trials trials"
        return _fieldtrip_to_epochdata(ft, filepath, layout, sample_rate)
    end
end


# === Helper Functions ===
"""
Remove BDF Status channel if present (last channel named "Status").
Returns updated channel labels and data matrix.
"""
function _remove_status_channel(ch_labels::Vector, data_matrix::AbstractMatrix)
    if !isempty(ch_labels) && lowercase(string(ch_labels[end])) == "status"
        @info "Detected BDF Status channel, removing it"
        return ch_labels[1:end-1], data_matrix[1:end-1, :]
    end
    return ch_labels, data_matrix
end

"""
Extract time vector from FieldTrip data.
Handles both cell array format (one vector per trial) and shared numeric matrix format.
"""
function _extract_time_vector(time_data)
    if time_data isa AbstractMatrix && eltype(time_data) <: Number && size(time_data, 1) == 1
        # Shared numeric time vector: (1 x n_timepoints)
        return vec(time_data)
    else
        # Cell array: extract first trial's time vector
        # This is needed for our custom stripped down version
        return vec(time_data[1])
    end
end

"""
Create metadata (condition name and analysis info) from filepath.
"""
function _create_metadata(filepath::String)
    condition_name = basename(filepath)
    analysis_info = AnalysisInfo()
    return condition_name, analysis_info
end

"""
Create trigger columns (zeros for ERP/continuous without events).
"""
function _create_trigger_columns(n_timepoints::Int)
    trigger_vec = zeros(Int, n_timepoints)
    trigger_info_vec = fill("", n_timepoints)
    return trigger_vec, trigger_info_vec
end


"""
Convert FieldTrip continuous data to EegFun ContinuousData.
"""
function _fieldtrip_to_continuousdata(ft::Dict, filepath::String, layout::Layout, sample_rate::Int)
    # Extract channel labels and data
    ch_labels = vec(ft["label"])
    data_matrix = ft["trial"][1]  # [channels × timepoints]
    time_vec = vec(ft["time"][1])

    # Remove Status channel if present
    ch_labels, data_matrix = _remove_status_channel(ch_labels, data_matrix)

    ch_names = Symbol.(ch_labels)
    n_channels, n_timepoints = size(data_matrix)

    @info "Data: $n_channels channels × $n_timepoints timepoints"

    # Create DataFrame with time column
    df = DataFrame(data_matrix', ch_names)
    insertcols!(df, 1, :time => time_vec)

    # Extract events from cfg.event if available (BDF files use struct-of-arrays format)
    trigger_vec, trigger_info_vec = _create_trigger_columns(n_timepoints)

    if haskey(ft, "cfg") && haskey(ft["cfg"], "event") && !isempty(ft["cfg"]["event"])

        events = ft["cfg"]["event"]

        # TODO: we have only tested this for our bdf data that is read into matlab
        # FieldTrip events are stored as struct-of-arrays (like EEGLAB)
        sample_idx = findfirst(==("sample"), events.names)
        value_idx = findfirst(==("value"), events.names)
        type_idx = findfirst(==("type"), events.names)

        if sample_idx !== nothing
            samples = events.values[sample_idx]
            @info "Extracting $(length(samples)) events from cfg.event"

            values = value_idx !== nothing ? events.values[value_idx] : nothing
            types = type_idx !== nothing ? events.values[type_idx] : nothing

            for i in eachindex(samples)
                sample = Int(round(samples[i]))

                if 1 <= sample <= n_timepoints
                    # Extract trigger value
                    if values !== nothing && i <= length(values)
                        val = values[i]
                        if val isa AbstractArray
                            !isempty(val) && (trigger_vec[sample] = Int(first(val)))
                        elseif !isnan(val)
                            trigger_vec[sample] = Int(val)
                        end
                    end

                    # Extract trigger type/info
                    if types !== nothing && i <= length(types)
                        typ = types[i]
                        if typ isa AbstractArray
                            !isempty(typ) && (trigger_info_vec[sample] = string(first(typ)))
                        else
                            trigger_info_vec[sample] = string(typ)
                        end
                    end
                end
            end
        end
    end

    insertcols!(df, 2, :triggers => trigger_vec)
    insertcols!(df, 3, :trigger_info => trigger_info_vec)

    # Create ContinuousData
    _, analysis_info = _create_metadata(filepath)
    continuous_data = ContinuousData(filepath, df, layout, sample_rate, analysis_info)

    @info "Successfully converted to ContinuousData"
    return continuous_data
end


"""
Convert FieldTrip epoched data to EegFun EpochData.
"""
function _fieldtrip_to_epochdata(ft::Dict, filepath::String, layout::Layout, sample_rate::Int)
    # Extract channel labels
    ch_labels = vec(ft["label"])
    ch_names = Symbol.(ch_labels)

    n_trials = length(ft["trial"])
    n_channels, n_timepoints = size(ft["trial"][1])

    @info "Data: $n_channels channels × $n_timepoints timepoints × $n_trials trials"

    # Extract time vector (handles both cell array and shared matrix formats)
    times = _extract_time_vector(ft["time"])

    # Extract trial info (trigger codes)
    if haskey(ft, "trialinfo") && !isempty(ft["trialinfo"])
        trial_triggers = vec(ft["trialinfo"])
    else
        trial_triggers = fill(0, n_trials)
    end

    # Convert trials to Vector{DataFrame}
    epoch_dataframes = Vector{DataFrame}(undef, n_trials)

    for trial_idx = 1:n_trials
        # Create DataFrame with transposed data and time column
        df = DataFrame(ft["trial"][trial_idx]', ch_names)
        insertcols!(df, 1, :time => times)

        # Create trigger vectors with trigger at time=0 (stimulus onset)
        trigger_vec, trigger_info_vec = _create_trigger_columns(n_timepoints)

        time_zero_idx = argmin(abs.(times))

        # Set trigger at time=0
        if trial_idx <= length(trial_triggers)
            trigger_code = Int(round(trial_triggers[trial_idx]))
            if trigger_code != 0
                trigger_vec[time_zero_idx] = trigger_code
                trigger_info_vec[time_zero_idx] = string(trigger_code)
            end
        end

        insertcols!(df, 2, :triggers => trigger_vec)
        insertcols!(df, 3, :trigger_info => trigger_info_vec)

        epoch_dataframes[trial_idx] = df
    end

    # Create EpochData
    condition_name, analysis_info = _create_metadata(filepath)
    epoch_data = EpochData(
        filepath,
        1,  # condition number (here, this is essentially a placeholder)
        condition_name,
        epoch_dataframes,
        layout,
        sample_rate,
        analysis_info,
    )

    @info "Successfully converted to EpochData with $n_trials epochs"
    return epoch_data
end

"""
Convert FieldTrip ERP data to EegFun ErpData.
"""
function _fieldtrip_to_erpdata(ft::Dict, filepath::String, layout::Layout, sample_rate::Int)
    # Extract channel labels and data
    ch_labels = vec(ft["label"])
    data_matrix = ft["avg"]  # [channels × timepoints]
    time_vec = vec(ft["time"])

    # Remove Status channel if present
    ch_labels, data_matrix = _remove_status_channel(ch_labels, data_matrix)

    ch_names = Symbol.(ch_labels)
    n_channels, n_timepoints = size(data_matrix)

    @info "Data: $n_channels channels × $n_timepoints timepoints"

    # Create DataFrame with time and trigger columns
    df = DataFrame(data_matrix', ch_names)
    insertcols!(df, 1, :time => time_vec)

    trigger_vec, trigger_info_vec = _create_trigger_columns(n_timepoints)
    insertcols!(df, 2, :triggers => trigger_vec)
    insertcols!(df, 3, :trigger_info => trigger_info_vec)

    # Create ErpData
    condition_name, analysis_info = _create_metadata(filepath)
    erp_data = ErpData(
        filepath,
        1,  # condition number
        condition_name,
        df,
        layout,
        sample_rate,
        analysis_info,
        1,  # n_epochs (unknown for FieldTrip ERP, default to 1)
    )

    @info "Successfully converted to ErpData"
    return erp_data
end
