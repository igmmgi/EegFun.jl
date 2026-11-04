"""
    read_fieldtrip_csv(data_dir::String;
                       file::String = "csv_data",
                       condition::Int = 1,
                       condition_name::String = "csv_data",
                       fsample::Union{Int,Nothing} = nothing) -> EpochData

Load FieldTrip epochs exported as `epochs_data.csv` (columns `[trial, time, ch1, ch2, ...]`) in `data_dir`
and return an `EpochData` object.
"""
function read_fieldtrip_csv(
    data_dir::String;
    file::String = "csv_data",
    condition::Int = 1,
    condition_name::String = "csv_data",
    fsample::Union{Int,Nothing} = nothing,
)
    data_file = joinpath(data_dir, "epochs_data.csv")
    isfile(data_file) || @minimal_error("Data file not found: $data_file")
    all_df = CSV.read(data_file, DataFrame)

    actual_cols = names(all_df)
    if length(actual_cols) < 2 || actual_cols[1] != "trial" || actual_cols[2] != "time"
        @minimal_error("Expected first two columns to be 'trial' and 'time', got: $(actual_cols[1:min(2, length(actual_cols))])")
    end

    ch_labels = [Symbol(col) for col in @view actual_cols[3:end]]
    grouped = groupby(all_df, :trial)
    n_trials = length(grouped)

    if isnothing(fsample)
        if n_trials > 0 && length(grouped[1].time) > 1
            fsample = round(Int, 1.0 / median(diff(sort(collect(grouped[1].time)))))
            @info "Estimated sample rate from time column: $fsample Hz"
        else
            @minimal_error("Cannot estimate sample rate. Please provide fsample parameter.")
        end
    end

    n_ch = length(ch_labels)
    layout = Layout(DataFrame(:label => ch_labels, :inc => zeros(n_ch), :azi => zeros(n_ch)), nothing, nothing, nothing)
    analysis_info = AnalysisInfo(sample_rate = Float64(fsample))

    col_names = vcat([:time, :sample, :epoch], ch_labels)
    trial_dfs = Vector{DataFrame}(undef, n_trials)
    for (trial_idx, trial_data) in enumerate(grouped)
        n_rows = nrow(trial_data)
        cols = Vector{Any}(undef, 3 + n_ch)
        cols[1] = collect(Float64, trial_data.time)
        cols[2] = collect(1:n_rows)
        cols[3] = fill(trial_idx, n_rows)
        for (j, ch) in enumerate(ch_labels)
            cols[3+j] = collect(Float64, trial_data[!, ch])
        end
        trial_dfs[trial_idx] = DataFrame(cols, col_names; copycols = false)
    end

    return EpochData(file, condition, condition_name, trial_dfs, layout, fsample, analysis_info)
end
