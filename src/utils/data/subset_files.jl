function _default_subset_output_dir(input_dir::String, pattern::String)
    joinpath(input_dir, "subset_$(clean_pattern(pattern))")
end

function _process_subset_file(filepath::String, output_path::String; kwargs...)
    filename = basename(filepath)
    dat = read_data(filepath)

    if isnothing(dat)
        return BatchResult(false, filename, "No data found in file")
    end

    # Exceptions are handled by _run_batch_operation
    subset_data = subset(dat; kwargs...)
    jldsave(output_path; data = subset_data)

    # Determine stats for logging
    if subset_data isa Vector && !isempty(subset_data) && hasproperty(subset_data[1], :data)
        n_items = sum(length(cond.data) for cond in subset_data)
        return BatchResult(true, filename, "Processed subset. Elements remaining: $(n_items)")
    elseif hasproperty(subset_data, :data) && subset_data.data isa Vector
        n_items = length(subset_data.data)
        return BatchResult(true, filename, "Processed subset. Elements remaining: $(n_items)")
    else
        return BatchResult(true, filename, "Processed subset.")
    end
end

"""
    subset(file_pattern::String; kwargs...)

Batch process JLD2 files that match `file_pattern`, applying the `subset` operation.
Supports all `subset` kwargs (`epoch_selection`, `channel_selection`, `sample_selection`, `interval_selection`).
"""
function subset(
    file_pattern::String;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    output_dir::Union{String,Nothing} = nothing,
    kwargs...,
)
    _with_batch_logging("subset", "subset.log", input_dir, output_dir) do ctx
        output_dir = _prepare_output_dir!(ctx, _default_subset_output_dir(input_dir, file_pattern))

        process_fn = (input_path, output_path) -> _process_subset_file(input_path, output_path; kwargs...)

        batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "Subsetting parameters")
    end
end
