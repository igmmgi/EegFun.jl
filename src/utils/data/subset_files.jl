_default_subset_output_dir(input_dir::String, pattern::String) = joinpath(input_dir, "subset_$(clean_pattern(pattern))")

function _process_subset_file(filepath::String, output_path::String; kwargs...)
    fname = basename(filepath)
    dat = read_data(filepath)
    isnothing(dat) && return BatchResult(false, fname, "No data found in file")

    subset_data = subset(dat; kwargs...)
    jldsave(output_path; data = subset_data, format_version = EEGFUN_FORMAT_VERSION)

    if subset_data isa Vector && !isempty(subset_data) && hasproperty(subset_data[1], :data)
        n_items = sum(length(cond.data) for cond in subset_data)
        return BatchResult(true, fname, "Processed subset. Elements remaining: $n_items")
    elseif hasproperty(subset_data, :data) && subset_data.data isa Vector
        return BatchResult(true, fname, "Processed subset. Elements remaining: $(length(subset_data.data))")
    else
        return BatchResult(true, fname, "Processed subset.")
    end
end

"""
    subset(file_pattern::String; input_dir = pwd(), participant_selection = participants(), output_dir = nothing, kwargs...)

Batch-subset JLD2 files matching `file_pattern` using any `subset` keyword arguments.
"""
function subset(
    file_pattern::String;
    input_dir::String = pwd(),
    participant_selection = participants(),
    output_dir::Union{String,Nothing} = nothing,
    kwargs...,
)
    _with_batch_logging("subset", "subset.log", input_dir, output_dir) do ctx
        out_dir = _prepare_output_dir!(ctx, _default_subset_output_dir(input_dir, file_pattern))
        process_fn = (in_path, out_path) -> _process_subset_file(in_path, out_path; kwargs...)
        batch_process(process_fn, file_pattern, input_dir, out_dir, participant_selection, "Subsetting parameters")
    end
end
