_default_grand_average_output_dir(input_dir::String, pattern::String) = joinpath(input_dir, "grand_average_$(clean_pattern(pattern))")

function _create_grand_average(erps::AbstractVector{<:ErpData}, cond_num::Int)
    isempty(erps) && @minimal_error("Cannot create grand average from empty ERP list")
    _have_same_structure(erps) || @minimal_error("ERPs have inconsistent structure")

    first_erp = erps[1]
    meta_set = Set{Symbol}(meta_labels(first_erp))
    common_chs = Set{Symbol}(c for c in propertynames(first_erp.data) if c ∉ meta_set)
    for m = 2:length(erps)
        intersect!(common_chs, propertynames(erps[m].data))
    end

    df1 = first_erp.data
    col_names = propertynames(df1)
    n_points = nrow(df1)
    n_erps = length(erps)
    inv_n = 1.0 / n_erps
    out_cols = Vector{Any}(undef, length(col_names))

    for (j, col) in enumerate(col_names)
        if col in common_chs
            out = Vector{Float64}(undef, n_points)
            c1 = df1[!, j]
            if c1 isa Vector{Float64}
                v1 = c1::Vector{Float64}
                @inbounds @simd for i = 1:n_points
                    out[i] = v1[i]
                end
            else
                @inbounds @simd for i = 1:n_points
                    out[i] = Float64(c1[i])
                end
            end
            for m = 2:n_erps
                cm = erps[m].data[!, col]
                if cm isa Vector{Float64}
                    vm = cm::Vector{Float64}
                    @inbounds @simd for i = 1:n_points
                        out[i] += vm[i]
                    end
                else
                    @inbounds @simd for i = 1:n_points
                        out[i] += Float64(cm[i])
                    end
                end
            end
            if n_erps > 1
                @inbounds @simd for i = 1:n_points
                    out[i] *= inv_n
                end
            end
            out_cols[j] = out
        else
            out_cols[j] = copy(df1[!, j])
        end
    end

    grand_avg_data = DataFrame(out_cols, col_names; copycols = false)
    total_epochs = sum(erp.n_epochs for erp in erps)
    return ErpData(
        "grand_avg",
        cond_num,
        "grand_avg_$(first_erp.condition_name)",
        grand_avg_data,
        copy(first_erp.layout),
        first_erp.sample_rate,
        copy(first_erp.analysis_info),
        total_epochs,
    )
end

function _create_grand_average(tfs::AbstractVector{<:TimeFreqData}, cond_num::Int)
    isempty(tfs) && @minimal_error("Cannot create grand average from empty TF list")
    _have_same_structure(tfs) || @minimal_error("TF data for condition $cond_num have inconsistent structure")

    ref_tf = tfs[1]
    electrodes = channel_labels(ref_tf)
    frequencies = Float64.(sort(unique(ref_tf.data_power.freq)))
    time_points = Float64.(sort(unique(ref_tf.data_power.time)))

    data_4d = _extract_tf_array(tfs, electrodes, frequencies, time_points)
    return _create_tf_grand_average(ref_tf, data_4d, electrodes, frequencies, time_points, cond_num)
end

_log_creating_grand_average(::Type{<:ErpData}, cond_num::Int, n::Int) =
    @info "Creating grand average for condition $cond_num (n=$n participants)"
_log_creating_grand_average(::Type{<:TimeFreqData}, cond_num::Int, n::Int) =
    @info "Creating TF grand average for condition $cond_num (n=$n participants)"
_log_created_grand_average(ga::ErpData) = @info "  ✓ Created grand average with $(nrow(ga.data)) time points"
function _log_created_grand_average(ga::TimeFreqData)
    n_freqs = length(unique(ga.data_power.freq))
    n_times = length(unique(ga.data_power.time))
    @info "  ✓ Created TF grand average: $(n_freqs) freqs × $(n_times) times"
end

function _create_all_grand_averages(data_by_condition::AbstractDict{<:Integer,<:AbstractVector{T}}) where {T<:Union{ErpData,TimeFreqData}}
    grand_averages = T[]
    for cond_num in sort(collect(keys(data_by_condition)))
        condition_data = data_by_condition[cond_num]
        _log_creating_grand_average(T, cond_num, length(condition_data))
        if length(condition_data) < 2
            @minimal_warning "Only $(length(condition_data)) participant(s) for condition $cond_num. Skipping."
            continue
        end
        grand_avg = _create_grand_average(condition_data, cond_num)
        push!(grand_averages, grand_avg)
        _log_created_grand_average(grand_avg)
    end
    return grand_averages
end

function _select_conditions_by_field(items::AbstractVector, condition_selection)
    isempty(items) && return Int[]
    (isnothing(condition_selection) || condition_selection === Colon()) && return collect(1:length(items))

    cond_func = _to_condition_func(condition_selection)
    cond_func === Colon() && return collect(1:length(items))

    if all(x -> hasfield(typeof(x), :condition), items)
        cond_nums = Int[x.condition for x in items]
        if allunique(cond_nums) && all(>(0), cond_nums)
            if hasfield(typeof(cond_func), :condition_indices)
                req_ids = getfield(cond_func, :condition_indices)
                return Int[i for (i, c) in enumerate(cond_nums) if c in req_ids]
            elseif hasfield(typeof(cond_func), :f) && hasfield(typeof(getfield(cond_func, :f)), :condition_indices)
                req_ids = getfield(getfield(cond_func, :f), :condition_indices)
                return Int[i for (i, c) in enumerate(cond_nums) if !(c in req_ids)]
            end
        end
    end

    items isa AbstractVector{<:EegFunData} && return get_selected_conditions(items, cond_func)
    res = cond_func(1:length(items))
    return res isa AbstractVector{Bool} ? findall(res) : collect(Int, res)
end

"""
    grand_average(data::AbstractVector{<:Union{ErpData, TimeFreqData}}; condition_selection = conditions()) -> Vector
    grand_average(rsa_data_list::AbstractVector{RsaData}; compute_noise_ceiling::Bool = true) -> RsaData
    grand_average(decoded_list::AbstractVector{DecodedData}) -> DecodedData
    grand_average(file_pattern::String; input_dir = pwd(), participant_selection = participants(),
                  condition_selection = conditions(), output_dir = nothing)

Average ERP, time-frequency, RSA, or decoding results across participants (per condition where
applicable), or batch-process matching JLD2 files into `output_dir`.
"""
function grand_average(data::AbstractVector{T}; condition_selection = conditions()) where {T<:Union{ErpData,TimeFreqData}}
    isempty(data) && @minimal_error("Cannot create grand average from empty $(T <: ErpData ? "ERP" : "TF") list")

    data_by_condition = group_by_condition(data)
    all_cond_nums = collect(keys(data_by_condition))
    rep_items = T[first(data_by_condition[num]) for num in all_cond_nums]
    selected_indices = _select_conditions_by_field(rep_items, condition_selection)
    selected_cond_nums = all_cond_nums[selected_indices]

    filtered_groups = OrderedDict(num => data_by_condition[num] for num in selected_cond_nums)
    if isempty(filtered_groups)
        @minimal_warning "No conditions selected or found for grand averaging"
        return T[]
    end
    return _create_all_grand_averages(filtered_groups)
end

function grand_average(
    file_pattern::String;
    input_dir::String = pwd(),
    participant_selection = participants(),
    condition_selection = conditions(),
    output_dir::Union{String,Nothing} = nothing,
)
    _with_batch_logging("grand averaging", "grand_average.log", input_dir, output_dir) do ctx
        @log_call "grand_average"
        output_dir = _prepare_output_dir!(ctx, _default_grand_average_output_dir(input_dir, file_pattern))

        files = _find_batch_files(file_pattern, input_dir, participant_selection)
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end
        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"

        all_data = _read_all_data_core(EegData, files, input_dir)
        if isempty(all_data)
            @minimal_warning "No valid data found in any files"
            return nothing
        end

        T = typeof(all_data[1])
        typed_data = T[x for x in all_data]
        grand_averages = grand_average(typed_data; condition_selection = condition_selection)
        if isempty(grand_averages)
            @minimal_warning "No grand averages created"
            return nothing
        end

        output_path = joinpath(output_dir, "grand_average_$(clean_pattern(file_pattern)).jld2")
        jldsave(output_path; data = grand_averages, format_version = EEGFUN_FORMAT_VERSION)

        @info "Grand averaging complete! Created $(length(grand_averages)) grand averages"
        @info "Output saved to: $output_path"
    end
end
