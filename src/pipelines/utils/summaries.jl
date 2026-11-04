function _file_percentage_summary(files::AbstractVector{S}, pcts::AbstractVector{<:Real}) where {S}
    file_idx = Dict{S,Int}()
    unique_files = S[]
    sums = Float64[]
    counts = Int[]
    @inbounds for i in eachindex(files, pcts)
        f = files[i]
        idx = get(file_idx, f, 0)
        if idx == 0
            push!(unique_files, f)
            push!(sums, Float64(pcts[i]))
            push!(counts, 1)
            file_idx[f] = length(unique_files)
        else
            sums[idx] += Float64(pcts[i])
            counts[idx] += 1
        end
    end
    @inbounds for k in eachindex(sums)
        sums[k] = round(sums[k] / counts[k]; digits = 1)
    end
    return DataFrame(AbstractVector[unique_files, sums], [:file, :percentage]; copycols = false)
end

function _epoch_and_file_summary(all_epoch_counts::Vector{DataFrame})
    epoch_summary = reduce(vcat, all_epoch_counts)
    file_summary = _file_percentage_summary(epoch_summary.file, epoch_summary.percentage)
    return epoch_summary, file_summary
end

function _merge_jld2_summary(path::String, new_df::DataFrame, new_files::Set)
    merged = if isfile(path)
        existing = JLD2.load(path, "data")
        keep = [!in(f, new_files) for f in existing.file]
        vcat(existing[keep, :], new_df)
    else
        new_df
    end
    return sort(merged, :file; by = _natural_sort_key)
end

function _merge_summaries(new_epoch_summary::DataFrame, new_file_summary::DataFrame, output_directory::String)
    new_files = Set(new_epoch_summary.file)
    merged_epoch = _merge_jld2_summary(joinpath(output_directory, "epoch_summary.jld2"), new_epoch_summary, new_files)
    merged_file = _merge_jld2_summary(joinpath(output_directory, "file_summary.jld2"), new_file_summary, new_files)
    return merged_epoch, merged_file
end

"""
    compare_rejections(rejection1::Vector{EpochRejectionInfo}, rejection2::Vector{EpochRejectionInfo}) -> DataFrame

Compare epoch rejection counts before (`rejection1`) and after (`rejection2`) channel repair per condition.
"""
function compare_rejections(rejection1::Vector{EpochRejectionInfo}, rejection2::Vector{EpochRejectionInfo})
    n = min(length(rejection1), length(rejection2))
    conds = Vector{Int}(undef, n)
    cond_names = Vector{String}(undef, n)
    n_orig = Vector{Int}(undef, n)
    rej_before = Vector{Int}(undef, n)
    rej_after = Vector{Int}(undef, n)
    rep_ok = Vector{Int}(undef, n)
    pct_rep = Vector{Float64}(undef, n)
    pct_kept = Vector{Float64}(undef, n)

    @inbounds for idx = 1:n
        info1 = rejection1[idx]
        info2 = rejection2[idx]
        s1 = unique_epochs(info1)
        s2 = unique_epochs(info2)

        n_original = info1.info.n
        n1 = length(s1)
        n2 = length(s2)
        s2_set = Set{Int}(s2)
        n_still_bad = count(in(s2_set), s1)
        n_repaired = n1 - n_still_bad

        conds[idx] = info1.info.number
        cond_names[idx] = info1.info.name
        n_orig[idx] = n_original
        rej_before[idx] = n1
        rej_after[idx] = n2
        rep_ok[idx] = n_repaired
        pct_rep[idx] = n1 > 0 ? round(100.0 * n_repaired / n1; digits = 1) : 0.0
        pct_kept[idx] = round(100.0 * (n_original - n2) / n_original; digits = 1)
    end

    return DataFrame(
        AbstractVector[conds, cond_names, n_orig, rej_before, rej_after, rep_ok, pct_rep, pct_kept],
        [
            :condition,
            :condition_name,
            :n_original,
            :rejected_before_repair,
            :rejected_after_repair,
            :successfully_repaired,
            :pct_repaired,
            :pct_final_kept,
        ];
        copycols = false,
    )
end

"""
    summarize_electrode_repairs(continuous_repairs::Vector{ContinuousRepairInfo}) -> DataFrame
    summarize_electrode_repairs(file_pattern::String; input_dir::String = pwd()) -> Tuple{DataFrame,DataFrame}

Summarize continuous-level electrode repairs across participants (either directly from a `Vector{ContinuousRepairInfo}`
or loaded from matching `*_artifact_info.jld2` files in `input_dir`).
"""
function summarize_electrode_repairs(continuous_repairs::Vector{ContinuousRepairInfo})
    counts = Dict{Symbol,Int}()
    last_seen = Dict{Symbol,Int}()

    for (participant_idx, repair_info) in enumerate(continuous_repairs)
        for ch in repair_info.repaired
            if get(last_seen, ch, 0) != participant_idx
                last_seen[ch] = participant_idx
                counts[ch] = get(counts, ch, 0) + 1
            end
        end
    end

    isempty(counts) && return DataFrame(electrode = Symbol[], n_participants = Int[])
    electrodes = collect(keys(counts))
    n_parts = [counts[ch] for ch in electrodes]
    perm = sortperm(n_parts; rev = true)
    return DataFrame(AbstractVector[electrodes[perm], n_parts[perm]], [:electrode, :n_participants]; copycols = false)
end

function summarize_electrode_repairs(file_pattern::String; input_dir::String = pwd())
    files = sort!(_find_batch_files(file_pattern, input_dir); by = _natural_sort_key)
    if isempty(files)
        @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
        return DataFrame(file = String[], repaired = String[], skipped = String[]), DataFrame(electrode = Symbol[], n_participants = Int[])
    end

    files_col = String[]
    rep_col = String[]
    skip_col = String[]
    continuous_repairs = ContinuousRepairInfo[]

    for file in files
        file_path = joinpath(input_dir, file)
        isfile(file_path) || continue
        artifact_info = try
            read_data(file_path)
        catch e
            @minimal_warning "Failed to load artifact info from $file_path: $(sprint(showerror, e))"
            continue
        end

        file_repaired = Symbol[]
        file_skipped = Symbol[]
        if !isnothing(artifact_info)
            for repair_info in artifact_info.continuous_repairs
                append!(file_repaired, repair_info.repaired)
                append!(file_skipped, repair_info.skipped)
                !isempty(repair_info.repaired) && push!(continuous_repairs, repair_info)
            end
        end
        push!(files_col, replace(file, "_artifact_info.jld2" => ""))
        push!(rep_col, isempty(file_repaired) ? "none" : join(sort!(string.(unique!(file_repaired))), ", "))
        push!(skip_col, isempty(file_skipped) ? "" : join(sort!(string.(unique!(file_skipped))), ", "))
    end

    per_file_df = DataFrame(AbstractVector[files_col, rep_col, skip_col], [:file, :repaired, :skipped]; copycols = false)
    return per_file_df, summarize_electrode_repairs(continuous_repairs)
end

function _empty_ica_summary()
    return (
        DataFrame(
            file = String[],
            total_components = Int[],
            vEOG = Int[],
            hEOG = Int[],
            ECG = Int[],
            line_noise = Int[],
            channel_noise = Int[],
        ),
        DataFrame(component_type = String[], avg_per_participant = Float64[]),
    )
end

"""
    summarize_ica_components(ica_components::Vector{ArtifactComponents}) -> Tuple{DataFrame,DataFrame}
    summarize_ica_components(file_pattern::String; input_dir::String = pwd()) -> Tuple{DataFrame,DataFrame}

Summarize removed ICA artifact components per participant and on average across participants.
"""
function summarize_ica_components(ica_components::Vector{ArtifactComponents})
    n = length(ica_components)
    n == 0 && return _empty_ica_summary()

    files_col = fill("", n)
    total_col = Vector{Int}(undef, n)
    veog_col = Vector{Int}(undef, n)
    heog_col = Vector{Int}(undef, n)
    ecg_col = Vector{Int}(undef, n)
    ln_col = Vector{Int}(undef, n)
    cn_col = Vector{Int}(undef, n)

    @inbounds for i = 1:n
        ca = ica_components[i]
        veog_col[i] = length(get(ca.eog, :vEOG, Int[]))
        heog_col[i] = length(get(ca.eog, :hEOG, Int[]))
        ecg_col[i] = length(ca.ecg)
        ln_col[i] = length(ca.line_noise)
        cn_col[i] = length(ca.channel_noise)
        total_col[i] = length(unique(get_all_ica_components(ca)))
    end

    per_file_df = DataFrame(
        AbstractVector[files_col, total_col, veog_col, heog_col, ecg_col, ln_col, cn_col],
        [:file, :total_components, :vEOG, :hEOG, :ECG, :line_noise, :channel_noise];
        copycols = false,
    )

    avg_vals = Float64[
        round(mean(veog_col); digits = 2),
        round(mean(heog_col); digits = 2),
        round(mean(ecg_col); digits = 2),
        round(mean(ln_col); digits = 2),
        round(mean(cn_col); digits = 2),
        round(mean(total_col); digits = 2),
    ]
    avg_df = DataFrame(
        AbstractVector[["vEOG", "hEOG", "ECG", "line_noise", "channel_noise", "total"], avg_vals],
        [:component_type, :avg_per_participant];
        copycols = false,
    )
    return per_file_df, avg_df
end

function summarize_ica_components(file_pattern::String; input_dir::String = pwd())
    files = _find_batch_files(file_pattern, input_dir)
    if isempty(files)
        @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
        return _empty_ica_summary()
    end

    ica_components = ArtifactComponents[]
    filenames = String[]
    for file in files
        file_path = joinpath(input_dir, file)
        isfile(file_path) || continue
        artifact_info = try
            read_data(file_path)
        catch e
            @minimal_warning "Failed to load artifact info from $file_path: $(sprint(showerror, e))"
            continue
        end
        if !isnothing(artifact_info) && !isnothing(artifact_info.ica_components)
            push!(ica_components, artifact_info.ica_components)
            push!(filenames, replace(file, "_artifact_info.jld2" => ""))
        end
    end

    per_file_df, avg_df = summarize_ica_components(ica_components)
    for i in eachindex(filenames)
        per_file_df.file[i] = filenames[i]
    end
    sort!(per_file_df, :file; by = _natural_sort_key)
    return per_file_df, avg_df
end
