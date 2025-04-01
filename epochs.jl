"""
search_sequence(array, sequence)
Return index of a sequence within an array.
### Examples:
```julia
idx = search_sequence([1, 2, 3, 4, 2, 3, 4, 2], 4)
idx = search_sequence([1, 2, 3, 4, 2, 3, 4, 2], [2, 3])
idx = search_sequence([1, 2, 3, 4, 2, 3, 4, 2], [2 99])
"""
function search_sequence(array, sequence::Array{Int})

    idx_start_positions = search_sequence(array, sequence[1])

    # sequence of values
    idx_positions = []
    for idx in idx_start_positions
        good_sequence = true
        for seq = 1:(length(sequence)-1)
            if ((idx + seq) > length(array)) || array[idx+seq] != sequence[seq+1]
                good_sequence = false
                break
            end
        end
        if good_sequence
            push!(idx_positions, idx)
        end
    end

    return idx_positions

end



function extract_epochs(dat::ContinuousData, condition, trigger_sequence, start_time, end_time; zero_position = 1)

    df = deepcopy(dat.data)

    # find t==0 positions
    zero_idx = search_sequence(df.triggers, trigger_sequence) .+ (zero_position - 1)
    isempty(zero_idx) && error("Trigger sequence not found!")

    # find number of samples pre/post epoch t = 0 position
    n_pre, n_post = find_idx_start_end(df.time, abs(start_time), abs(end_time))
    pre_idx = zero_idx .- n_pre .+ 1
    post_idx = zero_idx .+ n_post .- 1

    # extract and create array of dataframes
    epochs = []
    for (epoch, (pre, zero, post)) in enumerate(zip(pre_idx, zero_idx, post_idx))
        epoch_df = DataFrame(df[pre:post, :])
        epoch_df.time = epoch_df.time .- df.time[zero]
        insertcols!(epoch_df, 4, :condition => condition)
        insertcols!(epoch_df, 5, :epoch => epoch)
        push!(epochs, epoch_df)
    end

    return EpochData(epochs, dat.layout, dat.sample_rate, dat.analysis_info)

end

function remove_bad_epochs(dat::EpochData)
    dat_out = deepcopy(dat)
    dat_out.data = filter(x -> !any(x[!, :is_extreme]), dat_out.data)
    @info "Epochs remaining: $(length(dat_out.data)) from $(length(dat.data)) epochs"
    return dat_out
end

function average_epochs(dat::EpochData)
    erp = combine(
        groupby(reduce(vcat, dat.data), :time),
        Not([:time, :triggers, :epoch, :sample]) .=> mean .=> Not([:time, :triggers, :epoch, :sample]),
    )
    return ErpData(erp, dat.layout, dat.sample_rate, dat.analysis_info)
end
