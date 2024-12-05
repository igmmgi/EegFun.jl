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
search_sequence(array, sequence::Int) = intersect(findall(array .== sequence), findall(diff(vcat(0, array)) .>= 1))

find_idx_range(time, start_time, end_time) = findmin(abs.(time .- start_time))[2]:findmin(abs.(time .- end_time))[2]
find_idx_range(time, limits) = find_idx_range(time, limits[1], limits[end])
find_idx_start_end(time, start_time, end_time) = findmin(abs.(time .- start_time))[2], findmin(abs.(time .- end_time))[2]
find_idx_start_end(time, limits) = findmin(abs.(time .- limits[1]))[2], findmin(abs.(time .- limits[end]))[2]


function extract_epochs(dat::ContinuousData, trigger_sequence, start_time, end_time; zero_position=1)

  # find t==0 positions
  zero_idx = search_sequence(dat.data.triggers, trigger_sequence) .+ (zero_position - 1)
  isempty(zero_idx) && error("Trigger sequence not found!")

  # keep original sample index
  insertcols!(dat.data, 3, :sample => 1:nrow(dat.data))

  # find number of samples pre/post epoch t = 0 position
  n_pre, n_post = find_idx_start_end(dat.data.time, abs(start_time), abs(end_time))
  pre_idx = zero_idx .- n_pre .+ 1
  post_idx = zero_idx .+ n_post .- 1

  # extract and create array of dataframes
  epochs = []
  for (epoch, (pre, zero, post)) in enumerate(zip(pre_idx, zero_idx, post_idx))
    df = DataFrame(dat.data[pre:post, :])
    df.time = df.time .- dat.data.time[zero]
    insertcols!(df, 3, :epoch => epoch)
    push!(epochs, df)
  end

  return EpochData(epochs, dat.layout, dat.sample_rate)

end

function average_epochs(dat::EpochData)
  erp = combine(groupby(reduce(vcat, dat.data), :time), Not([:time, :events, :epoch, :sample]) .=> mean .=> Not([:time, :events, :epoch, :sample]))
  return ErpData(erp, dat.layout, dat.sample_rate)
end


