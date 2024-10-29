# channel difference
function calculate_channel_difference(dat::DataFrame, channels1::Vector, channels2::Vector)
  channels1 = sum([dat[!, c] for c in channels1]) / length(channels1)
  channels2 = sum([dat[!, c] for c in channels2]) / length(channels2)
  return channels1 .- channels2
end

function calculate_channel_difference(dat::DataFrame, channels1, channels2)
  return calculate_channel_difference(dat, [channels1], [channels2])
end

function diff_channel!(dat::DataFrame, channels1::Vector, channels2::Vector, difference_label)
  dat[:, difference_label] = calculate_channel_difference(dat, channels1, channels2)
end

function diff_channel!(dat::DataFrame, channels1, channels2, difference_label)
  diff_channel!(dat, [channels1], [channels2], difference_label)
end

function diff_channel!(dat::Union{ContinuousData,ErpData}, channels1::Vector, channels2::Vector, difference_label)
  dat.data[:, difference_label] = calculate_channel_difference(dat.data, channels1, channels2)
end

function diff_channel!(dat::Union{ContinuousData,ErpData}, channels1, channels2, difference_label)
  diff_channel!(dat, [channels1], [channels2], difference_label)
end

function diff_channel!(dat::EpochData, channels1::Vector, channels2::Vector, difference_label)
  for epoch in eachindex(dat.data)
    diff_channel!(dat.data[epoch], channels1, channels2, difference_label)
  end
end

function diff_channel!(dat::EpochData, channels1, channels2, difference_label)
  diff_channel!(dat, [channels1], [channels2], difference_label)
end


