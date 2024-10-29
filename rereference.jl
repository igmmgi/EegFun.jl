###############################################################
# re-rereference
function _apply_rereference!(dat::DataFrame, channel_labels, reference)
  for col in names(dat)
    if col in channel_labels
      dat[:, col] .-= reference
    end
  end
end

function rereference!(dat::DataFrame, channel_labels, reference_channels::Union{Int64,Vector{Int64}})
  reference_channels = channel_number_to_channel_label(channel_labels, reference_channels)
  reference = reduce(+, eachcol(dat[:, reference_channels])) ./ length(reference_channels)
  _apply_rereference!(dat, channel_labels, reference)
end

function rereference!(dat::DataFrame, channel_labels, reference_channels::Union{AbstractString,Symbol})
  rereference!(dat, channel_labels, [reference_channels])
end

function rereference!(dat::DataFrame, channel_labels, reference_channel::Union{Vector{<:AbstractString},Vector{Symbol}})
  reference = reduce(+, eachcol(dat[:, reference_channel])) ./ length(reference_channel)
  _apply_rereference!(dat, channel_labels, reference)
end

function rereference!(dat::Union{ContinuousData,ErpData}, channel_labels, reference_channel)
  rereference!(dat.data, channel_labels, reference_channel)
end

function rereference!(dat::Union{ContinuousData,ErpData}, reference_channel)
  rereference!(dat.data, dat.layout.label, reference_channel)
end

function rereference!(dat::EpochData, channel_labels, reference_channel)
  for epoch in eachindex(dat.data)
    rereference!(dat.data[epoch], channel_labels, reference_channel)
  end
end

function rereference!(dat::EpochData, reference_channel)
  for epoch in eachindex(dat.data)
    rereference!(dat.data[epoch], dat.layout.label, reference_channel)
  end
end

function rereference(dat::DataFrame, channel_labels, reference_channel::Union{Vector{<:AbstractString},Vector{Symbol}})
  dat_out = deepcopy(dat)
  reference = reduce(+, eachcol(dat[:, reference_channel])) ./ length(reference_channel)
  _apply_rereference!(dat_out, channel_labels, reference)
  return dat_out
end

function rereference(dat::DataFrame, channel_labels, reference_channels::Union{Int64,Vector{Int64}})
  dat_out = deepcopy(dat)
  rereference!(dat_out, channel_labels, reference_channels)
  return dat_out
end

function rereference(dat::DataFrame, channel_labels, reference_channel::Union{AbstractString,Symbol})
  dat_out = deepcopy(dat)
  reference_channel = [reference_channel]
  _apply_rereference!(dat_out, channel_labels, reference_channel)
  return dat_out
end

function rereference(dat::Union{ContinuousData,ErpData}, channel_labels, reference_channel)
  dat_out = deepcopy(dat)
  rereference!(dat_out.data, channel_labels, reference_channel)
  return dat_out
end

function rereference(dat::Union{ContinuousData,ErpData}, reference_channel)
  dat_out = deepcopy(dat)
  rereference!(dat_out.data, dat_out.layout.label, reference_channel)
  return dat_out
end

function rereference(dat::EpochData, channel_labels, reference_channel)
  dat_out = deepcopy(dat)
  rereference!(dat_out, channel_labels, reference_channel)
  return dat_out
end

function rereference(dat::EpochData, reference_channel)
  dat_out = deepcopy(dat)
  rereference!(dat_out, dat_out.layout.label, reference_channel)
  return dat_out
end


