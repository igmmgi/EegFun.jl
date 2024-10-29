###############################################################
function _apply_baseline!(dat::DataFrame, channel_labels, baseline_interval)
  for col in names(dat)
    if col in channel_labels
      dat[:, col] .-= mean(dat[baseline_interval[1]:baseline_interval[2], col])
    end
  end
end

function baseline!(dat::DataFrame, channel_labels, baseline_interval)
  if isempty(baseline_interval)
    baseline_interval = [dat.time[1], dat.time[end]]
  end
  baseline_interval = find_idx_start_end(dat.time, baseline_interval[1], baseline_interval[2])
  _apply_baseline!(dat, channel_labels, baseline_interval)
end

function baseline!(dat::Union{ContinuousData,ErpData}, channel_labels, baseline_interval)
  baseline!(dat.data, channel_labels, baseline_interval)
end

function baseline!(dat::Union{ContinuousData,ErpData}, baseline_interval)
  baseline!(dat.data, dat.layout.label, baseline_interval)
end

function baseline!(dat::EpochData, channel_labels, baseline_interval)
  for epoch in eachindex(dat.data)
    baseline!(dat.data[epoch], channel_labels, baseline_interval)
  end
end

function baseline!(dat::EpochData, baseline_interval)
  for epoch in eachindex(dat.data)
    baseline!(dat.data[epoch], dat.layout.label, baseline_interval)
  end
end

function baseline(dat::DataFrame, channel_labels, baseline_interval)
  dat_out = deepcopy(dat)
  baseline!(dat_out, channel_labels, baseline_interval)
  return dat_out
end

function baseline(dat::Union{ContinuousData,ErpData}, channel_labels, baseline_interval)
  dat_out = deepcopy(dat)
  baseline!(dat_out.data, channel_labels, baseline_interval)
  return dat_out
end

function baseline(dat::Union{ContinuousData,ErpData}, baseline_interval)
  dat_out = deepcopy(dat)
  baseline!(dat_out.data, dat_out.layout.label, baseline_interval)
  return dat_out
end

function baseline(dat::EpochData, channel_labels, baseline_interval)
  dat_out = deepcopy(dat)
  baseline!(dat_out, channel_labels, baseline_interval)
  return dat_out
end

function baseline(dat::EpochData, baseline_interval)
  dat_out = deepcopy(dat)
  baseline!(dat_out, dat_out.layout.label, baseline_interval)
  return dat_out
end


