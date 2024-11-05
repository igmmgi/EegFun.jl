function create_eeg_dataframe(data::BioSemiBDF.BioSemiData)
  return hcat(DataFrame(time=data.time, events=data.triggers.raw),
    DataFrame(data.data, Symbol.(data.header.channel_labels[1:end-1])))
end

function create_eeg_dataframe(dat::BioSemiBDF.BioSemiData, layout_file_name::String)
  return ContinuousData(create_eeg_dataframe(dat), DataFrame(CSV.File(layout_file_name)), dat.header.sample_rate[1])
end

function create_eeg_dataframe(dat::BioSemiBDF.BioSemiData, layout::DataFrame)
  return ContinuousData(create_eeg_dataframe(dat), layout, dat.header.sample_rate[1])
end





function channel_summary(dat, channel_labels)
  out_stats = OrderedDict(OrderedDict())
  for col in names(dat.data)
    if col in channel_labels
      out_stats[Symbol(col)] = OrderedDict(:min => minimum(dat.data[!, col]),
        :max => maximum(dat.data[!, col]),
        :range => datarange(dat.data[!, col]),
        :std => std(dat.data[!, col]),
        :mad => mad(dat.data[!, col]))
    end
  end
  return out_stats
end

function correlation_matrix(dat)
  return [dat.layout.label DataFrame(cor(Matrix(dat.data[!, dat.layout.label])), dat.layout.label)]
end

function detect_eog_onsets!(dat, criterion, channel_in, channel_out)
  step_size = div(dat.sample_rate, 20)
  eog_diff = diff(dat.data[!, channel_in][1:step_size:end])
  eog_idx = findall(abs.(eog_diff) .>= criterion)
  eog_idx = eog_idx[(diff([0; eog_idx]).>2)] .* step_size
  dat.data[!, channel_out] .= false
  dat.data[eog_idx, channel_out] .= true
end

function is_extreme_value(dat::DataFrame, columns, criterion)
  return any(x -> abs.(x) .>= criterion, Matrix(dat[!, columns]), dims=2)
end

function n_extreme_value(dat::DataFrame, columns, criterion)
  return sum(x -> abs.(x) .>= criterion, Matrix(dat[!, columns]), dims=2)
end


