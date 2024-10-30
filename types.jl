mutable struct ContinuousData
  data::DataFrame
  layout::DataFrame
  sample_rate::Int
end

function Base.show(io::IO, dat::ContinuousData)
  println(io, "Data: Rows ($(nrow(dat.data))) x Columns ($(ncol(dat.data))) ")
  println(io, "Labels: ", join(names(dat.data), ", "))
  println(io, "Sample Rate: ", dat.sample_rate)
end

mutable struct EpochData
  data::Array{DataFrame}
  layout::DataFrame
  sample_rate::Int
end

function Base.show(io::IO, dat::EpochData)
  println(io, "Number of trials: ", length(dat.data))
  println(io, "Timepoints ($(nrow(dat.data[1]))) x Channels ($(ncol(dat.data[1]))) ")
  println(io, "Channel Labels: ", join(names(dat.data[1]), ", "))
  println(io, "Sample Rate: ", dat.sample_rate)
end

mutable struct ErpData
  data::DataFrame
  layout::DataFrame
  sample_rate::Int
end

function Base.show(io::IO, dat::ErpData)
  println(io, "Data: Rows ($(nrow(dat.data))) x Columns ($(ncol(dat.data))) ")
  println(io, "Labels: ", join(names(dat.data), ", "))
  println(io, "Sample Rate: ", dat.sample_rate)
end

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




