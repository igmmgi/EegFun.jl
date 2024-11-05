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

