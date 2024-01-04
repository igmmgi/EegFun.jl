using BioSemiBDF
using CSV
using DSP
using DataFrames
using GLMakie
using JLD2
using ScatteredInterpolation
using StatsBase


function polar2cartXY(layout::DataFrame)
  inc = layout[!, :inc] .* (pi / 180)
  azi = layout[!, :azi] .* (pi / 180)
  layout[!, "X2"] = inc .* cos.(azi)
  layout[!, "Y2"] = inc .* sin.(azi)
end

mutable struct ContinuousData
  data::DataFrame
  layout::DataFrame
  sample_rate::Int
end

function Base.show(io::IO, dat::ContinuousData)
  println(io, "Timepoints ($(nrow(dat.data))) x Channels ($(ncol(dat.data))) ")
  println(io, "Channel Labels: ", join(names(dat.data), ", "))
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
  println(io, "Timepoints ($(nrow(dat.data))) x Channels ($(ncol(dat.data))) ")
  println(io, "Channel Labels: ", join(names(dat.data), ", "))
  println(io, "Sample Rate: ", dat.sample_rate)
end


########################################################################

function create_dataframe(data::BioSemiBDF.BioSemiData)
  df1 = DataFrame(time=data.time, events=data.triggers.raw)
  df2 = DataFrame(data.data, :auto)
  rename!(df2, data.header.channel_labels[1:end-1])
  return hcat(df1, df2)
end

function eeg_data(dat::BioSemiBDF.BioSemiData, layout_file_name::String)
  return ContinuousData(create_dataframe(dat), DataFrame(CSV.File(layout_file_name)), dat.header.sample_rate[1])
end


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
search_sequence(array, sequence::Int) = findall(array .== sequence)


find_idx_range(time, t1, t2) = findmin(abs.(time .- t1))[2]:findmin(abs.(time .- t2))[2]
find_idx_range(time, limits) = find_idx_range(time, limits[1], limits[end])
find_idx_start_end(time, t1, t2) = findmin(abs.(time .- t1))[2], findmin(abs.(time .- t2))[2]
find_idx_start_end(time, limits) = findmin(abs.(time .- limits[1]))[2], findmin(abs.(time .- limits[end]))[2]


function extract_epochs(dat::ContinuousData, trigger_sequence, start_time, end_time; zero_position=1)

  zero_idx = search_sequence(dat.data.events, trigger_sequence) .+ (zero_position - 1)
  isempty(zero_idx) && error("Trigger sequence not found!")

  n_pre, n_post = find_idx_start_end(dat.data.time, abs(start_time), abs(end_time))
  pre_idx = zero_idx .- n_pre .- 1
  post_idx = zero_idx .+ n_post .- 1

  epochs = []
  for (pre, zero, post) in zip(pre_idx, zero_idx, post_idx)
    df = DataFrame(dat.data[pre:post, :])
    df.time = df.time .- dat.data.time[zero]
    push!(epochs, df)
  end

  return EpochData(epochs, dat.layout, dat.sample_rate)

end

function average_epochs(dat::EpochData)
  erp = combine(groupby(reduce(vcat, dat.data), :time), Not([:time, :events]) .=> mean .=> Not([:time, :events]))
  return ErpData(erp, dat.layout, dat.sample_rate)
end

dat = read_bdf("../Flank_C_3.bdf")
dat = eeg_data(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
epochs = extract_epochs(dat, 1, -0.5, 2)
erp = average_epochs(epochs)


# Filter functions
function highpass_filter!(dat::ContinuousData, freq, order)
  filter = digitalfilter(Highpass(freq, fs=dat.sample_rate), Butterworth(order))
  for col in names(dat.data)[3:end]
    dat.data[:, col] .= filtfilt(filter, dat.data[:, col])
  end
end

function highpass_filter(dat::ContinuousData, freq, order)
  dat_out = deepcopy(dat)
  highpass_filter!(dat_out, freq, order)
  return dat_out
end

# Filter functions
function highpass_filter!(dat::ContinuousData, freq, order)
  filter = digitalfilter(Highpass(freq, fs=dat.sample_rate), Butterworth(order))
  for col in names(dat.data)[3:end]
    dat.data[:, col] .= filtfilt(filter, dat.data[:, col])
  end
end

function highpass_filter!(dat::EpochData, freq, order)
  filter = digitalfilter(Highpass(freq, fs=dat.sample_rate), Butterworth(order))
  for epoch in eachindex(dat.data)
    for col in names(dat.data[epoch])[3:end]
      dat.data[epoch][:, col] .= filtfilt(filter, dat.data[epoch][:, col])
    end
  end
end

function highass_filter(dat::EpochData, freq, order)
  dat_out = deepcopy(dat)
  highpass_filter!(dat_out, freq, order)
  return dat_out
end

function lowpass_filter!(dat::ContinuousData, freq, order)
  filter = digitalfilter(Lowpass(freq, fs=dat.sample_rate), Butterworth(order))
  for col in names(dat.data)[3:end]
    dat.data[:, col] .= filtfilt(filter, dat.data[:, col])
  end
end

function lowpass_filter(dat::ContinuousData, freq, order)
  dat_out = deepcopy(dat)
  lowpass_filter!(dat_out, freq, order)
  return dat_out
end

function lowpass_filter!(dat::EpochData, freq, order)
  filter = digitalfilter(Lowpass(freq, fs=dat.sample_rate), Butterworth(order))
  for epoch in eachindex(dat.data)
    for col in names(dat.data[epoch])[3:end]
      dat.data[epoch][:, col] .= filtfilt(filter, dat.data[epoch][:, col])
    end
  end
end




function lowpass_filter(dat::EpochData, freq, order)
  dat_out = deepcopy(dat)
  lowpass_filter!(dat_out, freq, order)
  return dat_out
end

dat = read_bdf("../Flank_C_3.bdf")
dat = eeg_data(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
epochs = extract_epochs(dat, 1, -0.5, 2)
erp = average_epochs(epochs)




# re-reference functions
function rereference!(dat::ContinuousData, reference_channel)
  reference = reduce(+, eachcol(dat.data[:, reference_channel])) ./ length(reference_channel)
  for col in names(dat.data)[3:end]
    dat.data[:, col] .-= reference
  end
end
rereference!(dat::ContinuousData, reference_channel::Symbol) = rereference!(dat, [reference_channel])
rereference!(dat::ContinuousData, reference_channel::Vector{UnitRange{Int64}}) = rereference!(dat, reference_channel[1])

function rereference(dat::ContinuousData, reference_channel)
  dat_out = deepcopy(dat)
  rereference!(dat_out, reference_channel)
  return dat_out
end
rereference(dat::ContinuousData, reference_channel::Symbol) = rereference(dat, [reference_channel])
rereference(dat::ContinuousData, reference_channel::Vector{UnitRange{Int64}}) = rereference(dat, reference_channel[1])


dat = read_bdf("../Flank_C_3.bdf")
dat = eeg_data(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")




function rereference(dat::ContinuousData, reference_channel)
  dat_out = deepcopy(dat)
  rereference!(dat_out, reference_channel)
  return dat_out
end
rereference(dat::ContinuousData, reference_channel::Symbol) = rereference(dat, String(reference_channel))






# baseline functions
function baseline!(dat::ContinuousData)
  for col in names(dat.data)[3:end]
    dat.data[:, col] .-= mean(dat.data[:, col])
  end
end

function baseline!(dat::EpochData)
  for epoch in eachindex(dat.data)
    for col in names(dat.data[epoch])[3:end]
      dat.data[epoch][:, col] .-= mean(dat.data[epoch][:, col])
    end
  end
end

# baseline!(dat)
# plot_databrowser(dat)
# baseline!(epochs)
# plot_databrowser(epochs)


function plot_databrowser(dat::ContinuousData)

  fig = Figure()
  ax = GLMakie.Axis(fig[1, 1])  # plot layout

  xrange = GLMakie.Observable(1:4000) # default xrange
  yrange = GLMakie.Observable(-1500:1500) # default yrange

  # keyboard events
  on(events(fig).keyboardbutton) do event
    if event.action in (Keyboard.press, Keyboard.repeat)
      event.key == Keyboard.left && step_back(ax, xrange)
      event.key == Keyboard.right && step_forward(ax, nrow(dat.data), xrange)
      event.key == Keyboard.down && chans_less(ax, yrange)
      event.key == Keyboard.up && chans_more(ax, yrange)
    end
    return Consume(false)
  end

  xlims!(ax, dat.data.time[xrange.val[1]], dat.data.time[xrange.val[end]])
  ylims!(ax, yrange.val[1], yrange.val[end])
  ax.xlabel = "Time (ms)"
  ax.ylabel = "Amplitude (mV)"

  function draw(ax::Axis, xrange::Observable)
    for i = 3:(size(dat.data)[2]) # for all channels
      lines!(ax, @lift(dat.data[$xrange, 1]), @lift(dat.data[$xrange, i]))
    end
  end

  function step_back(ax::Axis, xrange::Observable)
    xrange.val[1] - 100 < 0 && return
    xrange[] = xrange.val .- 100
    xlims!(ax, dat.data.time[xrange.val[1]], dat.data.time[xrange.val[end]])
  end

  function step_forward(ax::Axis, xmax, xrange::Observable)
    xrange.val[1] + 100 > xmax && return
    xrange[] = xrange.val .+ 100
    xlims!(ax, dat.data.time[xrange.val[1]], dat.data.time[xrange.val[end]])
  end

  function chans_less(ax::Axis, yrange::Observable)
    (yrange.val[1] + 100 >= 0 || yrange.val[end] - 100 <= 0) && return
    yrange.val = yrange[][1]+100:yrange[][end]-100
    ylims!(ax, yrange.val[1], yrange.val[end])
  end

  function chans_more(ax::Axis, yrange::Observable)
    yrange.val = yrange[][1]-100:yrange[][end]+100
    ylims!(ax, yrange.val[1], yrange.val[end])
  end

  draw(ax, xrange)
  display(fig)

end

function plot_databrowser(dat::EpochData)

  fig = Figure()
  ax = GLMakie.Axis(fig[1, 1])  # plot layout

  xrange = GLMakie.Observable(1:nrow(dat.data[1])) # default xrange
  yrange = GLMakie.Observable(-1500:1500) # default yrange
  trial = GLMakie.Observable(1) # first trial

  # keyboard events
  on(events(fig).keyboardbutton) do event
    if event.action in (Keyboard.press, Keyboard.repeat)
      event.key == Keyboard.right && step_epoch_forward(trial)
      event.key == Keyboard.left && step_epoch_backward(trial)
      event.key == Keyboard.down && chans_less(ax, yrange)
      event.key == Keyboard.up && chans_more(ax, yrange)
    end
    return Consume(false)
  end

  xlims!(ax, dat.data[1].time[xrange.val[1]], dat.data[1].time[xrange.val[end]])
  ylims!(ax, yrange.val[1], yrange.val[end])
  vlines!(0, color=:gray, linewidth=2)
  ax.title = "Epoch $(trial.val)/$(length(dat.data))"
  ax.xlabel = "Time (ms)"
  ax.ylabel = "Amplitude (mV)"

  function draw(ax::Axis, xrange::Observable)
    for i = 2:(size(dat.data[1])[2]) # for all channels
      lines!(ax, @lift(dat.data[$trial][$xrange, 1]), @lift(dat.data[$trial][$xrange, i]))
    end
  end

  function step_epoch_forward(trial::Observable)
    trial[] = min(length(dat.data), trial.val[1] + 1)
    ax.title = "Epoch $(trial.val)/$(length(dat.data))"
  end

  function step_epoch_backward(trial::Observable)
    trial[] = max(1, trial.val[1] - 1)
    ax.title = "Epoch $(trial.val)/$(length(dat.data))"
  end

  function chans_less(ax::Axis, yrange::Observable)
    (yrange.val[1] + 100 >= 0 || yrange.val[end] - 100 <= 0) && return
    yrange.val = yrange[][1]+100:yrange[][end]-100
    ylims!(ax, yrange.val[1], yrange.val[end])
  end

  function chans_more(ax::Axis, yrange::Observable)
    yrange.val = yrange[][1]-100:yrange[][end]+100
    ylims!(ax, yrange.val[1], yrange.val[end])
  end

  draw(ax, xrange)
  display(fig)

end
plot_databrowser(epochs)



# basic plotting functions
function plot_epoch(dat::EpochData, trial::Int, channel; xlim=nothing, ylim=nothing, legend=true, xlabel="Time (S)", ylabel="mV")

  f = Figure()
  ax = GLMakie.Axis(f[1, 1])

  GLMakie.lines!(dat.trials[trial][!, :time], dat.trials[trial][!, channel], label=string(channel))

  !isnothing(xlim) && xlims!(ax, xlim)
  !isnothing(ylim) && xlims!(ax, ylim)
  legend && axislegend()
  ax.xlabel = xlabel
  ax.ylabel = ylabel

  return f

end


# basic plotting functions
function plot_epoch(dat::EpochData, trials, channels; xlim=nothing, ylim=nothing, legend=true, xlabel="Time (S)", ylabel="mV")
  f = Figure()
  ax = GLMakie.Axis(f[1, 1])

  for trial in trials
    for channel in channels
      GLMakie.lines!(dat.trials[trial][!, :time], dat.trials[trial][!, channel], label=string(channel))
    end
  end

  !isnothing(xlim) && xlims!(ax, xlim)
  !isnothing(ylim) && xlims!(ax, ylim)
  legend && axislegend()
  ax.xlabel = xlabel
  ax.ylabel = ylabel

  return f
end




function data_interpolation_topo(dat, points; radius=2.5, grid_scale=300)
  x = y = range(-radius, radius, length=grid_scale)
  X, Y = repeat(x', grid_scale)[:], repeat(y', grid_scale)[:]
  grid = [X Y]'
  dat = interpolate(Multiquadratic(), points, dat)
  dat = ScatteredInterpolation.evaluate(dat, grid)
  dat = reshape(dat, grid_scale, grid_scale)
  circle_mask!(dat, grid_scale)
  return dat
end

function head_shape(f, ax, layout; radius=2.5, linewidth=2)
  # head shape
  arc!(ax, Point2f(0), radius, -π, π, color=:black, linewidth=linewidth)
  nose = (Point2f[(-0.05, 0.5), (0.0, 0.55), (0.05, 0.5)] .* radius * 2)
  ear_right = (Point2f[
    (0.497, 0.0555), (0.51, 0.0775), (0.518, 0.0783),
    (0.5299, 0.0746), (0.5419, 0.0555), (0.54, -0.0055),
    (0.547, -0.0932), (0.532, -0.1313), (0.51, -0.1384),
    (0.489, -0.1199)] .* radius * 2)
  ear_left = ear_right .* Point2f(-1, 1)
  lines!(ax, nose, color=:black, linewidth=linewidth)
  lines!(ax, ear_right, color=:black, linewidth=linewidth)
  lines!(ax, ear_left, color=:black, linewidth=linewidth)

  # points/labels
  scatter!(ax, layout[!, :Dx], layout[!, :Dy], marker=:circle, markersize=10, color=:black)
  xoffset = 0
  yoffset = 0
  foreach(i -> text!(ax, position=(layout[!, :Dx][i] + xoffset, layout[!, :Dy][i] + yoffset), layout.label[i]), 1:nrow(layout))

  # hide some plot stuff
  hidexdecorations!(ax; label=true, ticklabels=true, ticks=true, grid=true, minorgrid=true, minorticks=true)
  hideydecorations!(ax; label=true, ticklabels=true, ticks=true, grid=true, minorgrid=true, minorticks=true)
  hidespines!(ax, :t, :r, :l, :b)

  return f

end

function head_shape(layout; radius=2.5, linewidth=2)
  f = Figure()
  ax = GLMakie.Axis(f[1, 1])
  head_shape(f, ax, layout, radius=radius, linewidth=linewidth)
end

head_shape(dat.layout, radius=2, linewidth=5)






function plot_topoplot(dat; ylim=nothing, radius=2.5, grid_scale=300)

  points = Matrix(dat.layout[!, [:X2, :Y2]])'
  data = data_interpolation_topo(Vector(dat.data[1000, 3:end]), points)

  if isnothing(ylim)
    ylim = minimum(data[.!isnan.(data)]), maximum(data[.!isnan.(data)])
  end

  f = Figure()
  ax = GLMakie.Axis(f[1, 1])
  co = contourf!(range(-radius, radius, length=grid_scale), range(-radius, radius, length=grid_scale), data,
    levels=100, colormap=:jet)
  Colorbar(f[1, 2], co)

  # head shape
  head_shape(f, ax, dat.layout)

  return f
end



plot_topoplot(dat, ylim=(-100, 100))

function circle_mask!(dat, grid_scale)
  for col in 1:size(dat)[1]
    for row in 1:size(dat)[2]
      xcentre = (grid_scale / 2) - col
      ycenter = (grid_scale / 2) - row
      if sqrt((xcentre^2 + ycenter^2)) > (grid_scale / 2)
        dat[col, row] = NaN
      end
    end
  end
end










# average reference
rereference!(eeg, 1:72);

# high-pass/low-pass filter
highpass_filter!(eeg, 0.1, 2);
lowpass_filter!(eeg, 30, 6);

# extract epochs and baseline
epoched_data = extract_epochs(eeg, 1, -0.5, 2);
baseline!(epoched_data, 0, 0)

# plot some epoched data
plot_epoch(epoched_data, 1, [:PO7, :PO8]) # single trial/multiple eleectrodes
plot_epoch(epoched_data, collect(1:10), [:PO7])   # multiple trials/single electrode

# average epochs
erp = average_epochs(epoched_data);

# plot some erp data
plot_erp(erp, [:Fp1, :Fp2])













