using BioSemiBDF
using CSV
using DSP
using DataFrames
using GLMakie
using JLD2
using ScatteredInterpolation
using StatsBase
using MAT


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

function head_shape(f, ax, layout; radius=2.5, linewidth=2, plot_points=true, plot_labels=true, label_x_offset=0, label_y_offset=0)
  # head shape
  arc!(ax, Point2f(0), radius, -π, π, color=:black, linewidth=linewidth) # head
  arc!(Point2f(radius, 0), radius / 7, -π / 2, π / 2, color=:black, linewidth=linewidth) # ear right
  arc!(Point2f(-radius, 0), -radius / 7, π / 2, -π / 2, color=:black, linewidth=linewidth) # ear left
  lines!(ax, Point2f[(-0.05, 0.5), (0.0, 0.6), (0.05, 0.5)] .* radius * 2, color=:black, linewidth=linewidth) # nose

  # points
  if plot_points
    scatter!(ax, layout[!, :X2], layout[!, :Y2], marker=:circle, markersize=10, color=:black)
  end

  if plot_labels
    foreach(i -> text!(ax, position=(layout[!, :X2][i] + label_x_offset, layout[!, :Y2][i] + label_y_offset), layout.label[i]), 1:nrow(layout))
  end

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


function plot_topoplot(dat; ylim=nothing, radius=2.5, grid_scale=300, plot_points=true, plot_labels=true, label_x_offset=0, label_y_offset=0)

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
  head_shape(f, ax, dat.layout, plot_points=plot_points, plot_labels=plot_labels, label_x_offset=label_x_offset, label_y_offset=label_y_offset)

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



function read_mat_file(filename)
  file = matopen(filename)
  dat = read(file)
  close(file)
  return dat
end


struct InfoICA
  topo
  unmixing
  label
end




function pca_reduction!(dat, n)
  PCdat2 = dat'           # transpose data
  PCn = size(PCdat2, 1) # now p chans, n time points
  PCdat2 = PCdat2 / PCn
  PCout = dat * PCdat2
  PCD, PCV = eigen(PCout, sortby=x -> -abs(x))
  dat = PCV[:, 1:ncomps]' * dat
  return dat, PCV
end


function infoica(dat; extended=true)

  # define some default values
  l_rate = 0.001 # initial learning rate
  max_iter = 200 # maximum number of iterations
  w_change = 1e-12 # change to stop iteration
  use_bias = true
  anneal_deg = 60.0 # angle at which learning rate reduced
  anneal_step = 0.9 # factor by which learning rate reduced
  blowup = 1e4 # max difference allowed between two successive estimations of unmixing matrix
  blowup_fac = 0.5 # factor by which learning rate will be reduced if "blowup"
  n_small_angle = 20
  max_weight = 1e8     # larger than this have "blowup"
  restart_factor = 0.9 # if weights blowup, restart with lrate
  degconst = 180 ./ pi

  # for extended Infomax
  ext_blocks = 1
  n_subgauss = 1
  extmomentum = 0.5
  signsbias = 0.02
  signcount_threshold = 25
  signcount_step = 2
  kurt_size = 6000

  # check data shape and get block sizes
  n_features, n_samples = size(dat)
  n_features_square = n_features^2
  block = trunc(Int, sqrt(n_samples / 3.0))
  nblock = div(n_samples, block)
  lastt = (nblock - 1) * block + 1

  # start ICA
  @printf "Computing Infomax ICA\n"

  # initialize training
  weights = Matrix{Float64}(I, n_features, n_features) # identity matrix
  BI = block * Matrix{Float64}(I, n_features, n_features)
  bias = zeros(n_features, 1)
  onesrow = ones((1, block))
  startweights = copy(weights)
  oldweights = copy(startweights)
  step = 0
  count_small_angle = 0
  wts_blowup = false
  blockno = 0
  signcount = 0
  initial_ext_blocks = ext_blocks  # save the initial value in case of reset

  if extended
    signs = ones(n_features)
    for k in 1:n_subgauss
      signs[k] = -1
    end
    kurt_size = min(kurt_size, n_samples)
    old_kurt = zeros(n_features)
    oldsigns = zeros(n_features)
  end

  # trainings loop
  olddelta = 1.0
  oldchange = 0.0
  while step < max_iter

    permute = shuffle(1:n_samples)

    for t = 1:block:lastt

      u = *(dat[:, permute[t:t+block-1]]', weights) + *(bias, onesrow)'

      if extended
        # extended ICA update
        y = tanh.(u)
        weights += l_rate * *(weights, BI .- signs[:]' * *(u', y) - *(u', u))
        if use_bias
          bias += (l_rate .* sum(y, dims=1) .* -2.0)'
        end
      else
        y = 1 ./ (1 .+ exp.(-u))
        weights += l_rate * *(weights, BI + *(u', (1.0 .- 2.0 .* y)))
        if use_bias
          bias += (l_rate .* sum((1.0 .- 2.0 .* y), dims=1))'
        end
      end

      # check change limit
      if maximum(abs.(weights)) > max_weight
        wts_blowup = true
        break
      end

      blockno += 1
      # ICA kurtosis estimation
      if extended
        if ext_blocks > 0 & blockno % ext_blocks == 0
          if kurt_size < n_samples
            rp = trunc.(Int, (rand(kurt_size) .* (n_samples - 1))) .+ 1
            tpartact = *(dat[:, rp]', weights)'
          else
            tpartact = *(dat', weights)'
          end
          # estimate kurtosis
          kurt = kurtosis.(eachrow(tpartact))
          if extmomentum != 0
            kurt = extmomentum .* old_kurt .+ (1.0 .- extmomentum) .* kurt
            old_kurt = kurt
          end
          # estimate weighted signs
          signs = sign.(kurt .+ signsbias)
          ndiff = sum(signs .- oldsigns .!= 0)
          ndiff == 0 ? signcount += 1 : signcount = 0
          oldsigns = signs
          if signcount >= signcount_threshold
            ext_blocks = trunc.(ext_blocks .* signcount_step)
            signcount = 0
          end
        end
      end

    end

    if !wts_blowup
      oldwtchange = weights .- oldweights
      step += 1
      angledelta = 0.0
      delta = oldwtchange[:]
      change = sum(delta .* delta)
      if step > 2
        angledelta = acos(sum(delta .* olddelta) / sqrt(change * oldchange))
        angledelta *= degconst
      end

      @printf "step %d: lrate %5f, wchange %8.8f, angledelta %4.1f\n" step l_rate change angledelta

      # anneal learning rate
      oldweights = copy(weights)
      if angledelta > anneal_deg
        l_rate *= anneal_step  # anneal learning rate
        # accumulate angledelta until anneal_deg reaches l_rate
        olddelta = delta
        oldchange = change
        count_small_angle = 0  # reset count when angledelta is large
      else
        if step == 1  # on first step only
          olddelta = delta  # initialize
          oldchange = change
        end
        if !isnothing(n_small_angle)
          count_small_angle += 1
          if count_small_angle > n_small_angle
            max_iter = step
          end
        end
      end

      # apply stopping rule
      if step > 2 && change < w_change
        step = max_iter
      elseif change > blowup
        l_rate *= blowup_fac
      end

      # restart if weights blow up (for lowering l_rate)
    else
      step = 0  # start again
      wts_blowup = false  # re-initialize variables
      blockno = 1
      l_rate *= restart_factor  # with lower learning rate
      weights = copy(startweights)
      oldweights = copy(startweights)
      olddelta = zeros((1, n_features_square))
      bias = zeros((n_features, 1))

      ext_blocks = initial_ext_blocks

      # for extended Infomax
      if extended
        signs = ones(n_features)
        for k in 1:n_subgauss
          signs[k] = -1
        end
        oldsigns = zeros(n_features)
      end

      @printf "lowering learning rate to %g ... re-starting ...\n" l_rate

    end
  end

  return weights

end

dat = read_mat_file("../dat.mat")["dat"]
@time weights = infoica(dat)












