using BioSemiBDF
using CSV
using DSP
using DataFrames
using JLD2


mutable struct CoordXY
  label
  x
  y
end

mutable struct Coord
  coord::Array{CoordXY}
end

function polar2cartXY(inc, azi)
  inc = inc .* (pi / 180)
  azi = azi .* (pi / 180)
  return CoordXY(inc .* cos.(azi), inc .* sin.(azi))
end

function polar2cartXY(layout::DataFrame)
  inc = layout.inc .* (pi / 180)
  azi = layout.azi .* (pi / 180)
  return CoordXY(inc .* cos.(azi), inc .* sin.(azi))
end

mutable struct ContinuousData
  data::DataFrame
  layout::CoordXY
  sample_rate::Int
end

function create_dataframe(data::BioSemiBDF.BioSemiData)
  df1 = DataFrame(time=data.time, events=data.triggers.raw)
  df2 = DataFrame(data.data, :auto)
  rename!(df2, data.header.channel_labels[1:end-1])
  return hcat(df1, df2)
end

function eeg_data(dat::BioSemiBDF.BioSemiData, layout_file_name::String)
  return ContinuousData(create_dataframe(dat), polar2cartXY(DataFrame(CSV.File(layout_file_name))), dat.header.sample_rate[1])
end

function Base.show(io::IO, dat::ContinuousData)
  println("Timepoints x Channels: ", size(dat.data))
  println("Channel Labels: ", join(names(dat.data)[3:end], ", "))
  println("Sample Rate: ", dat.sample_rate)
end

mutable struct EpochData
  trials::Array{DataFrame}
  layout::CoordXY
  sample_rate::Int
end

function Base.show(io::IO, dat::EpochData)
  println("Number of trials: ", length(dat.trials))
  println("Timepoints x Channels: ", size(dat.trials[1]))
  println("Channel Labels: ", join(names(dat.trials[1])[3:end], ", "))
  println("Sample Rate: ", dat.sample_rate)
end

mutable struct ErpData
  avg::DataFrame
  layout::CoordXY
  sample_rate::Int
end

function Base.show(io::IO, dat::ErpData)
  println("Timepoints x Channels: $(nrow(dat.avg)) x $(ncol(dat.avg)-1)")
  println("Channel Labels: ", join(names(dat.avg)[2:end], ", "))
  println("Sample Rate: ", dat.sample_rate)
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
function search_sequence(array, sequence::Int)
  return findall(array .== sequence)
end

function search_sequence(array, sequence::Array{Int})

  idx_start_positions = findall(array .== sequence[1])

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

function extract_epochs(dat::ContinuousData, trigger_sequence, start_time, end_time; zero_position=1)
  idx_positions = search_sequence(dat.data.events, trigger_sequence) .+ (zero_position - 1)
  isempty(idx_positions) && error("No epochs matching triggerr sequence found!")
  n_pre = findmin(abs.(dat.data.time .- abs(start_time)))[2] - 1
  n_post = findmin(abs.(dat.data.time .- abs(end_time)))[2] - 1
  epochs = []
  for epoch in eachindex(idx_positions)
    df = DataFrame(dat.data[idx_positions[epoch]-n_pre:idx_positions[epoch]+n_post, :])
    df.time = df.time .- dat.data.time[idx_positions[epoch]]
    push!(epochs, df)
  end
  return EpochData(epochs, dat.layout, dat.sample_rate)
end


function average_epochs(dat::EpochData)
  vdf = reduce(vcat, dat.trials)
  erp = combine(groupby(vdf, :time), Not([:time, :events]) .=> mean .=> Not([:time, :events]))
  return ErpData(erp, dat.layout, dat.sample_rate)
end


dat = read_bdf("../Flank_C_3.bdf")
dat = eeg_data(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
epochs = extract_epochs(dat, 1, -0.5, 2)
erp = average_epochs(epochs)

using JLD2
save_object("test.jld2", dat)
dat = load_object("test.jld2")

using GLMakie

# basic plotting functions
function plot_epoch(dat::EpochData, trial::Int, channel; xlim=nothing)
  f = Figure()
  ax = GLMakie.Axis(f[1, 1])

  if isnothing(xlim)
    xlim = 1:nrow(dat.trials[trial])
  else
    xlim = findmin(abs.(dat.trials[trial].time .- abs(xlim[1])))[2]-1:findmin(abs.(dat.trials[trial].time .- abs(xlim[2])))[2]-1
  end

  GLMakie.lines!(dat.trials[trial][xlim, :time], dat.trials[trial][xlim, channel], label=string(channel))
  ax.xlabel = "Time (S)"
  ax.ylabel = "mV"
  axislegend()
  return f
end




# basic plotting functions
function plot_epoch(dat::EpochData, trials, channels; xlim=nothing)
  f = Figure()
  ax = GLMakie.Axis(f[1, 1])

  if isnothing(xlim)
    xlim = 1:nrow(dat.trials[trials[1]])
  else
    xlim = findmin(abs.(dat.trials[trials[1]].time .- abs(xlim[1])))[2]-1:findmin(abs.(dat.trials[trials[1]].time .- abs(xlim[2])))[2]-1
  end

  for trial in trials
    for channel in channels
      GLMakie.lines!(dat.trials[trial][xlim, :time], dat.trials[trial][xlim, channel], label=channel)
    end
  end
  ax.xlabel = "Time (S)"
  ax.ylabel = "mV"

  return f
end

plot_epoch(epochs, 1:20, [:Fp2])

# re-reference functions
function rereference!(dat::ContinuousData, reference_channel)
  reference = reduce(+, eachcol(dat.data[:, reference_channel])) ./ length(reference_channel)
  for col in names(dat.data)[3:end]
    dat.data[:, col] .= dat.data[:, col] .- reference
  end
end

@btime rereference!(dat, [:Fp1, :Fp2])


# baseline functions
function baseline!(dat::ContinuousData)
  for col in names(dat.data)[3:end]
    dat.data[:, col] .= dat.data[:, col] .- mean(dat.data[:, col])
  end
end

@btime baseline!(dat)



# Filter functions
function lowpass_filter!(dat::ContinuousData, freq, order)
  filter = digitalfilter(Lowpass(freq, fs=dat.sample_rate), Butterworth(order))
  for col in names(dat.data)[3:end]
    dat.data[:, col] .= filtfilt(filter, dat.data[:, col])
  end
end



# Filter functions
function lowpass_filter!(dat::ContinuousData, freq, order)
  filter = digitalfilter(Lowpass(freq, fs=dat.sample_rate), Butterworth(order))
  for col in names(dat.data)[3:end]
    dat.data[:, col] .= filtfilt(filter, dat.data[:, col])
  end
end

# Filter functions
function highpass_filter!(dat::ContinuousData, freq, order)
  filter = digitalfilter(Highpass(freq, fs=dat.sample_rate), Butterworth(order))
  for col in names(dat.data)[3:end]
    dat.data[:, col] .= filtfilt(filter, dat.data[:, col])
  end
end

dat = read_bdf("../Flank_C_3.bdf")
dat = eeg_data(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
lines(dat.data.Fp1)
highpass_filter!(dat, 0.5, 2)
lowpass_filter!(dat, 5, 6)
lines!(dat.data.Fp1)


resample(dat.data.time, 0.5)


function downsample!(dat::ContinuousData, dec::Int)

  for col in names(dat.data)[3:end]
    dat.data[:, col] .= resample(dat.data[:, col], 1 / dec)
  end


  Fp1 =
    t = range(dat.data.time[1], dat.data.time[end], length(Fp1))
  dat.sample_rate = div(dat.sample_rate, dec)

  # update triggers
  dat.events.raw = zeros(Int16, 1, size(dat.data, 1))
  dat.events.idx = convert(Array{Int64}, round.(events.idx / dec))
  dat.events.raw[events.idx] = events.val

  return dat

end

dat.data.time

dat.data.time[1]::dat.data.time[end]








function shrink_coordinates(x, y, limit)
  while sum(sqrt.((x .^ 2) + (y .^ 2)) .> limit) != 0
    x *= 0.99
    y *= 0.99
  end
  return (x, y)
end

function circle_shape(radius=1, xpos=0, ypos=0)
  pos = range(0, stop=2 * pi, length=180)
  return Shape((cos.(pos) .* radius) .+ xpos, (sin.(pos) * radius) .+ ypos)
end

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

function plot_topography_v1(dat, start_time, end_time; clim=(-30, 30), grid_scale::Int=300)
  coords = shrink_coordinates(dat.layout.x, dat.layout.y, 1.05)
  samples = Statistics.mean(dat.avg[(start_time)Unitful.s..(end_time)Unitful.s, :].avg, dims=1)'
  points = hcat(coords[1], coords[2])'
  itp = interpolate(Multiquadratic(), points, samples)
  x = y = range(-1, 1, length=grid_scale)
  X = repeat(x', grid_scale)[:]
  Y = repeat(y, grid_scale)[:]
  gridPoints = [X Y]'
  interpolated = ScatteredInterpolation.evaluate(itp, gridPoints)
  gridded = reshape(interpolated, grid_scale, grid_scale)
  circle_mask!(gridded, grid_scale)
  p = contourf(x, y, levels=50, gridded, legend=:false, border=:none, colorbar=:true, c=:jet, clims=(-30, 30), linewidth=0)
  # p = heatmap(x, y, levels = 50, gridded, legend = :false, border = :none, colorbar = :true, linewidth=:none, c=:jet, clims = clim)
  plot!(circle_shape(1), fillcolour=:white, linecolor=:black, lw=2, aspect_ratio=1, fillalpha=0)
  scatter!(coords[1], coords[2], markercolor=:black, series_annotations=text.(AxisArrays.axes(erp.data, Axis{2})[:], 8, "monaco", :top))
  return p
end

plot_topography_v1(erp, 0.2, 0.3, clim=(-10, 10))

function plot_topography_v2(dat, start_time, end_time; clim=(-30, 30), grid_scale::Int=300)
  coords = shrink_coordinates(dat.layout.x, dat.layout.y, 1.05)
  samples = Statistics.mean(dat.avg[(start_time)Unitful.S..(end_time)Unitful.S, :], dims=1)[:]
  spl = Spline2D(coords[2], coords[1], samples; kx=5, ky=5, s=50)
  xs2 = ys2 = range(-1, 1, length=grid_scale)
  data = evalgrid(spl, xs2, ys2)
  circle_mask!(data, grid_scale)
  p = contourf(xs2, ys2, levels=50, data, legend=:false, border=:none, colorbar=:true, c=:jet, clims=(-30, 30))
  # p = heatmap(x, y, levels = 50, gridded, legend = :false, border = :none, colorbar = :true, linewidth=:none, c=:jet, clims = clim)
  plot!(circle_shape(1), fillcolour=:white, linecolor=:black, lw=2, c=:jet, aspect_ratio=1, fillalpha=0)
  scatter!(coords[1], coords[2], markercolor=:black, series_annotations=text.(AxisArrays.axes(erp.avg, Axis{2})[:], 8, "monaco", :top))
  return p
end

plot_topography_v2(erp, 0, 1, clim=(-10, 10))

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

# Used biosemipy PyQt dataviewer via PyCall
function plot_bdf()
  try
    dataviewer = pyimport("biosemipy.dataviewer")
    dataviewer.run()
  catch
    nothing
  end
end

function plot_bdf(filename::String)
  try
    dataviewer = pyimport("biosemipy.dataviewer")
    dataviewer.run(filename)
  catch
    nothing
  end
end

# using biosemipy
plot_bdf()
plot_bdf("/home/ian/Documents/Julia/Flank_C_3.bdf")

function find_nearest_time_index(time::Array{Float64}, start_point, end_point)
  idx1 = findmin(abs.(time .- start_point))[2]
  idx2 = findmin(abs.(time .- end_point))[2]
  return [idx1, idx2]
end

function find_nearest_time_index(time::StepRangeLen, start_point, end_point)
  time = ustrip(time)
  idx1 = findmin(abs.(time .- start_point))[2]
  idx2 = findmin(abs.(time .- end_point))[2]
  return [idx1, idx2]
end









function downsample(dat::AxisArray, dec::Int)
  time = ustrip(get_axis(dat, 1))
  chan = get_axis(dat, 2)
  dat = mapslices(x -> resample(x, 1 // dec), dat.data; dims=1)
  time = time_range(time[1], time[end], size(dat, 1))

  return AxisArray(dat, time=time, chan=chan)
end

function downsample(dat_in::ContinuousData, dec::Int)

  data_out = deepcopy(dat_in)
  data_out.data = downsample(dat_in.data, dec)
  data_out.sample_rate = div(dat_in.sample_rate, dec)

  # update triggers
  data_out.events = deepcopy(dat_in.events)

  data_out.events.raw = zeros(Int16, 1, size(data_out.data, 1))
  data_out.events.idx = convert(Array{Int64}, round.(data_out.events.idx / dec))
  data_out.events.raw[data_out.events.idx] = data_out.events.val

  return data_out

end

function downsample!(dat::ContinuousData, dec::Int)

  dat.data = downsample(dat.data, dec)
  dat.sample_rate = div(dat.sample_rate, dec)

  # update triggers
  dat.events.raw = zeros(Int16, 1, size(dat.data, 1))
  dat.events.idx = convert(Array{Int64}, round.(events.idx / dec))
  dat.events.raw[events.idx] = events.val

  return dat

end

eeg1 = downsample(eeg, 2)
eeg2 = downsample(eeg1, 2)
eeg3 = downsample(eeg2, 2)

downsample(eeg, 2)
downsample(eeg1, 2)
downsample(eeg2, 2)

function downsample(dat_in::EpochData, dec::Int)
  data_out = deepcopy(dat_in)
  for trl = 1:length(dat_in.data)
    data_out.data[trl] = downsample(data_out.data[trl], dec)
  end
  data_out.sample_rate = div(dat_in.sample_rate, dec)
  return data_out
end

function downsample!(dat::EpochData, dec::Int)
  for trl = 1:length(dat.data)
    dat.data[trl] = downsample(dat.data[trl], dat.sample_rate, dec)
  end
  dat.sample_rate = div(dat.sample_rate, dec)
  return dat
end

function downsample(dat::ErpData, dec::Int)
  dat = deepcopy(dat)
  dat.data = downsample(dat.data, dec)
  dat.sample_rate = div(dat.sample_rate, dec)
  return dat
end

function downsample!(dat::ErpData, dec::Int)
  dat.data = downsample(dat.data, dec)
  dat.sample_rate = div(dat.sample_rate, dec)
  return dat
end

@time downsample!(erp, 2)

erp.data

typeof(erp)

using GLMakie
x = range(0, 10, length=100)
y = sin.(x)
scatter(x, y)

data = cumsum(randn(4, 101), dims=2)

dat.data'

using GLMakie

fig, ax, sp = series(dat.time, dat.data[:, 1:7]', linewidth=4, labels=["label $i" for i in 1:72])
xlims!(0, 8000)
ylims!(-20000, 20000)
axislegend(ax)
#fig

f = Figure()
ax = GLMakie.Axis(f[1, 1])
for y in 1:100
  for i in 1:72
    # lines!(dat.time[1:2000], dat.data[1:2000, i], linewidth=3, )
    lines!(dat.time[y:2000+y], dat.data[y:2000+y, i], linewidth=3)
  end
  sleep(0.1)
  empty!(ax)
end


xlims!(0, 8)



using GLMakie

f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98), resolution=(1500, 700))
ax = GLMakie.Axis(f[1, 1], xlabel="Time [s]", ylabel="Voltage amplitude [µV]")

#N = 1:length(pos)
positions = Observable(rand(Point2f, 10))

xs = 0:0.01:10
ys = 0.5 .* sin.(xs)

lines!(xs, ys)
lines!(xs, ys * 2)

hidedecorations!(ax, label=false, ticks=false, ticklabels=false)
hidespines!(ax, :t, :r)
hlines!(0, color=:gray, linewidth=1)
vlines!(0, color=:gray, linewidth=1)

register_interaction!(ax, :my_interaction) do event, axis
  if event.type === MouseEventTypes.leftclick
    println("Graph axis position: $(event.data)")
  end
end

i = Observable(0)
on(events(f).mousebutton, priority=2) do event
  if event.button == Mouse.left && event.action == Mouse.press
    plt, i[] = pick(f)
    str = lift(i -> "$(i)", i)
    text!(ax, 1, -0.5, text=str, align=(:center, :center))
    @show mouseposition(f)
  end
end
f





using CairoMakie


f = Figure()
CairoMakie.Axis(f[1, 1])
xs = 0:0.01:10
ys = 0.5 .* sin.(xs)
for (i, lw) in enumerate([1, 2, 3])
  lines!(xs, ys .- i / 6, linestyle=nothing, linewidth=lw)
  lines!(xs, ys .- i / 6 .- 1, linestyle=:dash, linewidth=lw)
  lines!(xs, ys .- i / 6 .- 2, linestyle=:dot, linewidth=lw)
  lines!(xs, ys .- i / 6 .- 3, linestyle=:dashdot, linewidth=lw)
  lines!(xs, ys .- i / 6 .- 4, linestyle=:dashdotdot, linewidth=lw)
  lines!(xs, ys .- i / 6 .- 5, linestyle=Linestyle([0.5, 1.0, 1.5, 2.5]), linewidth=lw)
end
f


using GLMakie, MousetrapMakie
canvas = GLMakieArea()
window = Mousetrap.Window()
set_child!(window, canvas) # can be used like any other widget

screen = create_glmakie_screen(canvas)
display(screen, scatter(rand(123)))





@time a = downsample(epoched_data, 2)


@time downsample!(epoched_data, 2)


dat = read_bdf("../Flank_C_3.bdf")

f = Figure()
GLMakie.Axis(f[1, 1])
for i in 1:71
  lines!(dat.data[:, i], linestyle=nothing, linewidth=2)
end
xlims!(0, 2000)
ylims!(-200, 200)

lines!(dat.data[1:1000, 1], linestyle=nothing, linewidth=2)
f



f = Figure()
Axis(f[1, 1])
xs = 0:0.01:10
ys = 0.5 .* sin.(xs)
for (i, lw) in enumerate([1, 2, 3])
  lines!(xs, ys .- i / 6, linestyle=nothing, linewidth=lw)
  lines!(xs, ys .- i / 6 .- 1, linestyle=:dash, linewidth=lw)
  lines!(xs, ys .- i / 6 .- 2, linestyle=:dot, linewidth=lw)
  lines!(xs, ys .- i / 6 .- 3, linestyle=:dashdot, linewidth=lw)
  lines!(xs, ys .- i / 6 .- 4, linestyle=:dashdotdot, linewidth=lw)
  lines!(xs, ys .- i / 6 .- 5, linestyle=[0.5, 1.0, 1.5, 2.5], linewidth=lw)
end
f


using GLMakie
xs = 1:100
ys = Node(rand(length(xs)) .- 0.5);
plot(xs, ys)
for _ in 1:30
  ys[] = ys[] .+ 0.1 * rand(length(xs)) .- 0.05
  sleep(0.1)
end


#################### ICA ##################
using Random
using MAT
using Statistics
using LinearAlgebra
using Printf
using BenchmarkTools


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


# data = read_mat_file("data.mat")
# data["data"]["trial"][1]

# # baseline
# for i = 1:length(data["data"]["trial"])
#     data["data"]["trial"][i] = data["data"]["trial"][i] .- mean(data["data"]["trial"][i], dims = 2)
# end

# # scale
# scale = sqrt(opnorm((data["data"]["trial"][1]*transpose(data["data"]["trial"][1])) ./ size(data["data"]["trial"][1], 2)))
# for i = 1:length(data["data"]["trial"])
#     data["data"]["trial"][i] = data["data"]["trial"][i] ./ scale
# end




function pca_reduction!(dat, n)
  PCdat2 = dat'           # transpose data
  PCn = size(PCdat2, 1) # now p chans, n time points
  PCdat2 = PCdat2 / PCn
  PCout = dat * PCdat2
  PCD, PCV = eigen(PCout, sortby=x -> -abs(x))
  dat = PCV[:, 1:ncomps]' * dat
  return dat, PCV
end


function infoica(dat, ncomps; maxsteps=512, lrate=0.001, minlrate=0.000001, sphering=true)

  # define some default values
  MAX_WEIGHT = 1e8          # guess that weights larger than this have blown up
  DEFAULT_RESTART_FAC = 0.9          # if weights blowup, restart with lrate
  degconst = 180 ./ pi

  nchans, datalength = size(dat)

  DEFAULT_BLOWUP = 1000000000.0 # learning rate has 'blown up'
  DEFAULT_BLOWUP_FAC = 0.8          # when lrate 'blows up,' anneal by this fac
  ANNEALDEG = 60.0
  ANNEALSTEP = 0.90

  wts_blowup = false
  pcaflag = true
  sphering = true
  blocksize = Int(ceil(min(5 * log(datalength), 0.3 * datalength))) # heuristic default - may need adjustment!
  minchange = 0.000001

  # remove overall row means of data
  dat .-= mean(dat, dims=2)

  # Perform PCA reduction
  if ncomps < nchans
    dat, PCV = pca_reduction!(dat, ncomps)
  end

  # Perform sphering
  if sphering
    @printf("Computing the sphering matrix and sphering data...\n")
    sphere = 2.0 * inv(sqrt(cov(dat')))
    dat = sphere * dat
  end

  # initial weights
  weights = Matrix{Float64}(I, ncomps, ncomps)

  # initialize ICA training

  lastt = Int(floor((datalength / blocksize - 1) * blocksize + 1))

  BI = Matrix{Float64}(I, ncomps, ncomps) * blocksize
  delta = zeros(1, ncomps * ncomps)
  startweights = weights
  oldweights = weights
  prevweights = weights
  onesrow = ones(1, Int(blocksize))
  bias = zeros(ncomps, 1)

  olddelta = 1.0
  oldchange = 0.0

  u = Array{Float64,2}(undef, ncomps, blocksize)
  y = Array{Float64,2}(undef, ncomps, blocksize)
  yu = Array{Float64,2}(undef, ncomps, ncomps)
  qq = Array{Float64,2}(undef, ncomps, ncomps)

  timeperm = randperm(datalength)
  timepermv = Array{Int}(undef, blocksize)
  datv = Array{Float64,2}(undef, ncomps, blocksize)
  wts_blowup = false

  nstep = 0
  while nstep < maxsteps

    for t = 1:blocksize:lastt

      timepermv .= @view(timeperm[t:t+blocksize-1])
      datv .= @view(dat[:, @view(timeperm[t:t+blocksize-1])])

      mul!(u, weights, datv) .+= bias .* onesrow
      y .= 1.0 .- ((1.0 ./ (1.0 .+ exp.(-u))) .* 2.0)
      weights += mul!(qq, ((BI .+ mul!(yu, y, u'))), weights) .* lrate
      bias += sum(y, dims=2) .* lrate'

      if any(abs.(weights) .> MAX_WEIGHT)
        wts_blowup .= true
        break
      end

    end

    if !wts_blowup
      oldwtchange = weights .- oldweights
      nstep += 1
      # compute angle changes
      angledelta = 0.0
      delta = reshape(oldwtchange, 1, ncomps * ncomps)
      change = (delta*delta')[1]

      if nstep > 2
        angledelta = ((acos((delta * olddelta') / sqrt(change .* oldchange)))[1]) * degconst
      end

      @printf "step: %d, lrate: %.10f, wchange: %.10f, angledelta %.2f\n" nstep lrate change angledelta
      #
      #%%%%%%%%%%%%%%%%%%% save current values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      oldweights .= weights

      #%%%%%%%%%%%%%%%%%%% anneal learning rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if angledelta > ANNEALDEG
        lrate = lrate * ANNEALSTEP     # anneal learning rate
        olddelta = delta                # accumulate angledelta until
        oldchange = change               #  annealdeg is reached
      elseif nstep === 1                    # on first step only
        olddelta = delta                # initialize
        oldchange = change
      end

      #%%%%%%%%%%%%%%%%%%% apply stopping rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if nstep > 2 && change < minchange    #% apply stopping rule
        nstep = maxsteps                   # stop when weights stabilize
      elseif change > DEFAULT_BLOWUP       # if weights blow up,
        lrate = lrate * DEFAULT_BLOWUP_FAC  # keep trying
      end                                  # with a smaller learning rate

    else

      # start again/re-initialize variables with lower learing rate
      step = 0
      wts_blowup = false
      lrate = lrate * DEFAULT_RESTART_FAC

      weights = Matrix{Float64}(I, ncomps, ncomps)
      oldweights = startweights
      # delta        = zeros(1, ncomps*ncomps);
      olddelta = delta
      # prevweights  = startweights;
      bias = zeros(ncomps, 1)

      if lrate > minlrate
        r = rank(dat) # determine if data rank is too low
        if rand(dat) < ncomps
          @printf "Data has rank %d. Cannot compute %d components.\n" r ncomps
          return
        else
          @printf "Lowering learning rate to %f and starting again\n" lrate
        end
      else
        @printf "QUITTING - weight matrix may not be invertible!\n"
        return
      end
    end

  end

  if pcaflag
    weights = weights * sphere * transpose(PCV[:, 1:ncomps])
    sphere = Diagonal(ones(nchans, nchans))
  end

  # %%%%%%%%%%%%%%%%%%%% Find mean variances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ncomps == nchans # if weights are square . . .
    winv = inv(weights * sphere)
  else
    # fprintf('Info (runica): Using pseudo-inverse of weight matrix to rank order component projections.\n');
    winv = pinv(weights * sphere)
  end

  # println(size(weights))
  # compute variances without backprojecting to save time and memory -sm 7/05
  meanvar = sum(winv .^ 2, dims=1) .* sum((dat') .^ 2, dims=1) / ((ncomps * datalength) - 1)

  windex = sortperm(meanvar[1, :])
  windex = windex[ncomps:-1:1] # order large to small
  weights = weights[windex, :]   # final weights

  return weights, sphere

end

dat = read_mat_file("dat.mat")["dat"]
@btime weight, sphere = infoica(dat, 68)
weight, sphere = infoica(dat)



using Plots

plot(dat[1, :])
for i in 2:60
  plot!(dat[i, :])
end



using GLMakie
points = [Point2f0(cos(t), sin(t)) for t in LinRange(0, 2pi, 20)]
colors = 1:20
figure, axis, scatterobject = scatter(points, color=colors, markersize=20)
figure

fig = Figure(resolution=(1200, 900))
lines(dat[1, :], axis=(frame=(linewidth=0,),))
fig

for i in 1:10
  lines(fig[i, 1], dat[i, :], axis=(frame=(linewidth=0,),))
end
fig
using AbstractPlotting.MakieLayout
using AbstractPlotting
scene, layout = layoutscene(resolution=(1200, 900))

ax1 = layout[1, 1] = LAxis(scene, title="Axis 1")
ax2 = layout[1, 2] = LAxis(scene, title="Axis 2")

hidespines!(ax1)
hidespines!(ax2, :t, :r) # only top and right


for (i, n) in enumerate([2, 5, 9])
  lines(fig[i, 1], 0 .. 20, sin, axis=(frame=(linewidth=0,),))
end
fig

using Makie

fs = [sin, cos, sinc]
colors = [:red, :green, :blue]

xspan = 0:0.01:2*pi
y = [fs[i](x) for x in xspan, i in 1:size(fs, 1)]

scene = Scene()

for i in 1:size(fs, 1)
  lines!(scene, xspan, y[:, i], color=colors[i])
end

scene

lines(dat[1, :], dat[2, :])


using GLMakie, Colors
scene = lines([0.25, 0.75], [0.25, 0.75], color=RGBA(0, 0, 0, 0.5), linewidth=30)
lines!([0.75, 0.25], [0.25, 0.75], color=RGBA(0, 0, 0, 0.5), linewidth=10)
scene # or save scene...
## This shows the same thing...
scene = lines([0.25, 0.75, NaN, 0.75, 0.25], [0.25, 0.75, NaN, 0.25, 0.75], color=RGBA(0, 0, 0, 0.5), linewidth=3)
scene #



scene = Scene()
for i in 1:3
  lines!(scene, 1:size(dat)[2], dat[i, :] .+ i .* 1000, color=colors[i])
end
scene

function example_plot()
  f = Figure(resolution=(1000, 800))
  for i in 1:2, j in 1:2
    lines(f[i, j], cumsum(randn(1000)))
  end
  f
end
example_plot()


unmixing = weights * sphere

mixing = [];

if (isempty(mixing) && !isempty(unmixing))
  if (size(unmixing, 1) == size(unmixing, 2)) && rank(unmixing) == size(unmixing, 1)
  else
    mixing = pinv(unmixing)
  end
end

unmixing
mixing

labels = []
for i = 1:size(mixing, 2)
  push!(labels, @sprintf("infoica_%d", i))
end


plot(dat[1, 1:2000])


out = InfoICA(mixing, unmixing, labels)
out.topo
out.unmixing

tra = Diagonal(ones(69, 69)) - (out.topo[:, 1] * transpose(out.unmixing[1, :]))

dat1 = tra * dat

plot(dat[1, 642:1400])
plot!(dat1[1, 642:1400])
