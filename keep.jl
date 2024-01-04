mutable struct CoordXY
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

using JLD2
save_object("test.jld2", dat)
dat = load_object("test.jld2")


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

function downsample(dat::AxisArray, dec::Int)
  time = ustrip(get_axis(dat, 1))
  chan = get_axis(dat, 2)
  dat = mapslices(x -> resample(x, 1 // dec), dat.data; dims=1)
  time = time_range(time[1], time[end], size(dat, 1))

  return AxisArray(dat, time=time, chan=chan)
end




# hidedecorations!(ax, label = false, ticks = false, ticklabels = false) 
# hidespines!(ax, :t, :r) 
# hlines!(0, color = :gray, linewidth = 1)
# vlines!(0, color = :gray, linewidth = 1)




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
#

display(GLMakie.Screen(), lines(epochs.data[1].Fp1))
display(GLMakie.Screen(), lines(epochs1.data[1].Fp1))


mapcols(x -> filtfilt(filter, x), dat.data[:, 3:end])
mapcols!(x -> filtfilt(filter, x), dat.data[:, 3:end])


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


