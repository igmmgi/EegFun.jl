using BioSemiBDF
using CSV
using DSP
using DataFrames
using GLMakie
using JLD2
using LinearAlgebra
using MAT
using OrderedCollections
using Printf
using Random
using ScatteredInterpolation
using StatsBase

include("types.jl")
include("utils.jl")
include("analyse.jl")
include("baseline.jl")
include("channel_difference.jl")
include("epochs.jl")
include("filter.jl")
include("layout.jl")
include("plot.jl")
include("rereference.jl")
include("topo.jl")
include("utils.jl")


# basic layouts
layout = read_layout("/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv");
# head_shape_2d(layout)
# head_shape_3d(layout)

# read bdf file
subject = 3
dat = read_bdf("../Flank_C_$(subject).bdf")
dat = create_eeg_dataframe(dat, layout)
# basic bdf plot
# plot_databrowser(dat)
filter_data!(dat, "hp", 0.1, 2)
# filter_data!(dat, "lp", 30, 6)
include("plot.jl")
# calculate EOG channels
diff_channel!(dat, ["Fp1", "Fp2"], ["IO1", "IO2"], "vEOG");
diff_channel!(dat, "F9", "F10", "hEOG");
# autodetect EOG signals
detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)
dat.data[!, "is_extreme"] .= is_extreme_value(dat.data, dat.layout.label, 100);
plot_databrowser(dat)
# plot_databrowser(dat, [dat.layout.label; "hEOG"; "vEOG"])
# plot_databrowser(dat, ["vEOG", "hEOG"])
#plot_databrowser(dat, ["hEOG"])





# extract epochs
epochs = extract_epochs(dat, 1, -0.5, 2)

# plot epochs
plot_databrowser(epochs)
plot_epochs(epochs, ["PO7", "PO8"])
plot_epochs(epochs, [:PO7])

# average epochs
erp = average_epochs(epochs)




plot_databrowser(erp)


save_object("$(subject)_$(cond)_epochs.jld2", epochs)
save_object("$(subject)_$(cond)_erp.jld2", erp)



diff_channel!(dat, "F9", "F10", "hEOG");
diff_channel!(dat, ["Fp1", "Fp2"], ["IO1", "IO2"], "vEOG");



detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)

function test_plot_eog_detection(dat, xlim, channel, detected)
  fig = Figure()
  ax = GLMakie.Axis(fig[1, 1])  # plot layout
  lines!(ax, dat.data.time[xlim], dat.data[!, channel][xlim])
  vlines!(ax, dat.data.time[xlim][dat.data[!, detected][xlim]], color=:black)
  display(fig)
end
test_plot_eog_detection(dat, 1000:14000, "vEOG", "is_vEOG")
test_plot_eog_detection(dat, 1000:4000, "hEOG", "is_hEOG")



function test_analysis()
  for subject in 3:4
    for condition in [1, 3]
      println("Reading file: Flank_C_$(subject).bdf")
      dat = read_bdf("../test_data/Flank_C_$(subject).bdf")
      dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
      filter_data!(dat, "hp", 1, 2)
      epochs = extract_epochs(dat, condition, -0.5, 2)
      erp = average_epochs(epochs)
      save_object("$(subject)_$(condition)_epochs.jld2", epochs)
      save_object("$(subject)_$(condition)_erp.jld2", erp)
    end
  end
end

@time test_analysis()




function grand_average_erps(subjects, conditions)

  file_problem = check_files_exist(subjects, conditions, "erp")
  if file_problem
    return
  end

  sample_rate = nothing
  layout = nothing

  for condition in conditions

    ind_subject_condition_array = []

    for subject in subjects

      # load individual subject data
      erp = load_object("$(subject)_$(condition)_erp.jld2")

      # basic data checks to make sure sample layout and sample rate
      if isnothing(sample_rate)
        sample_rate = erp.sample_rate
      end
      if isnothing(layout)
        layout = erp.layout
      end
      if !isnothing(sample_rate)
        if sample_rate != erp.sample_rate
          throw(DomainError([sample_rate, erp.sample_rate], "sample rates across files do not match!"))
        end
      end
      if !isnothing(layout)
        if layout != erp.layout
          throw(DomainError([layout, erp.layout], "layout across files do not match!"))
        end
      end
      # update sample_rate/layout from current file
      sample_rate = erp.sample_rate
      layout = erp.layout

      # perform subject average
      for subject in subjects
        push!(ind_subject_condition_array, erp.data)
      end
      grand_average = reduce(.+, ind_subject_condition_array) ./ length(ind_subject_condition_array)
      save_object("$(cond)_ga_erp.jld2", grand_average)

    end
  end
end
grand_average_erps([3, 4], 1)




dat = read_bdf("../Flank_C_3.bdf")
dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
epochs = extract_epochs(dat, 1, -0.5, 2)
plot_epoch(epochs, 1:10, ["Cz", "CPz"], legend=false)


########################################################################


function data_interpolation_topo(dat, points; grid_scale=300)

  radius = 88 # mm
  x = y = range(-radius, radius, length=grid_scale)
  X, Y = repeat(x', grid_scale)[:], repeat(y', grid_scale)[:]
  grid = [X Y]'
  dat = interpolate(Multiquadratic(), points, dat)
  dat = ScatteredInterpolation.evaluate(dat, grid)
  dat = reshape(dat, grid_scale, grid_scale)
  circle_mask!(dat, grid_scale)
  return dat
end

















# average reference
# rereference!(eeg, 1:72);
# 
# # high-pass/low-pass filter
# highpass_filter!(eeg, 0.1, 2);
# lowpass_filter!(eeg, 30, 6);
# 
# # extract epochs and baseline
# epoched_data = extract_epochs(eeg, 1, -0.5, 2);
# baseline!(epoched_data, 0, 0)
# 
# # plot some epoched data
# plot_epoch(epoched_data, 1, [:PO7, :PO8]) # single trial/multiple eleectrodes
# plot_epoch(epoched_data, collect(1:10), [:PO7])   # multiple trials/single electrode
# 
# # average epochs
# erp = average_epochs(epoched_data);
# 
# # plot some erp data
# plot_erp(erp, [:Fp1, :Fp2])



# function read_mat_file(filename)
#   file = matopen(filename)
#   dat = read(file)
#   close(file)
#   return dat
# end
# 
# 
# struct InfoICA
#   topo
#   unmixing
#   label
# end
# 
# 
# function infomax_ica(dat; extended=true)
# 
#   # define some default values
#   l_rate = 0.001 # initial learning rate
#   max_iter = 200 # maximum number of iterations
#   w_change = 1e-12 # change to stop iteration
#   use_bias = true
#   anneal_deg = 60.0 # angle at which learning rate reduced
#   anneal_step = 0.9 # factor by which learning rate reduced
#   blowup = 1e4 # max difference allowed between two successive estimations of unmixing matrix
#   blowup_fac = 0.5 # factor by which learning rate will be reduced if "blowup"
#   n_small_angle = 20
#   max_weight = 1e8     # larger than this have "blowup"
#   restart_factor = 0.9 # if weights blowup, restart with lrate
#   degconst = 180 ./ pi
# 
#   # for extended Infomax
#   ext_blocks = 1
#   n_subgauss = 1
#   extmomentum = 0.5
#   signsbias = 0.02
#   signcount_threshold = 25
#   signcount_step = 2
#   kurt_size = 6000
# 
#   # check data shape and get block sizes
#   n_features, n_samples = size(dat)
#   n_features_square = n_features^2
#   block = trunc(Int, sqrt(n_samples / 3.0))
#   nblock = div(n_samples, block)
#   lastt = (nblock - 1) * block + 1
# 
#   # start ICA
#   @printf "Computing Infomax ICA\n"
# 
#   # initialize training
#   weights = Matrix{Float64}(I, n_features, n_features) # identity matrix
#   BI = block * Matrix{Float64}(I, n_features, n_features)
#   bias = zeros(n_features, 1)
#   onesrow = ones((1, block))
#   startweights = copy(weights)
#   oldweights = copy(startweights)
#   step = 0
#   count_small_angle = 0
#   wts_blowup = false
#   blockno = 0
#   signcount = 0
#   initial_ext_blocks = ext_blocks  # save the initial value in case of reset
# 
#   if extended
#     signs = ones(n_features)
#     signs[1:n_subgauss] .= -1
#     kurt_size = min(kurt_size, n_samples)
#     old_kurt = zeros(n_features)
#     oldsigns = zeros(n_features)
#   end
# 
#   # trainings loop
#   olddelta = 1.0
#   oldchange = 0.0
#   while step < max_iter
# 
#     permute = shuffle(1:n_samples)
#     permute = 1:n_samples
# 
#     for t = 1:block:lastt
# 
#       u = *(weights, dat[:, permute[t:t+block-1]]) .+ *(bias, onesrow)
# 
#       if extended
#         # extended ICA update
#         y = tanh.(u)
#         weights += l_rate * *(weights, BI .- signs[:] .* *(y, u') - *(u, u'))
#         if use_bias
#           bias += (l_rate .* sum(y, dims=2) .* -2.0)
#         end
#       else
#         y = 1 ./ (1 .+ exp.(-u))
#         weights += l_rate * *(weights, BI + *(u', (1.0 .- 2.0 .* y)))
#         if use_bias
#           bias += (l_rate .* sum((1.0 .- 2.0 .* y), dims=1))'
#         end
#       end
# 
#       # check change limit
#       if maximum(abs.(weights)) > max_weight
#         wts_blowup = true
#         break
#       end
# 
#       blockno += 1
#       # ICA kurtosis estimation
#       if extended
#         if ext_blocks > 0 & blockno % ext_blocks == 0
#           if kurt_size < n_samples
#             rp = trunc.(Int, (rand(kurt_size) .* (n_samples - 1))) .+ 1
#             tpartact = *(weights, dat[:, rp])'
#           else
#             tpartact = *(weights, dat)'
#           end
#           # estimate kurtosis
#           kurt = kurtosis.(eachcol(tpartact))
#           if extmomentum != 0
#             kurt = extmomentum .* old_kurt .+ (1.0 .- extmomentum) .* kurt
#             old_kurt = kurt
#           end
#           # estimate weighted signs
#           signs = sign.(kurt .+ signsbias)
#           ndiff = sum(signs .- oldsigns .!= 0)
#           ndiff == 0 ? signcount += 1 : signcount = 0
#           oldsigns = signs
#           if signcount >= signcount_threshold
#             ext_blocks = trunc.(ext_blocks .* signcount_step)
#             signcount = 0
#           end
#         end
#       end
# 
#     end
# 
#     if !wts_blowup
#       oldwtchange = weights .- oldweights
#       step += 1
#       angledelta = 0.0
#       delta = oldwtchange[:]
#       change = sum(delta .* delta)
#       if step > 2
#         angledelta = acos(sum(delta .* olddelta) / sqrt(change * oldchange))
#         angledelta *= degconst
#       end
# 
#       @printf "step %d: lrate %5f, wchange %8.8f, angledelta %4.1f\n" step l_rate change angledelta
# 
#       # anneal learning rate
#       oldweights = copy(weights)
#       if angledelta > anneal_deg
#         l_rate *= anneal_step  # anneal learning rate
#         # accumulate angledelta until anneal_deg reaches l_rate
#         olddelta = delta
#         oldchange = change
#         count_small_angle = 0  # reset count when angledelta is large
#       else
#         if step == 1  # on first step only
#           olddelta = delta  # initialize
#           oldchange = change
#         end
#         if !isnothing(n_small_angle)
#           count_small_angle += 1
#           if count_small_angle > n_small_angle
#             max_iter = step
#           end
#         end
#       end
# 
#       # apply stopping rule
#       if step > 2 && change < w_change
#         step = max_iter
#       elseif change > blowup
#         l_rate *= blowup_fac
#       end
# 
#       # restart if weights blow up (for lowering l_rate)
#     else
#       step = 0  # start again
#       wts_blowup = false  # re-initialize variables
#       blockno = 1
#       l_rate *= restart_factor  # with lower learning rate
#       weights = copy(startweights)
#       oldweights = copy(startweights)
#       olddelta = zeros((1, n_features_square))
#       bias = zeros((n_features, 1))
# 
#       ext_blocks = initial_ext_blocks
# 
#       # for extended Infomax
#       if extended
#         signs = ones(n_features)
#         signs[1:n_subgauss] .= -1
#         oldsigns = zeros(n_features)
#       end
# 
#       @printf "lowering learning rate to %g ... re-starting ...\n" l_rate
# 
#     end
#   end
# 
#   return weights
# 
# end
# 
# 
# 
# 
# 
# function pre_whiten(dat)
#   return dat ./ std(dat)
# end
# 
# 
# # pca
# function pca_reduction(dat)
#   dat = copy(dat)
#   dat .-= mean(dat, dims=1)
# 
#   U, S, V = LinearAlgebra.svd(dat, full=false)
#   max_abs_cols = argmax(abs.(U), dims=1)
#   signs = sign.(U[max_abs_cols])
#   U .*= signs
#   V .*= signs
# 
#   explained_variance = (S .^ 2) ./ (size(dat)[1] - 1)
#   # total_var = sum(explained_variance)
#   # explained_variance_ratio_ = explained_variance / total_var
# 
#   U .*= sqrt(size(dat)[1] - 1)
# 
#   return U, explained_variance, V
# 
# end





dat = read_bdf("../Flank_C_3.bdf")
dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
filter_data!(dat, "hp", 2, 2)
rereference!(dat, dat.layout.label, collect(1:72));

# dat = Matrix(dat.data[:, 3:end])
# 
# dat = pre_whiten(dat)
# 
# dat = read_mat_file("../dat_ica.mat")
# dat = dat["dat_ica"]
# dat = dat'
# 
# dat, explained_variance = pca_reduction(dat)
# 
# @time unmixing_matrix = infomax_ica(dat[1:20000, 1:68]', extended=true)
# @time unmixing_matrix = infomax_ica(dat', extended=true)

using MultivariateStats

model = fit(ICA, Matrix(dat.data[:, 3:end])', 71, maxiter=2, tol=0.1)

ic_mw = 68 == size(dat', 1) ? inv(model.W)' : pinv(model.W)'
ic = MultivariateStats.predict(model, dat')

var(ic, dims=2)

# f = Figure()
# ax = GLMakie.Axis(f[1, 1])
# for i in 1:71
#   GLMakie.lines!(unmixing_matrix[:, i])
# end

stable = explained_variance ./ explained_variance[1] .> 1e-6
norms = explained_variance[1:10]
norms = sqrt.(norms)
norms[norms.==0] .= 1.0

unmixing_matrix ./= norms

mixing_matrix = LinearAlgebra.pinv(unmixing_matrix)

# sort
source_data = dat[1:20000, 1:10]
var = sum(mixing_matrix^2, dims=1) .* sum(source_data .^ 2, dims=1) / (10 * 20000 - 1)
var ./= sum(var)

order = sortperm(var, dims=2, rev=true)
var[order]

unmixing_matrix = unmixing_matrix[order[:], :]
mixing_matrix = mixing_matrix[order[:], :]

lines(mixing_matrix[1, :])

unmixing_matrix = unmixing_matrix * eye(69, 69) * eigen(:)






