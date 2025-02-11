using BioSemiBDF
using CSV
using DSP
using DataFrames
using GLMakie
# using CairoMakie
using JLD2
using LibGEOS
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
include("ica.jl")

# include("test/runtests.jl")
# test_baseline()
# test_filter()

# using Logging
# # Show all messages
# global_logger(ConsoleLogger(stderr, Logging.Debug))
# # Show only info and above
# global_logger(ConsoleLogger(stderr, Logging.Info))
# # Show only warnings and errors
# global_logger(ConsoleLogger(stderr, Logging.Warn))

# basic layouts
layout = read_layout("./layouts/biosemi72.csv");
polar_to_cartesian_xy!(layout)
polar_to_cartesian_xyz!(layout)
neighbours, nneighbours = get_electrode_neighbours_xy(layout, 50)
neighbours, nneighbours = get_electrode_neighbours_xyz(layout, 40)

#head_shape_2d(layout)
#head_shape_2d(layout, point_kwargs = Dict(:markersize => 30), label_kwargs = Dict(:fontsize => 30, :xoffset => 1))
#layout = filter(row -> row.label in ["PO7", "PO8"], layout)
head_shape_2d(layout)
#head_shape_3d(layout);

# read bdf file
subject = 3
dat = read_bdf("../Flank_C_$(subject).bdf");

# save / load
# save_object("$(subject)_continuous_raw.jld2", dat)
# dat = load_object("3_continuous_raw.jld2")

# plot_events(dat)
dat = create_eeg_dataframe(dat, layout);
# save_object("$(subject)_continuous_raw_eegfun.jld2", dat)
# dat = load_object("3_continuous_raw.jld2")

# rereference!(dat.data, dat.layout.label, :Fp1)
rereference!(dat.data, dat.layout.label, dat.layout.label)
plot_databrowser(dat)

filter_data!(dat, "hp", 0.1, 2)
plot_databrowser(dat)


# plot_events(dat)

# # # search for some bad channels
channel_data = channel_summary(dat.data, dat.layout.label[1:66])
channel_data = channel_summary(dat.data, dat.layout.label)


# # # bad channels zscore variance
# # c[!, :channel][c[!, :zvar].>3]
# c = channel_joint_probability(dat.data, dat.layout.label[1:66])
# # lines(c[!, :jp])
cm = correlation_matrix(dat.data, dat.layout.label)
plot_correlation_heatmap(cm)
plot_correlation_heatmap(cm, (-0.2, 0.2))

# TODO: add neighbour correlation values

# filter_data!(dat, "lp", 10, 6)
# include("plot.jl")
# calculate EOG channels
diff_channel!(dat, ["Fp1", "Fp2"], ["IO1", "IO2"], "vEOG");
diff_channel!(dat, "F9", "F10", "hEOG");
## # autodetect EOG signals
detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)

# TODO: this seems a bit inconsistent with functions above
dat.data[!, "is_extreme"] .= is_extreme_value(dat.data, dat.layout.label, 500);


# ICA "continuous" data
dat_ica = filter_data(dat, "hp", 1, 2)

dat_ica = create_ica_data_matrix(dat_ica)
dat_ica = create_ica_data_matrix(dat_ica; channels_to_include = setdiff(dat_ica.layout.label, ["PO9"]))
# dat_ica = create_ica_data_matrix(dat, dat.layout.label, 1:10)


ica_result = infomax_ica(dat_ica, setdiff(dat.layout.label, ["PO9"]), n_components = 70)
plot_ica_topoplot(ica_result, dat.layout)

# or with bad electrode removed
tmp = dat_ica.data[dat_ica.data[!, "is_extreme"].==false, [4:65; 67:75]]
ica_result = infomax_ica(permutedims(Float64.(Matrix(tmp))), dat.layout.label[[1:64; 66:72]], n_components = 10)
plot_ica_topoplot(ica_result, dat.layout)


# Continuous Data Browser
# TODO: Labels position when changing x-range
# TODO: Improve logic of plotting marker (triggers/EOG) lines?
# plot_databrowser(dat_ica)
# plot_databrowser(dat, [dat.layout.label; "hEOG"; "vEOG"])
# plot_databrowser(dat, ["vEOG", "hEOG"])
# plot_databrowser(dat, "hEOG")

# extract epochs
epochs = []
for (idx, epoch) in enumerate([1, 4, 5, 3])
    push!(epochs, extract_epochs(dat_ica, idx, epoch, -0.5, 2))
end

# concat data from all epochs
tmp = vcat([vcat(epochs[i].data[:]...) for i in eachindex(epochs)]...)

output = infomax_ica(permutedims(Float64.(Matrix(tmp[:, 6:77]))), dat.layout.label, n_components = 71)
plot_ica_topoplot(output, dat.layout)

# epochs_cleaned = remove_bad_epochs(epochs)

# Epoch Data Browser
plot_databrowser(epochs[1])
plot_databrowser(epochs[1], [epochs[1].layout.label; "hEOG"; "vEOG"])

plot_epochs(epochs[1], [:Fp1])
plot_epochs(epochs[1], :Fp1)
# plot_epochs(epochs, "Fp1")
plot_epochs(epochs[1], ["PO7", "PO8"])

# average epochs
erps = []
for (idx, epoch) in enumerate(epochs)
    push!(erps, average_epochs(epochs[idx]))
end

# ERP Plot
# f, ax = plot_erp(erp, :Fp1)
plot_erp(erps[1], ["Fp1", "Fp2"])
# plot_erp(erp, [:Fp1, :Fp2])
# plot_erp(erp, ["Fp1", "Fp2"])
plot_erp(erp, [:Fp1, :Fp2], kwargs = Dict(:yreversed => true))

# Topoplot
include("topo.jl")
plot_topoplot(erp, xlim = [0.3, 0.3], ylim = [-10, 10])
plot_topoplot(erp; method = :spherical_splines)

# ERP Image
plot_erp_image(epochs, :Fp1)
plot_erp_image(epochs, "Fp1")
plot_erp_image(epochs, [:Fp1, :Fp2], colorrange = [-50, 50])


save_object("$(subject)_$(cond)_epochs.jld2", epochs)
save_object("$(subject)_$(cond)_erp.jld2", erp)

diff_channel!(dat, "F9", "F10", "hEOG");
diff_channel!(dat, ["Fp1", "Fp2"], ["IO1", "IO2"], "vEOG");


detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)

test_plot_eog_detection(dat, 1000:14000, "vEOG", "is_vEOG")
test_plot_eog_detection(dat, 1000:4000, "hEOG", "is_hEOG")



function test_analysis()
    for subject = 3:4
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




########################################################################



















