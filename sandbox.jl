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

using Logging
# Show all messages
global_logger(ConsoleLogger(stderr, Logging.Debug))
# Show only info and above
global_logger(ConsoleLogger(stderr, Logging.Info))
# Show only warnings and errors
global_logger(ConsoleLogger(stderr, Logging.Warn))


# basic layouts
layout = read_layout("./layouts/biosemi72.csv");
head_shape_2d(layout);
# head_shape_3d(layout);

# read bdf file
subject = 3
dat = read_bdf("../Flank_C_$(subject).bdf");
dat = create_eeg_dataframe(dat, layout);
# basic bdf plot
# plot_databrowser(dat)
filter_data!(dat, "hp", 1, 2)
# filter_data!(dat, "lp", 10, 6)
# include("plot.jl")
# calculate EOG channels
diff_channel!(dat, ["Fp1", "Fp2"], ["IO1", "IO2"], "vEOG");
diff_channel!(dat, "F9", "F10", "hEOG");
## # autodetect EOG signals
detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)
dat.data[!, "is_extreme"] .= is_extreme_value(dat.data, dat.layout.label, 100);

# data_whitened = pre_whiten(Float64.(transpose(Matrix(dat.data[!, 3:end-4]))))
# Run ICA
# weights = infomax_ica(data_whitened, extended=true)
# Get independent components
# components = weights * data_whitened


# Continuous Data Browser
# TODO: Labels position when changing x-range
# TODO: Improve logic of plotting marker (triggers/EOG) lines?
# plot_databrowser(dat)
# plot_databrowser(dat, [dat.layout.label; "hEOG"; "vEOG"])
# plot_databrowser(dat, ["vEOG", "hEOG"])
# plot_databrowser(dat, "hEOG")

# extract epochs
epochs = extract_epochs(dat, 1, -0.5, 2)

# Epoch Data Browser
# plot_databrowser(epochs)
# plot_databrowser(epochs, [epochs.layout.label; "hEOG"; "vEOG"])
# plot_databrowser(epochs, ["hEOG", "vEOG"])
# plot_databrowser(epochs, "hEOG")

# # Plot Epochs (all)
# plot_epochs(epochs, :Fp1)
# plot_epochs(epochs, "Fp1")
# plot_epochs(epochs, ["PO7", "PO8"])

# average epochs
erp = average_epochs(epochs)

# ERP Plot
# f, ax = plot_erp(erp, :Fp1)
# plot_erp(erp, ["Fp1"])
# plot_erp(erp, [:Fp1, :Fp2])
# plot_erp(erp, ["Fp1", "Fp2"])
# plot_erp(erp, [:Fp1, :Fp2], kwargs = Dict(:yreversed => true))

# Topoplot
include("topo.jl")
plot_topoplot(erp )
plot_topoplot(erp; method = :spherical_splines )

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

function test_plot_eog_detection(dat, xlim, channel, detected)
    fig = Figure()
    ax = Axis(fig[1, 1])  # plot layout
    lines!(ax, dat.data.time[xlim], dat.data[!, channel][xlim])
    vlines!(ax, dat.data.time[xlim][dat.data[!, detected][xlim]], color = :black)
    display(fig)
    return fig, ax
end

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




dat = read_bdf("../Flank_C_3.bdf")
dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
epochs = extract_epochs(dat, 1, -0.5, 2)
plot_epoch(epochs, 1:10, ["Cz", "CPz"], legend = false)


########################################################################



















# function read_mat_file(filename)
#   file = matopen(filename)
#   dat = read(file)
#   close(file)
#   return dat
# end
