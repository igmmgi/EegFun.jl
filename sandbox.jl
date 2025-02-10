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
#head_shape_2d(layout)
#head_shape_2d(layout, point_kwargs = Dict(:markersize => 30), label_kwargs = Dict(:fontsize => 30, :xoffset => 1))
#layout = filter(row -> row.label in ["PO7", "PO8"], layout)
#head_shape_2d(layout)
#head_shape_3d(layout);

# read bdf file
subject = 3
dat = read_bdf("../Flank_C_$(subject).bdf");
plot_events(dat)
dat = create_eeg_dataframe(dat, layout);
filter_data!(dat, "hp", 0.1, 2)
plot_events(dat)



# basic bdf plot
# plot_databrowser(dat)
rereference!(dat.data, dat.layout.label, dat.layout.label)
# rereference!(dat.data, dat.layout.label, :Fp1)
# plot_databrowser(dat)

# filter_data!(dat, "lp", 10, 6)
# include("plot.jl")
# calculate EOG channels
diff_channel!(dat, ["Fp1", "Fp2"], ["IO1", "IO2"], "vEOG");
diff_channel!(dat, "F9", "F10", "hEOG");
## # autodetect EOG signals
detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)

dat.data[!, "is_extreme"] .= is_extreme_value(dat.data, dat.layout.label, 100);

dat_ica = filter_data(dat, "hp", 1, 2)
output = infomax_ica(permutedims(Float64.(Matrix(dat_ica.data[:, 3:74]))), dat.layout.label, n_components = 71)





plot_ica_topoplot(output, dat.layout, ncomps = 4)

function plot_ica_topoplot(
    ica,
    layout;
    ncomps = nothing,
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict(),
    topo_kwargs = Dict(),
    colorbar_kwargs = Dict(),
)
    if (:x2 ∉ names(layout) || :y2 ∉ names(layout))
        polar_to_cartesian_xy!(layout)
    end
    if isnothing(ncomps)
        ncomps = size(ica.mixing)[2]
    end
    head_default_kwargs = Dict(:color => :black, :linewidth => 2)
    head_kwargs = merge(head_default_kwargs, head_kwargs)
    point_default_kwargs = Dict(:plot_points => false, :marker => :circle, :markersize => 12, :color => :black)
    point_kwargs = merge(point_default_kwargs, point_kwargs)
    label_default_kwargs =
        Dict(:plot_labels => false, :fontsize => 20, :color => :black, :color => :black, :xoffset => 0, :yoffset => 0)
    label_kwargs = merge(label_default_kwargs, label_kwargs)
    xoffset = pop!(label_kwargs, :xoffset)
    yoffset = pop!(label_kwargs, :yoffset)
    topo_default_kwargs = Dict(:colormap => :jet, :gridscale => 300)
    topo_kwargs = merge(topo_default_kwargs, topo_kwargs)
    gridscale = pop!(topo_kwargs, :gridscale)
    colorbar_default_kwargs = Dict(:plot_colorbar => true, :width => 30)
    colorbar_kwargs = merge(colorbar_default_kwargs, colorbar_kwargs)
    plot_colorbar = pop!(colorbar_kwargs, :plot_colorbar)
    fig = Figure()
    dims = best_rect(ncomps)
    count = 1
    axs = []
    for dim1 = 1:dims[1]
        for dim2 = 1:dims[2]
            ax = Axis(fig[dim1, dim2])
            push!(axs, ax)
            count += 1
            if count > ncomps
                break
            end
        end
    end
    count = 1
    for ax in axs
        ax.title = ica.ica_label[count]
        data = data_interpolation_topo(ica.mixing[:, count], permutedims(Matrix(layout[!, [:x2, :y2]])), gridscale)
        gridscale = gridscale
        radius = 88 # mm
        co = contourf!(
            ax,
            range(-radius * 2, radius * 2, length = gridscale),
            range(-radius * 2, radius * 2, length = gridscale),
            data,
            colormap = :jet,
        )
        # TODO: improve colorbar stuff
        # if plot_colorbar
        #     Colorbar(ax, co; colorbar_kwargs...)
        # end
        # head shape
        head_shape_2d(
            fig,
            ax,
            layout,
            head_kwargs = head_kwargs,
            point_kwargs = point_kwargs,
            label_kwargs = label_kwargs,
        )
        count += 1
        if count > ncomps
            break
        end
    end
    return fig
end

# interpolate data







# size(dat_ica.data)

# data_whitened = pre_whiten(Float64.(transpose(Matrix(dat.data[!, 3:end-4]))))
# Run ICA
# weights = infomax_ica(data_whitened, extended=true)
# Get independent components
# components = weights * data_whitened


# Continuous Data Browser
# TODO: Labels position when changing x-range
# TODO: Improve logic of plotting marker (triggers/EOG) lines?
plot_databrowser(dat_ica)
# plot_databrowser(dat, [dat.layout.label; "hEOG"; "vEOG"])
# plot_databrowser(dat, ["vEOG", "hEOG"])
# plot_databrowser(dat, "hEOG")

# extract epochs
epochs = extract_epochs(dat, 1, -0.5, 2)
epochs_cleaned = remove_bad_epochs(epochs)

# Epoch Data Browser
plot_databrowser(epochs)
plot_databrowser(epochs, [epochs.layout.label; "hEOG"; "vEOG"])
plot_databrowser(epochs, ["hEOG", "vEOG"])
# plot_databrowser(epochs, "hEOG")

plot_epochs(epochs, [:Fp1])
# # Plot Epochs (all)
plot_epochs(epochs, :Fp1)
# plot_epochs(epochs, "Fp1")
plot_epochs(epochs, ["PO7", "PO8"])

# average epochs
erp = average_epochs(epochs)
erp_cleaned = average_epochs(epochs_cleaned)

plot_erp(erp)
plot_erp(erp_cleaned)

# ERP Plot
# f, ax = plot_erp(erp, :Fp1)
plot_erp(erp, ["Fp1"])
# plot_erp(erp, [:Fp1, :Fp2])
# plot_erp(erp, ["Fp1", "Fp2"])
# plot_erp(erp, [:Fp1, :Fp2], kwargs = Dict(:yreversed => true))

# Topoplot
include("topo.jl")
plot_topoplot(erp)
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




dat = read_bdf("../Flank_C_3.bdf")
dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
epochs = extract_epochs(dat, 1, -0.5, 2)
plot_epoch(epochs, 1:10, ["Cz", "CPz"], legend = false)


########################################################################



















