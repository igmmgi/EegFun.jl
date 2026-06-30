using CairoMakie
using DelimitedFiles

function plot_comparison()
    if length(ARGS) < 1
        error("data directory requiored as a command line argument (e.g., julia plot_benchmark_comparison.jl /path/to/data)")
    end
    data_dir = ARGS[1]

    if !isdir(data_dir)
        error("Data directory not found: $data_dir")
    end

    results_dir = joinpath(data_dir, "benchmarks")

    # Define paths based on results_dir
    python_file = joinpath(results_dir, "python_data.csv")
    matlab_file = joinpath(results_dir, "matlab_data.csv")
    julia_file = joinpath(results_dir, "julia_data.csv")

    # Load data (skipping the header row)
    python_data = isfile(python_file) ? readdlm(python_file, ',', skipstart = 1) : missing
    matlab_data = isfile(matlab_file) ? readdlm(matlab_file, ',', skipstart = 1) : missing
    julia_data = isfile(julia_file) ? readdlm(julia_file, ',', skipstart = 1) : missing

    # Load timings if they exist
    time_file_j = joinpath(results_dir, "julia_time.txt")
    time_file_p = joinpath(results_dir, "python_time.txt")
    time_file_m = joinpath(results_dir, "matlab_time.txt")

    julia_time = isfile(time_file_j) ? parse(Float64, read(time_file_j, String)) : NaN
    python_time = isfile(time_file_p) ? parse(Float64, read(time_file_p, String)) : NaN
    matlab_time = isfile(time_file_m) ? parse(Float64, read(time_file_m, String)) : NaN

    lbl_julia = isnan(julia_time) ? "Julia" : "Julia [$(round(julia_time, digits=2))s]"
    lbl_python = isnan(python_time) ? "Python" : "MNE-Python [$(round(python_time, digits=2))s]"
    lbl_matlab = isnan(matlab_time) ? "MATLAB" : "EEGLAB [$(round(matlab_time, digits=2))s]"

    # Initialize combined figure
    fig_combined = Figure(size = (800, 600))
    ax1 = Axis(
        fig_combined[1, 1],
        xlabel = "Time (s)",
        ylabel = "Amplitude (µV)",
        title = "Combined Cross-Pipeline Comparison (PO7/PO8)",
        xgridvisible = false,
        ygridvisible = false,
    )

    # Initialize individual figure
    fig_indiv = Figure(size = (800, 900))
    ax2 = Axis(fig_indiv[1, 1], ylabel = "Amplitude (µV)", title = "Julia (EegFun.jl)", xgridvisible = false, ygridvisible = false)
    ax3 = Axis(fig_indiv[2, 1], ylabel = "Amplitude (µV)", title = "Python (MNE)", xgridvisible = false, ygridvisible = false)
    ax4 = Axis(
        fig_indiv[3, 1],
        xlabel = "Time (s)",
        ylabel = "Amplitude (µV)",
        title = "MATLAB (EEGLAB)",
        xgridvisible = false,
        ygridvisible = false,
    )

    all_axes = [ax1, ax2, ax3, ax4]

    # Add reference lines and x limits to all axes
    for ax in all_axes
        vlines!(ax, [0.0], color = :black, linewidth = 1)
        hlines!(ax, [0.0], color = :black, linewidth = 1)
        xlims!(ax, -0.2, 1.0)
    end

    # Colors for pipelines
    c_julia = :blue
    c_python = :red
    c_matlab = :green

    # Plot Julia
    if !ismissing(julia_data)
        # Combined plot
        lines!(ax1, julia_data[:, 1], julia_data[:, 2], color = c_julia, linewidth = 2, label = "$lbl_julia (Valid)")
        lines!(ax1, julia_data[:, 1], julia_data[:, 3], color = c_julia, linestyle = :dash, linewidth = 2, label = "$lbl_julia (Invalid)")
        # Julia only plot
        lines!(ax2, julia_data[:, 1], julia_data[:, 2], color = c_julia, linewidth = 2, label = "Valid")
        lines!(ax2, julia_data[:, 1], julia_data[:, 3], color = c_julia, linestyle = :dash, linewidth = 2, label = "Invalid")
    end

    # Plot Python
    if !ismissing(python_data)
        # Combined plot
        lines!(ax1, python_data[:, 1], python_data[:, 2], color = c_python, linewidth = 2, label = "$lbl_python (Valid)")
        lines!(
            ax1,
            python_data[:, 1],
            python_data[:, 3],
            color = c_python,
            linestyle = :dash,
            linewidth = 2,
            label = "$lbl_python (Invalid)",
        )
        # Python only plot
        lines!(ax3, python_data[:, 1], python_data[:, 2], color = c_python, linewidth = 2, label = "Valid")
        lines!(ax3, python_data[:, 1], python_data[:, 3], color = c_python, linestyle = :dash, linewidth = 2, label = "Invalid")
    end

    # Plot MATLAB
    if !ismissing(matlab_data)
        # Combined plot
        lines!(ax1, matlab_data[:, 1], matlab_data[:, 2], color = c_matlab, linewidth = 2, label = "$lbl_matlab (Valid)")
        lines!(
            ax1,
            matlab_data[:, 1],
            matlab_data[:, 3],
            color = c_matlab,
            linestyle = :dash,
            linewidth = 2,
            label = "$lbl_matlab (Invalid)",
        )
        # MATLAB only plot
        lines!(ax4, matlab_data[:, 1], matlab_data[:, 2], color = c_matlab, linewidth = 2, label = "Valid")
        lines!(ax4, matlab_data[:, 1], matlab_data[:, 3], color = c_matlab, linestyle = :dash, linewidth = 2, label = "Invalid")
    end

    # Link axes and add legends
    linkaxes!(ax2, ax3, ax4)

    for ax in all_axes
        try
            axislegend(ax, position = :rt)
        catch
            # Ignore if there are no plots with labels in this axis
        end
    end

    # Save to PDF
    out_file1 = joinpath(results_dir, "cross_pipeline_comparison.pdf")
    save(out_file1, fig_combined)

    out_file2 = joinpath(results_dir, "individual_pipelines.pdf")
    save(out_file2, fig_indiv)
    println("Successfully generated plots: \n  - ", out_file1, "\n  - ", out_file2)
end

plot_comparison()
