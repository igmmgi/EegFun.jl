using CairoMakie
using DelimitedFiles

function plot_comparison()
    if length(ARGS) < 1
        error("Please provide the data directory as a command line argument (e.g., julia plot_benchmark_comparison.jl /path/to/data)")
    end
    data_dir = ARGS[1]
    
    if !isdir(data_dir)
        error("Data directory not found: $data_dir")
    end

    # Define paths based on data_dir
    python_file = joinpath(data_dir, "benchmark_python_data.csv")
    matlab_file = joinpath(data_dir, "benchmark_matlab_data.csv")
    julia_file = joinpath(data_dir, "benchmark_julia_data.csv")

    # Load data (skipping the header row)
    python_data = isfile(python_file) ? readdlm(python_file, ',', skipstart=1) : missing
    matlab_data = isfile(matlab_file) ? readdlm(matlab_file, ',', skipstart=1) : missing
    julia_data = isfile(julia_file) ? readdlm(julia_file, ',', skipstart=1) : missing
    
    # Load timings if they exist
    time_file_j = joinpath(data_dir, "benchmark_julia_time.txt")
    time_file_p = joinpath(data_dir, "benchmark_python_time.txt")
    time_file_m = joinpath(data_dir, "benchmark_matlab_time.txt")
    
    julia_time = isfile(time_file_j) ? parse(Float64, read(time_file_j, String)) : NaN
    python_time = isfile(time_file_p) ? parse(Float64, read(time_file_p, String)) : NaN
    matlab_time = isfile(time_file_m) ? parse(Float64, read(time_file_m, String)) : NaN
    
    lbl_julia = isnan(julia_time) ? "Julia" : "Julia [$(round(julia_time, digits=2))s]"
    lbl_python = isnan(python_time) ? "Python" : "MNE-Python [$(round(python_time, digits=2))s]"
    lbl_matlab = isnan(matlab_time) ? "MATLAB" : "EEGLAB [$(round(matlab_time, digits=2))s]"

    # Initialize figure
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1],
        xlabel="Time (s)",
        ylabel="Amplitude (µV)",
        title="Cross-Pipeline ERP Comparison (PO7/PO8)",
        xgridvisible=false,
        ygridvisible=false
    )

    # Add reference lines (x=0, y=0)
    vlines!(ax, [0.0], color=:black, linewidth=1)
    hlines!(ax, [0.0], color=:black, linewidth=1)

    # Colors for pipelines
    c_julia = :blue
    c_python = :red
    c_matlab = :green

    # Plot Julia
    if !ismissing(julia_data)
        lines!(ax, julia_data[:, 1], julia_data[:, 2], color=c_julia, linewidth=2, label="$lbl_julia (Valid)")
        lines!(ax, julia_data[:, 1], julia_data[:, 3], color=c_julia, linestyle=:dash, linewidth=2, label="$lbl_julia (Invalid)")
    end

    # Plot Python
    if !ismissing(python_data)
        lines!(ax, python_data[:, 1], python_data[:, 2], color=c_python, linewidth=2, label="$lbl_python (Valid)")
        lines!(ax, python_data[:, 1], python_data[:, 3], color=c_python, linestyle=:dash, linewidth=2, label="$lbl_python (Invalid)")
    end

    # Plot MATLAB
    if !ismissing(matlab_data)
        lines!(ax, matlab_data[:, 1], matlab_data[:, 2], color=c_matlab, linewidth=2, label="$lbl_matlab (Valid)")
        lines!(ax, matlab_data[:, 1], matlab_data[:, 3], color=c_matlab, linestyle=:dash, linewidth=2, label="$lbl_matlab (Invalid)")
    end

    # Set x limits to match the standard bounds
    xlims!(ax, -0.2, 1.0)

    # Add legend
    axislegend(ax, position=:rt)

    # Save to PDF
    out_file = joinpath(data_dir, "cross_pipeline_comparison.pdf")
    save(out_file, fig)
    println("Successfully generated cross-pipeline comparison plot: ", out_file)
end

plot_comparison()
