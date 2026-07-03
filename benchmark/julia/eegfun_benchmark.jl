using EegFun
using BenchmarkTools
using CairoMakie
using DelimitedFiles
using Statistics
CairoMakie.activate!()

function process_eeg_to_grandaverage(raw_files::Vector{String}, data_dir::String, run_ica::Bool)

    all_subject_erps = EegFun.ErpData[]

    # Define the 4 epoch conditions as specified in the Visual Attention tutorial
    epoch_conditions = [
        EegFun.EpochCondition(name = "Valid", trigger_sequences = [[1, 5], [3, 5]], reference_index = 2),
        EegFun.EpochCondition(name = "Invalid", trigger_sequences = [[2, 5], [4, 5]], reference_index = 2),
    ]

    individual_times = Float64[]
    for (i, file) in enumerate(raw_files)
        file_time = @elapsed begin
        # 1. Read data, layout, and create EegFun.jl data structure
        raw_data = EegFun.read_raw_data(file)

        layout_path = joinpath(@__DIR__, "..", "..", "resources", "layouts", "biosemi", "biosemi72.csv")
        layout = EegFun.read_layout(layout_path)
        EegFun.polar_to_cartesian_xy!(layout)

        dat = EegFun.create_eegfun_data(raw_data, copy(layout))

        # 2. Initial preprocessing options: rereference, high-pass filter, and mark bad data sections
        EegFun.rereference!(dat, :avg)
        EegFun.highpass_filter!(dat, 0.1)
        EegFun.is_extreme_value!(dat, 250.0, channel_out = :bad_segments)

        if run_ica
            # 3. Apply extended Infomax ICA (additional high-pass filter at 1 Hz)
            dat_ica = deepcopy(dat)
            EegFun.highpass_filter!(dat_ica, 1.0)
            ica_result = EegFun.run_ica(
                dat_ica;
                algorithm = :infomax_extended,
                n_components = EegFun.n_channels(dat_ica) - 1,
                sample_selection = EegFun.samples_not(:bad_segments),
            )

            # Sanity check: Plot the first 30 ICA components for the first subject
            if i == 1
                fig_ica = EegFun.plot_topography(
                    ica_result,
                    component_selection = EegFun.components(1:30),
                    colorbar_plot = false,
                    label_plot = false,
                    point_plot = false,
                )
                resize!(fig_ica.fig, 1200, 1000)
                save(joinpath(results_dir, "julia_ica.pdf"), fig_ica.fig)
            end

            # 4. Remove ICA component
            EegFun.remove_ica_components!(dat, ica_result, component_selection = EegFun.components([1]))
        end

        # 5. Epoch and baseline data
        epochs = EegFun.extract_epochs(dat, epoch_conditions, (-0.5, 1.0))
        EegFun.baseline!(epochs, (-0.2, 0.0))

        # Average epochs        
        participant_erps = EegFun.average_epochs(epochs)

        append!(all_subject_erps, participant_erps)
        end
        push!(individual_times, file_time)
    end

    # 6. Lowpass filter before grand averaging
    EegFun.lowpass_filter!.(all_subject_erps, 30.0)

    # 7. Grand Average (or just use the single subject if n=1)
    grand_avg = length(raw_files) == 1 ? all_subject_erps : EegFun.grand_average(all_subject_erps)

    # 8. Plot result
    fig_res = EegFun.plot_erp(
        grand_avg;
        channel_selection = EegFun.channels([:PO7, :PO8]),
        average_channels = true,
        interval_selection = (-0.2, 1.0),
        legend_position = :rt,
        legend_labels = ["Valid", "Invalid"],
    )
    save(joinpath(results_dir, "julia_erp.pdf"), fig_res.fig)

    # 9. Export data to CSV
    valid_po = EegFun.subset(grand_avg[1]; channel_selection = EegFun.channels([:PO7, :PO8]))
    invalid_po = EegFun.subset(grand_avg[2]; channel_selection = EegFun.channels([:PO7, :PO8]))
    valid_data = vec(mean(Matrix(valid_po.data[!, [:PO7, :PO8]]), dims = 2))
    invalid_data = vec(mean(Matrix(invalid_po.data[!, [:PO7, :PO8]]), dims = 2))

    open(joinpath(results_dir, "julia_data.csv"), "w") do io
        write(io, "time,valid,invalid\n")
        writedlm(io, [valid_po.data.time valid_data invalid_data], ',')
    end

    return grand_avg, individual_times
end

if length(ARGS) < 1
    error("Please provide the data directory as a command line argument (e.g., julia eegfun_benchmark.jl /path/to/data)")
end
data_dir = ARGS[1]

if !isdir(data_dir)
    error("Data directory not found: $data_dir")
end
test_files = sort(filter(f -> endswith(lowercase(f), ".bdf"), readdir(data_dir, join = true)))

n_files_to_process = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 0
if n_files_to_process > 0
    test_files = test_files[1:min(length(test_files), n_files_to_process)]
end

run_ica_flag = length(ARGS) >= 3 ? parse(Bool, lowercase(ARGS[3])) : true

results_dir = joinpath(data_dir, "benchmarks")
mkpath(results_dir)

# println("Compiling Julia pipeline (first run)...")
warmup_files = min(length(test_files), 2)
Base.invokelatest(process_eeg_to_grandaverage, test_files[1:warmup_files], data_dir, run_ica_flag) # Warm-up to trigger JIT compilation

println("Benchmarking...")
val, t_elapsed = @timed Base.invokelatest(process_eeg_to_grandaverage, test_files, data_dir, run_ica_flag)
grand_avg, individual_times = val

# Save execution time
import Pkg
eegfun_version = Pkg.project().version

open(joinpath(results_dir, "julia_time.txt"), "w") do io
    println(io, t_elapsed)
    println(io, "Julia Version: ", VERSION)
    println(io, "EegFun Version: ", eegfun_version)
    println(io, "Individual Times: ", join(round.(individual_times, digits=2), ", "))
end
println("Julia execution time: ", round(t_elapsed, digits = 2), " seconds")
