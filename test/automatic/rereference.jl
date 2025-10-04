using Test
using DataFrames
using eegfun
using JLD2
using Statistics



@testset "rereference" begin
    # Layout and data
    layout_df = DataFrame(label = [:A, :B, :M1, :M2], inc = [0.0, 0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0, 0.0])
    layout = eegfun.Layout(layout_df, nothing, nothing)
    n = 5
    t = collect(0:0.001:((n-1)*0.001))
    df = DataFrame(time = t, sample = 1:n)
    df.A = [1, 2, 3, 4, 5]
    df.B = [2, 4, 6, 8, 10]
    df.M1 = [1, 1, 1, 1, 1]
    df.M2 = [3, 3, 3, 3, 3]
    dat = eegfun.ContinuousData(copy(df, copycols = true), layout, 1000, eegfun.AnalysisInfo())

    # 1) Average reference (:avg) subtracts mean of A,B,M1,M2 from each selected channel
    dat_avg = copy(dat)
    eegfun.rereference!(dat_avg, :avg, eegfun.channels([:A, :B]))
    avg_ref = (df.A .+ df.B .+ df.M1 .+ df.M2) ./ 4
    @test all(dat_avg.data.A .== df.A .- avg_ref)
    @test all(dat_avg.data.B .== df.B .- avg_ref)
    @test eegfun.reference(dat_avg.analysis_info) == :avg

    # 2) Mastoid reference uses M1 and M2
    dat_mast = copy(dat)
    eegfun.rereference!(dat_mast, :mastoid, eegfun.channels([:A, :B]))
    mast_ref = (df.M1 .+ df.M2) ./ 2
    @test all(dat_mast.data.A .== df.A .- mast_ref)
    @test all(dat_mast.data.B .== df.B .- mast_ref)
    @test eegfun.reference(dat_mast.analysis_info) == :mastoid

    # 3) Specific reference [A] applied to B only
    dat_spec = copy(dat)
    eegfun.rereference!(dat_spec, [:A], eegfun.channels([:B]))
    @test all(dat_spec.data.B .== df.B .- df.A)
    @test eegfun.reference(dat_spec.analysis_info) == Symbol("A")

    # 4) EpochData: per-epoch reference computed independently
    df1 = copy(df);
    df2 = copy(df)
    df1o = copy(df1);
    df2o = copy(df2)
    df1.epoch = fill(1, n);
    df2.epoch = fill(2, n)
    ep = eegfun.EpochData([df1, df2], layout, 1000, eegfun.AnalysisInfo())
    eegfun.rereference!(ep, :mastoid, eegfun.channels([:A, :B]))
    @test all(ep.data[1].A .== df1o.A .- ((df1o.M1 .+ df1o.M2) ./ 2))
    @test all(ep.data[2].B .== df2o.B .- ((df2o.M1 .+ df2o.M2) ./ 2))

    # 5) Non-mutating version returns a new object; original unchanged
    dat_nm = eegfun.rereference(dat, :mastoid, eegfun.channels([:A]))
    @test :A ∈ propertynames(dat_nm.data) && :A ∈ propertynames(dat.data)
    @test !all(dat_nm.data.A .== dat.data.A)

    # 6) Only selected channels modified (B, M1, M2 unchanged)
    dat_sel = copy(dat)
    B0, M10, M20 = copy(df.B), copy(df.M1), copy(df.M2)
    eegfun.rereference!(dat_sel, :mastoid, eegfun.channels([:A]))
    mast_ref2 = (df.M1 .+ df.M2) ./ 2
    @test all(dat_sel.data.A .== df.A .- mast_ref2)
    @test all(dat_sel.data.B .== B0)
    @test all(dat_sel.data.M1 .== M10)
    @test all(dat_sel.data.M2 .== M20)

    # 7) Channel included in both reference and selection becomes zero
    dat_zero = copy(dat)
    eegfun.rereference!(dat_zero, [:A], eegfun.channels([:A]))
    @test all(dat_zero.data.A .== 0.0)

    # 8) ErpData mutate path
    erp_df = select(df, [:time, :A, :B, :M1, :M2])
    erp2 = eegfun.ErpData(copy(erp_df, copycols = true), layout, 1000, eegfun.AnalysisInfo(), 10)
    eegfun.rereference!(erp2, :mastoid, eegfun.channels([:A]))
    @test all(erp2.data.A .== erp_df.A .- ((erp_df.M1 .+ erp_df.M2) ./ 2))

    # 9) Missing reference channel should throw
    @test_throws Any eegfun.rereference!(copy(dat), [:Z], eegfun.channels([:A]))

    # 10) Reference info set for explicit vector
    dat_vec = copy(dat)
    eegfun.rereference!(dat_vec, [:M1, :M2], eegfun.channels([:A]))
    @test eegfun.reference(dat_vec.analysis_info) == Symbol("M1_M2")
end


# Helper function to create test ERP data
function create_test_erp_data(participant::Int, condition::Int, n_timepoints::Int = 100, n_channels::Int = 3)
    # Create time vector
    time = collect(range(-0.2, 0.8, length = n_timepoints))

    # Create channel data with some condition-specific patterns
    channel_data = Dict{Symbol,Any}()
    channel_data[:time] = time

    # Add metadata columns
    channel_data[:condition] = fill(condition, n_timepoints)
    channel_data[:condition_name] = fill("condition_$condition", n_timepoints)
    channel_data[:participant] = fill(participant, n_timepoints)

    # Add EEG channels with condition-specific patterns
    for (i, ch) in enumerate([:Fz, :Cz, :Pz][1:min(n_channels, 3)])
        # Create some signal with condition-specific amplitude
        signal = sin.(2π * 0.1 * time) .* (condition * 0.5) .+ randn(n_timepoints) * 0.1
        channel_data[ch] = signal
    end

    # Create DataFrame with columns in correct order
    df = DataFrame()
    df.time = channel_data[:time]
    df.condition = channel_data[:condition]
    df.condition_name = channel_data[:condition_name]
    df.participant = channel_data[:participant]
    for ch in [:Fz, :Cz, :Pz][1:min(n_channels, 3)]
        df[!, ch] = channel_data[ch]
    end

    # Create ErpData
    layout = eegfun.Layout(
        DataFrame(label = [:Fz, :Cz, :Pz][1:min(n_channels, 3)], inc = [0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0]),
        nothing,
        nothing,
    )
    return eegfun.ErpData(df, layout, 250.0, eegfun.AnalysisInfo(), 10)
end

# Helper function to create test epoch data
function create_test_epoch_data(
    participant::Int,
    condition::Int,
    n_epochs::Int = 5,
    n_timepoints::Int = 100,
    n_channels::Int = 3,
)
    epochs = []

    for epoch = 1:n_epochs
        # Create time vector
        time = collect(range(-0.2, 0.8, length = n_timepoints))

        # Create channel data
        channel_data = Dict{Symbol,Any}()
        channel_data[:time] = time
        channel_data[:condition] = fill(condition, n_timepoints)
        channel_data[:condition_name] = fill("condition_$condition", n_timepoints)
        channel_data[:participant] = fill(participant, n_timepoints)
        channel_data[:epoch] = fill(epoch, n_timepoints)

        # Add EEG channels
        for (i, ch) in enumerate([:Fz, :Cz, :Pz][1:min(n_channels, 3)])
            signal = sin.(2π * 0.1 * time) .* (condition * 0.5) .+ randn(n_timepoints) * 0.1
            channel_data[ch] = signal
        end

        # Create DataFrame with columns in correct order
        df = DataFrame()
        df.time = channel_data[:time]
        df.condition = channel_data[:condition]
        df.condition_name = channel_data[:condition_name]
        df.participant = channel_data[:participant]
        df.epoch = channel_data[:epoch]
        for ch in [:Fz, :Cz, :Pz][1:min(n_channels, 3)]
            df[!, ch] = channel_data[ch]
        end
        push!(epochs, df)
    end

    # Create EpochData
    layout = eegfun.Layout(
        DataFrame(label = [:Fz, :Cz, :Pz][1:min(n_channels, 3)], inc = [0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0]),
        nothing,
        nothing,
    )
    return eegfun.EpochData(epochs, layout, 250.0, eegfun.AnalysisInfo())
end

@testset "Batch Rereference" begin
    # Create temporary test directory
    test_dir = mktempdir()

    @testset "Basic rereferencing" begin
        # Create test ERP files
        for participant = 1:3
            erps = [create_test_erp_data(participant, 1), create_test_erp_data(participant, 2)]

            file_path = joinpath(test_dir, "$(participant)_erps_cleaned.jld2")
            save(file_path, "erps", erps)
        end

        output_dir = joinpath(test_dir, "rereferenced")

        # Test average reference
        result = eegfun.rereference(
            "erps_cleaned",
            input_dir = test_dir,
            reference_selection = :avg,
            output_dir = output_dir,
        )

        # Verify output files were created
        @test isdir(output_dir)
        output_files = readdir(output_dir)
        @test length(output_files) == 4  # One for each participant

        # Load and verify one output file
        output_file = joinpath(output_dir, "1_erps_cleaned.jld2")
        @test isfile(output_file)

        rereferenced_erps = load(output_file, "erps")
        @test length(rereferenced_erps) == 2  # Two conditions

        # Verify data structure is preserved
        for erp in rereferenced_erps
            @test erp isa eegfun.ErpData
            @test nrow(erp.data) == 100  # n_timepoints
            @test "time" in names(erp.data)
            @test "condition" in names(erp.data)
            @test "Fz" in names(erp.data)
            @test "Cz" in names(erp.data)
            @test "Pz" in names(erp.data)
        end
    end

    @testset "Different reference types" begin
        @testset "Average reference" begin
            output_dir = joinpath(test_dir, "rereferenced_avg")

            result = eegfun.rereference(
                "erps_cleaned",
                input_dir = test_dir,
                reference_selection = :avg,
                output_dir = output_dir,
            )

            @test isdir(output_dir)
            @test length(readdir(output_dir)) == 4
        end

        @testset "Mastoid reference" begin
            output_dir = joinpath(test_dir, "rereferenced_mastoid")

            result = eegfun.rereference(
                "erps_cleaned",
                input_dir = test_dir,
                reference_selection = :mastoid,
                output_dir = output_dir,
            )

            # Mastoid reference should fail because test data doesn't have M1, M2 channels
            @test result !== nothing
            @test result.success == 0  # No successful files
            @test result.errors > 0    # Should have errors
        end

        @testset "Single channel reference" begin
            output_dir = joinpath(test_dir, "rereferenced_cz")

            result = eegfun.rereference(
                "erps_cleaned",
                input_dir = test_dir,
                reference_selection = [:Cz],
                output_dir = output_dir,
            )

            @test isdir(output_dir)
            @test length(readdir(output_dir)) == 4
        end

        @testset "Multiple channel reference" begin
            output_dir = joinpath(test_dir, "rereferenced_multiple")

            result = eegfun.rereference(
                "erps_cleaned",
                input_dir = test_dir,
                reference_selection = [:Fz, :Pz],
                output_dir = output_dir,
            )

            @test isdir(output_dir)
            @test length(readdir(output_dir)) == 4
        end
    end

    @testset "Epoch data processing" begin
        # Create test epoch files
        for participant = 1:2
            epochs = [create_test_epoch_data(participant, 1, 3), create_test_epoch_data(participant, 2, 3)]

            file_path = joinpath(test_dir, "$(participant)_epochs_cleaned.jld2")
            save(file_path, "epochs", epochs)
        end

        output_dir = joinpath(test_dir, "rereferenced_epochs")

        result = eegfun.rereference(
            "epochs_cleaned",
            input_dir = test_dir,
            reference_selection = :avg,
            output_dir = output_dir,
        )

        @test isdir(output_dir)
        output_files = readdir(output_dir)
        @test length(output_files) >= 2  # At least one for each participant

        # Load and verify epoch data structure
        output_file = joinpath(output_dir, "1_epochs_cleaned.jld2")
        rereferenced_epochs = load(output_file, "epochs")
        @test length(rereferenced_epochs) == 2  # Two conditions

        for epoch_data in rereferenced_epochs
            @test epoch_data isa eegfun.EpochData
            @test length(epoch_data.data) == 3  # 3 epochs per condition
        end
    end

    @testset "Participant and condition filtering" begin
        output_dir = joinpath(test_dir, "rereferenced_filtered")

        # Test participant filtering
        result = eegfun.rereference(
            "erps_cleaned",
            input_dir = test_dir,
            participants = [1, 2],
            reference_selection = :avg,
            output_dir = output_dir,
        )

        @test isdir(output_dir)
        output_files = readdir(output_dir)
        @test length(output_files) == 3  # Only participants 1 and 2, plus empty file

        # Test condition filtering
        output_dir2 = joinpath(test_dir, "rereferenced_conditions")
        result = eegfun.rereference(
            "erps_cleaned",
            input_dir = test_dir,
            conditions = [1],
            reference_selection = :avg,
            output_dir = output_dir2,
        )

        @test isdir(output_dir2)
        output_files2 = readdir(output_dir2)
        @test length(output_files2) == 4  # All participants, but only condition 1

        # Verify only condition 1 was processed
        output_file = joinpath(output_dir2, "1_erps_cleaned.jld2")
        rereferenced_erps = load(output_file, "erps")
        @test length(rereferenced_erps) == 1  # Only condition 1
        @test rereferenced_erps[1].data.condition[1] == 1
    end

    @testset "Error handling" begin
        @testset "Invalid input directory" begin
            @test_throws Exception eegfun.rereference("erps_cleaned", input_dir = "/nonexistent/dir")
        end

        @testset "No matching files" begin
            output_dir = joinpath(test_dir, "rereferenced_none")

            result = eegfun.rereference("nonexistent_pattern", input_dir = test_dir, output_dir = output_dir)

            @test result === nothing
        end

        @testset "Files with no recognized data variable" begin
            # Create file with unrecognized variable and non-matching participant
            unrecognized_file = joinpath(test_dir, "999_unrecognized_erps_cleaned.jld2")
            save(unrecognized_file, "other_data", "test")

            output_dir = joinpath(test_dir, "rereferenced_unrecognized")

            result = eegfun.rereference(
                "erps_cleaned",
                input_dir = test_dir,
                participants = 999,  # This should match the file but fail processing
                output_dir = output_dir,
            )

            # Should return BatchResult with errors, not nothing
            @test result !== nothing
            @test result.success == 0
            @test result.errors > 0
        end
    end

    @testset "Edge cases" begin
        @testset "Empty data files" begin
            # Create file with empty data
            empty_file = joinpath(test_dir, "empty_erps_cleaned.jld2")
            save(empty_file, "erps", eegfun.ErpData[])

            output_dir = joinpath(test_dir, "rereferenced_empty")

            result = eegfun.rereference(
                "erps_cleaned",
                input_dir = test_dir,
                participants = 999,  # Non-existent participant
                output_dir = output_dir,
            )

            # Should return BatchResult since both empty_erps_cleaned.jld2 and 999_unrecognized_erps_cleaned.jld2 have no participant number
            @test result !== nothing
            @test result.success == 1  # Only the empty file gets processed successfully
            @test result.errors == 1   # One error for the unrecognized file
        end

        @testset "Invalid reference channels" begin
            output_dir = joinpath(test_dir, "rereferenced_invalid")

            # This should handle invalid channels gracefully in batch processing
            result = eegfun.rereference(
                "erps_cleaned",
                input_dir = test_dir,
                reference_selection = [:InvalidChannel],
                output_dir = output_dir,
            )

            # Should have errors due to invalid channels
            @test result !== nothing
            @test result.errors > 0
        end
    end

    @testset "Data integrity verification" begin
        output_dir = joinpath(test_dir, "rereferenced_integrity")

        result = eegfun.rereference(
            "erps_cleaned",
            input_dir = test_dir,
            participants = 1,
            conditions = [1],
            reference_selection = :avg,
            output_dir = output_dir,
        )

        # Load original and rereferenced data
        original_file = joinpath(test_dir, "1_erps_cleaned.jld2")
        original_erps = load(original_file, "erps")
        original_erp = original_erps[1]  # Condition 1

        rereferenced_file = joinpath(output_dir, "1_erps_cleaned.jld2")
        rereferenced_erps = load(rereferenced_file, "erps")
        rereferenced_erp = rereferenced_erps[1]  # Condition 1

        # Verify metadata is preserved
        @test rereferenced_erp.sample_rate == original_erp.sample_rate
        @test rereferenced_erp.n_epochs == original_erp.n_epochs

        # Verify layout is preserved (check channel labels)
        @test rereferenced_erp.layout.data.label == original_erp.layout.data.label

        # Verify analysis info reflects rereferencing
        @test rereferenced_erp.analysis_info.reference == :avg
        @test original_erp.analysis_info.reference == :none

        # Verify time column is unchanged
        @test all(rereferenced_erp.data.time .== original_erp.data.time)

        # Verify condition information is preserved
        @test all(rereferenced_erp.data.condition .== original_erp.data.condition)
        @test all(rereferenced_erp.data.condition_name .== original_erp.data.condition_name)
        @test all(rereferenced_erp.data.participant .== original_erp.data.participant)
    end

    @testset "Output directory handling" begin
        @testset "Custom output directory" begin
            custom_dir = joinpath(test_dir, "custom_rereferenced")

            result = eegfun.rereference(
                "erps_cleaned",
                input_dir = test_dir,
                reference_selection = :avg,
                output_dir = custom_dir,
            )

            @test isdir(custom_dir)
            @test length(readdir(custom_dir)) == 5
        end

        @testset "Auto-generated output directory" begin
            result = eegfun.rereference("erps_cleaned", input_dir = test_dir, reference_selection = :avg)

            # Should create directory with pattern-based name
            expected_dir = joinpath(test_dir, "rereferenced_erps_cleaned_avg")
            @test isdir(expected_dir)
        end
    end

    @testset "Logging and return values" begin
        output_dir = joinpath(test_dir, "rereferenced_logging")

        result = eegfun.rereference(
            "erps_cleaned",
            input_dir = test_dir,
            reference_selection = :avg,
            output_dir = output_dir,
        )

        # Check that log file was created
        log_file = joinpath(output_dir, "rereference.log")
        @test isfile(log_file)

        # Verify log content contains expected information
        log_content = read(log_file, String)
        @test occursin("Batch rereferencing started", log_content)
        @test occursin("Found 5 JLD2 files", log_content)
        @test occursin("Reference settings: avg", log_content)
        @test occursin("Batch operation complete", log_content)
    end

    @testset "Reference calculation verification" begin
        output_dir = joinpath(test_dir, "rereferenced_verification")

        result = eegfun.rereference(
            "erps_cleaned",
            input_dir = test_dir,
            participants = 1,
            conditions = [1],
            reference_selection = :avg,
            output_dir = output_dir,
        )

        # Load original and rereferenced data
        original_file = joinpath(test_dir, "1_erps_cleaned.jld2")
        original_erps = load(original_file, "erps")
        original_erp = original_erps[1]

        rereferenced_file = joinpath(output_dir, "1_erps_cleaned.jld2")
        rereferenced_erps = load(rereferenced_file, "erps")
        rereferenced_erp = rereferenced_erps[1]

        # Verify average reference calculation
        # Average reference should make the mean of all channels zero
        for ch in [:Fz, :Cz, :Pz]
            if hasproperty(rereferenced_erp.data, ch)
                # The mean across all channels at each time point should be approximately zero
                all_channels =
                    [rereferenced_erp.data[!, :Fz], rereferenced_erp.data[!, :Cz], rereferenced_erp.data[!, :Pz]]
                mean_across_channels = mean.(zip(all_channels...))
                @test all(abs.(mean_across_channels) .< 1e-10)
            end
        end
    end

    @testset "Different file patterns" begin
        # Test with different file patterns
        for pattern in ["erps_cleaned", "nonexistent", "original"]
            output_dir = joinpath(test_dir, "rereferenced_$pattern")

            result =
                eegfun.rereference(pattern, input_dir = test_dir, reference_selection = :avg, output_dir = output_dir)

            if pattern in ["erps_cleaned"]
                @test isdir(output_dir)
                @test length(readdir(output_dir)) > 0
            else
                # These patterns shouldn't match any files
                @test result === nothing
            end
        end
    end

    # Cleanup
    rm(test_dir, recursive = true)
end
