using Test
using DataFrames
using EegFun
using Makie
using Statistics
using JLD2


@testset "filter" begin


    @testset "create_filter" begin
        fs = 1000.0
        # IIR low-pass
        fi_iir = EegFun.create_lowpass_filter(40.0, fs; order = 4, transition_width = 0.1)
        @test fi_iir isa EegFun.FilterInfo
        @test fi_iir.filter_type == "lp"
        @test fi_iir.filter_method == "iir"
        @test fi_iir.cutoff_freq == 40.0
        @test fi_iir.sample_rate == fs
        @test fi_iir.order == 4
        @test fi_iir.n_taps === nothing

        # FIR high-pass
        fi_fir = EegFun.create_highpass_filter(1.0, fs; filter_method = "fir", transition_width = 0.25)
        @test fi_fir isa EegFun.FilterInfo
        @test fi_fir.filter_type == "hp"
        @test fi_fir.filter_method == "fir"
        @test fi_fir.n_taps !== nothing
        @test fi_fir.n_taps % 2 == 1  # odd taps
        @test fi_fir.n_taps >= 101

        # FIR low-pass and tap sizing monotonicity
        fi_lp_wide = EegFun.create_lowpass_filter(40.0, fs; filter_method = "fir", transition_width = 0.2)
        fi_lp_narrow = EegFun.create_lowpass_filter(40.0, fs; filter_method = "fir", transition_width = 0.05)
        @test fi_lp_wide.n_taps !== nothing && fi_lp_narrow.n_taps !== nothing
        @test fi_lp_narrow.n_taps > fi_lp_wide.n_taps

    end

    @testset "filter_data! application and metadata" begin
        dat = create_test_data()
        dat_orig = copy(dat)

        # High-pass to remove DC; check mean is reduced towards ~0 for channel Ch1
        EegFun.highpass_filter!(dat, 1.0; order = 1, filter_method = "iir", channel_selection = EegFun.channels([:Ch1]))
        @test abs(mean(dat.data.Ch1)) < abs(mean(dat_orig.data.Ch1))
        # Only selected channel modified
        @test !all(dat.data.Ch1 .== dat_orig.data.Ch1)
        @test all(dat.data.Ch2 .== dat_orig.data.Ch2)
        # Analysis info updated
        @test dat.analysis_info.hp_filter == 1.0

        # Low-pass; update lp field and modify both channels when selecting both
        EegFun.lowpass_filter!(dat, 30.0; order = 3, filter_method = "iir", channel_selection = EegFun.channels([:Ch1, :Ch2]))
        @test dat.analysis_info.lp_filter == 30.0
        @test !all(dat.data.Ch2 .== dat_orig.data.Ch2)
    end

    @testset "non-mutating lowpass_filter" begin
        dat = create_test_data()
        dat_orig = copy(dat)
        dat2 = EegFun.lowpass_filter(dat, 30.0; order = 3)
        # Original unchanged
        @test all(dat.data.Ch1 .== dat_orig.data.Ch1)
        # Copy modified
        @test !all(dat2.data.Ch1 .== dat.data.Ch1)
    end

    @testset "no channels selected returns early" begin
        dat = create_test_data()
        dat_orig = copy(dat)
        # channel_selection picks none
        result = EegFun.highpass_filter!(dat, 1.0; channel_selection = EegFun.channels(Symbol[]))
        @test result === nothing
        # Data and analysis_info unchanged
        @test all(dat.data.Ch1 .== dat_orig.data.Ch1)
        @test dat.analysis_info.hp_filter == dat_orig.analysis_info.hp_filter
        @test dat.analysis_info.lp_filter == dat_orig.analysis_info.lp_filter
    end

    @testset "single-pass vs zero-phase" begin
        dat1 = create_test_data()
        dat2 = copy(dat1)
        # Single-pass introduces phase; zero-phase differs from single-pass
        EegFun.lowpass_filter!(dat1, 20.0; filter_func = "filt")
        EegFun.lowpass_filter!(dat2, 20.0; filter_func = "filtfilt")
        @test !all(dat1.data.Ch1 .== dat2.data.Ch1)
    end

    @testset "EpochData filtering" begin
        # Build two epochs from the same base
        base = create_test_data(; n = 1000, fs = 500)
        df1 = copy(base.data, copycols = true)
        df2 = copy(base.data, copycols = true)
        # Keep originals for comparison (distinct objects)
        df1o = copy(df1, copycols = true)
        df2o = copy(df2, copycols = true)
        # Add epoch identifiers to the versions that will be filtered
        df1.epoch = fill(1, nrow(df1))
        df2.epoch = fill(2, nrow(df2))
        ep = EegFun.EpochData(base.file, 1, "condition_1", [df1, df2], base.layout, base.sample_rate, EegFun.AnalysisInfo())
        EegFun.highpass_filter!(ep, 0.5)
        @test ep.analysis_info.hp_filter == 0.5
        @test !all(ep.data[1].Ch1 .== df1o.Ch1)
        @test !all(ep.data[2].Ch1 .== df2o.Ch1)
    end

    @testset "ErpData filtering" begin
        # Build ERP from base
        base = create_test_data(; n = 2000, fs = 1000)
        erp_df = select(base.data, [:time, :Ch1, :Ch2])
        erp = EegFun.ErpData(
            base.file,
            1,
            "condition_1",
            copy(erp_df, copycols = true),
            base.layout,
            base.sample_rate,
            EegFun.AnalysisInfo(),
            25,
        )
        erp_orig = copy(erp)
        EegFun.lowpass_filter!(erp, 30.0; order = 3)
        @test erp.analysis_info.lp_filter == 30.0
        @test !all(erp.data.Ch1 .== erp_orig.data.Ch1)
        # Non-mutating path
        erp2 = EegFun.highpass_filter(erp_orig, 0.5)
        @test erp_orig.analysis_info.hp_filter == 0.0  # unchanged
        @test erp2.analysis_info.hp_filter == 0.5
    end

    @testset "filter characteristics" begin
        fs = 1000.0
        fi = EegFun.create_lowpass_filter(40.0, fs; order = 4, transition_width = 0.1)
        chars = EegFun.get_filter_characteristics(fi; npoints = 256)
        @test chars.filter_type == "lp"
        @test isapprox(chars.transition_band, 4.0; atol = 1e-6)  # 40 Hz * 0.1 = 4.0 Hz
        @test any(abs.(chars.cutoff_freq_3db .- 40.0) .< 5.0)  # near cutoff
        @test chars.stopband_atten < -10  # should be attenuated
        # print helper should not error
        @test EegFun.print_filter_characteristics(fi; npoints = 128) === nothing
    end

    @testset "plot_filter_response (no display)" begin
        fs = 1000.0
        fi = EegFun.create_highpass_filter(1.0, fs)
        fig, axes = EegFun.plot_filter_response(fi; xscale = :linear, display_plot = false)
        @test fig isa Figure
        @test length(axes) == 3
        # log scale path
        fig2, axes2 = EegFun.plot_filter_response(fi; xscale = :log, display_plot = false)
        @test fig2 isa Figure
        @test length(axes2) == 3
    end

end


@testset "Batch Filter" begin

    # Create a temporary directory for test files
    test_dir = mktempdir()

    try
        # Create test data files
        @testset "Setup test files" begin
            for participant in [1, 2]
                erps = create_batch_test_erp_data(2)
                # Use filename format consistent with codebase (numeric participant ID)
                filename = joinpath(test_dir, "$(participant)_erps.jld2")
                jldsave(filename; data = erps)
                @test isfile(filename)
            end
        end

        @testset "Basic filtering" begin
            output_dir = joinpath(test_dir, "filtered_output")

            # Test low-pass filtering
            result = EegFun.lowpass_filter("erps", 30.0, input_dir = test_dir, output_dir = output_dir)

            @test result !== nothing
            @test result.success == 2
            @test result.errors == 0
            @test isdir(output_dir)

            # Check that filtered files exist
            @test isfile(joinpath(output_dir, "1_erps.jld2"))
            @test isfile(joinpath(output_dir, "2_erps.jld2"))

            # Load and verify filtered data
            filtered_data = load(joinpath(output_dir, "1_erps.jld2"), "data")
            @test length(filtered_data) == 2  # 2 conditions
            @test filtered_data[1] isa EegFun.ErpData
            @test hasproperty(filtered_data[1].data, :Ch1)
            @test hasproperty(filtered_data[1].data, :Ch2)

            # Verify high frequencies are attenuated (not a rigorous test, just sanity check)
            @test std(filtered_data[1].data.Ch1) < 2.0  # Should be reduced
        end

        @testset "Filter specific participants" begin
            output_dir = joinpath(test_dir, "filtered_participant")

            result = EegFun.lowpass_filter(
                "erps",
                30.0,
                input_dir = test_dir,
                output_dir = output_dir,
                participant_selection = EegFun.participants(1),
            )

            @test result.success == 1
            @test result.errors == 0
            @test isfile(joinpath(output_dir, "1_erps.jld2"))
            @test !isfile(joinpath(output_dir, "2_erps.jld2"))
        end

        @testset "Filter specific conditions" begin
            output_dir = joinpath(test_dir, "filtered_condition")

            result = EegFun.lowpass_filter(
                "erps",
                30.0,
                input_dir = test_dir,
                output_dir = output_dir,
                condition_selection = EegFun.conditions(1),
            )

            @test result.success == 2

            # Load and verify only one condition
            filtered_data = load(joinpath(output_dir, "1_erps.jld2"), "data")
            @test length(filtered_data) == 1
            @test filtered_data[1].data[1, :condition] == 1
        end

        @testset "High-pass filter" begin
            output_dir = joinpath(test_dir, "filtered_hp")

            result = EegFun.highpass_filter("erps", 1.0, input_dir = test_dir, output_dir = output_dir)

            @test result.success == 2
            @test result.errors == 0

            # Load and verify
            filtered_data = load(joinpath(output_dir, "1_erps.jld2"), "data")
            @test filtered_data[1] isa EegFun.ErpData
        end

        @testset "Error handling" begin
            # Non-existent directory
            @test_throws Exception EegFun.lowpass_filter("erps", 30.0, input_dir = "/nonexistent/path")

            # Invalid cutoff frequency
            @test_throws Exception EegFun.lowpass_filter("erps", -10.0, input_dir = test_dir)
        end

        @testset "No matching files" begin
            output_dir = joinpath(test_dir, "filtered_nomatch")

            # Pattern that won't match any files
            result = EegFun.lowpass_filter("nonexistent_pattern", 30.0, input_dir = test_dir, output_dir = output_dir)

            @test result === nothing  # Function returns nothing when no files found
        end

        @testset "Logging" begin
            output_dir = joinpath(test_dir, "filtered_with_log")

            # Use lowpass_filter explicitly
            result = EegFun.lowpass_filter("erps", 30.0, input_dir = test_dir, output_dir = output_dir)

            # Check log file exists
            # Check log file exists (it should be filter_lp.log for lowpass)
            log_file = joinpath(output_dir, "filter_lp.log")
            @test isfile(log_file)

            # Verify log contains expected information
            log_contents = read(log_file, String)
            @test contains(log_contents, "filtering started")
            # The log message now says "lowpass" or "highpass" instead of "filter_type=\"lp\""
            @test contains(log_contents, "lowpass")
            @test contains(log_contents, "cutoff: 30.0 Hz")
        end

        @testset "Multiple participants filter" begin
            output_dir = joinpath(test_dir, "filtered_multi_participants")

            # Filter both participants
            result = EegFun.lowpass_filter(
                "erps",
                30.0,
                input_dir = test_dir,
                output_dir = output_dir,
                participant_selection = EegFun.participants([1, 2]),
            )

            @test result.success == 2
            @test result.errors == 0
            @test isfile(joinpath(output_dir, "1_erps.jld2"))
            @test isfile(joinpath(output_dir, "2_erps.jld2"))
        end

        @testset "Multiple conditions filter" begin
            output_dir = joinpath(test_dir, "filtered_multi_conditions")

            # Filter both conditions
            result = EegFun.lowpass_filter(
                "erps",
                30.0,
                input_dir = test_dir,
                output_dir = output_dir,
                condition_selection = EegFun.conditions([1, 2]),
            )

            @test result.success == 2

            # Load and verify both conditions are present
            filtered_data = load(joinpath(output_dir, "1_erps.jld2"), "data")
            @test length(filtered_data) == 2
            @test filtered_data[1].data[1, :condition] == 1
            @test filtered_data[2].data[1, :condition] == 2
        end

        @testset "EpochData support" begin
            # Create EpochData test files
            epochs_dir = joinpath(test_dir, "epochs_test")
            mkpath(epochs_dir)

            # Use generic create_test_epoch_data from test_utils.jl
            # create_test_epoch_data(participant, condition, n_timepoints, n_channels)

            # Save epoch data - create a vector of EpochData for batch processing
            epochs = [create_test_epoch_data(conditions = 1), create_test_epoch_data(conditions = 1)]
            jldsave(joinpath(epochs_dir, "1_epochs.jld2"); data = epochs)

            # Filter epoch data
            output_dir = joinpath(test_dir, "filtered_epochs")
            result = EegFun.lowpass_filter("epochs", 30.0, input_dir = epochs_dir, output_dir = output_dir)

            @test result.success == 1
            @test result.errors == 0
            @test isfile(joinpath(output_dir, "1_epochs.jld2"))

            # Load and verify
            filtered_epochs = load(joinpath(output_dir, "1_epochs.jld2"), "data")
            @test filtered_epochs[1] isa EegFun.EpochData
            @test length(filtered_epochs) == 2  # 2 conditions
        end

        @testset "Existing output directory" begin
            output_dir = joinpath(test_dir, "existing_output")
            mkpath(output_dir)

            # Create a dummy file in the output directory
            touch(joinpath(output_dir, "dummy.txt"))
            @test isfile(joinpath(output_dir, "dummy.txt"))

            # Run filter - should work fine with existing directory
            result = EegFun.lowpass_filter("erps", 30.0, input_dir = test_dir, output_dir = output_dir)

            @test result.success == 2
            @test isfile(joinpath(output_dir, "dummy.txt"))  # Original file preserved
            @test isfile(joinpath(output_dir, "1_erps.jld2"))
        end

        @testset "Partial failures" begin
            partial_dir = joinpath(test_dir, "partial_test")
            mkpath(partial_dir)

            # Create one valid file
            erps = create_batch_test_erp_data(2)
            jldsave(joinpath(partial_dir, "1_erps.jld2"); data = erps)

            # Create one malformed file (invalid data type - String instead of Vector{ErpData})
            jldsave(joinpath(partial_dir, "2_erps.jld2"); data = "invalid_data")

            output_dir = joinpath(test_dir, "filtered_partial")
            result = EegFun.lowpass_filter("erps", 30.0, input_dir = partial_dir, output_dir = output_dir)

            @test result.success == 1
            @test result.errors == 1
            @test isfile(joinpath(output_dir, "1_erps.jld2"))
            @test !isfile(joinpath(output_dir, "2_erps.jld2"))  # Failed file not saved
        end

        @testset "Condition out of range" begin
            output_dir = joinpath(test_dir, "filtered_invalid_condition")

            # Request condition 5 when only 2 exist
            # With predicate-based selection, this results in empty selection but successful processing
            result = EegFun.lowpass_filter(
                "erps",
                30.0,
                input_dir = test_dir,
                output_dir = output_dir,
                condition_selection = EegFun.conditions(5),
            )

            # Files are processed successfully but with empty condition selection
            @test result.success == 2
            @test result.errors == 0
        end

        @testset "Empty pattern match" begin
            empty_dir = joinpath(test_dir, "empty_match")
            mkpath(empty_dir)

            # Directory exists but has no JLD2 files
            result = EegFun.lowpass_filter("erps", 30.0, input_dir = empty_dir)

            @test result === nothing  # No files to process
        end

        @testset "Return value structure" begin
            output_dir = joinpath(test_dir, "filtered_return_check")

            result = EegFun.lowpass_filter("erps", 30.0, input_dir = test_dir, output_dir = output_dir)

            # Check result structure
            @test hasfield(typeof(result), :success)
            @test hasfield(typeof(result), :errors)
            @test result.success isa Integer
            @test result.errors isa Integer
            @test result.success >= 0
            @test result.errors >= 0
            @test result.success + result.errors == 2  # Total files processed
        end

        @testset "Data integrity - frequency attenuation" begin
            output_dir = joinpath(test_dir, "filtered_integrity")

            # Get original data statistics
            original_data = load(joinpath(test_dir, "1_erps.jld2"), "data")
            original_signal = original_data[1].data.Ch1

            # Apply low-pass filter
            EegFun.lowpass_filter("erps", 30.0, input_dir = test_dir, output_dir = output_dir)

            # Load filtered data
            filtered_data = load(joinpath(output_dir, "1_erps.jld2"), "data")
            filtered_signal = filtered_data[1].data.Ch1

            # Check that high-frequency noise is reduced (lower std deviation)
            @test std(filtered_signal) < std(original_signal)

            # Check that data length is preserved
            @test length(filtered_signal) == length(original_signal)

            # Check that time vector is preserved
            @test filtered_data[1].data.time == original_data[1].data.time
        end

        @testset "Combined filters" begin
            output_dir = joinpath(test_dir, "filtered_combined")

            # Filter specific participant AND condition
            result = EegFun.lowpass_filter(
                "erps",
                30.0,
                input_dir = test_dir,
                output_dir = output_dir,
                participant_selection = EegFun.participants(1),
                condition_selection = EegFun.conditions(1),
            )

            @test result.success == 1

            # Load and verify
            filtered_data = load(joinpath(output_dir, "1_erps.jld2"), "data")
            @test length(filtered_data) == 1  # Only one condition
            @test filtered_data[1].data[1, :condition] == 1
            @test !isfile(joinpath(output_dir, "2_erps.jld2"))  # Participant 2 not processed
        end

    finally
        # Clean up
        rm(test_dir, recursive = true, force = true)
    end
end
