using Test
using JLD2
using DataFrames

@testset "Time-Frequency Batch Processing" begin

    @testset "Batch tf_morlet" begin
        mktempdir() do tmpdir
            # Create test epoch files (2 participants, 2 conditions each)
            for i = 1:2
                epochs = EegFun.create_test_epoch_data_vector(conditions = 1:2, n_epochs = 5, n_channels = 3, fs = 256, n = 256)
                jldsave(joinpath(tmpdir, "$(i)_epochs.jld2"); data = epochs)
            end

            output_dir = joinpath(tmpdir, "tf_morlet_output")

            # Run batch tf_morlet
            EegFun.tf_morlet("epochs"; input_dir = tmpdir, output_dir = output_dir, cycles = 7, frequencies = range(2, 20, length = 10))

            # Check output files exist
            @test isfile(joinpath(output_dir, "1_epochs.jld2"))
            @test isfile(joinpath(output_dir, "2_epochs.jld2"))

            # Load and verify structure
            result = load(joinpath(output_dir, "1_epochs.jld2"), "data")
            @test result isa Vector{EegFun.TimeFreqData}
            @test length(result) == 2  # 2 conditions
        end
    end

    @testset "Batch tf_stft" begin
        mktempdir() do tmpdir
            # Create test epoch files
            for i = 1:2
                epochs = EegFun.create_test_epoch_data_vector(conditions = 1:2, n_epochs = 5, n_channels = 3, fs = 256, n = 256)
                jldsave(joinpath(tmpdir, "$(i)_epochs.jld2"); data = epochs)
            end

            output_dir = joinpath(tmpdir, "tf_stft_output")

            # Run batch tf_stft with fixed window
            EegFun.tf_stft(
                "epochs";
                input_dir = tmpdir,
                output_dir = output_dir,
                window_length = 0.3,
                frequencies = range(2, 20, length = 10),
            )

            # Check output files exist
            @test isfile(joinpath(output_dir, "1_epochs.jld2"))
            @test isfile(joinpath(output_dir, "2_epochs.jld2"))

            # Load and verify structure
            result = load(joinpath(output_dir, "1_epochs.jld2"), "data")
            @test result isa Vector{EegFun.TimeFreqData}
            @test length(result) == 2
        end
    end

    @testset "Batch tf_multitaper" begin
        mktempdir() do tmpdir
            # Create test epoch files
            for i = 1:2
                epochs = EegFun.create_test_epoch_data_vector(conditions = 1:2, n_epochs = 5, n_channels = 3, fs = 256, n = 256)
                jldsave(joinpath(tmpdir, "$(i)_epochs.jld2"); data = epochs)
            end

            output_dir = joinpath(tmpdir, "tf_multitaper_output")

            # Run batch tf_multitaper
            EegFun.tf_multitaper("epochs"; input_dir = tmpdir, output_dir = output_dir, cycles = 5, frequencies = range(2, 20, length = 10))

            # Check output files exist
            @test isfile(joinpath(output_dir, "1_epochs.jld2"))
            @test isfile(joinpath(output_dir, "2_epochs.jld2"))

            # Load and verify structure
            result = load(joinpath(output_dir, "1_epochs.jld2"), "data")
            @test result isa Vector{EegFun.TimeFreqData}
            @test length(result) == 2
        end
    end

    @testset "Participant filtering" begin
        mktempdir() do tmpdir
            # Create files for 3 participants
            for i = 1:3
                epochs = EegFun.create_test_epoch_data_vector(conditions = 1:2, n_epochs = 5, n_channels = 3, fs = 256, n = 256)
                jldsave(joinpath(tmpdir, "$(i)_epochs.jld2"); data = epochs)
            end

            output_dir = joinpath(tmpdir, "tf_filtered")

            # Process only participants 1 and 3
            EegFun.tf_morlet(
                "epochs";
                input_dir = tmpdir,
                output_dir = output_dir,
                participant_selection = EegFun.participants([1, 3]),
                cycles = 7,
                frequencies = range(2, 20, length = 10),
            )

            # Only participants 1 and 3 should be processed
            @test isfile(joinpath(output_dir, "1_epochs.jld2"))
            @test !isfile(joinpath(output_dir, "2_epochs.jld2"))
            @test isfile(joinpath(output_dir, "3_epochs.jld2"))
        end
    end

    @testset "Condition selection" begin
        mktempdir() do tmpdir
            # Create file with 3 conditions
            epochs = EegFun.create_test_epoch_data_vector(conditions = 1:3, n_epochs = 5, n_channels = 3, fs = 256, n = 256)
            jldsave(joinpath(tmpdir, "1_epochs.jld2"); data = epochs)

            output_dir = joinpath(tmpdir, "tf_cond_select")

            # Select only conditions 1 and 3
            EegFun.tf_morlet(
                "epochs";
                input_dir = tmpdir,
                output_dir = output_dir,
                condition_selection = EegFun.conditions([1, 3]),
                cycles = 7,
                frequencies = range(2, 20, length = 10),
            )

            result = load(joinpath(output_dir, "1_epochs.jld2"), "data")
            @test length(result) == 2  # Only 2 conditions
        end
    end

    @testset "Direct tf_stft and tf_multitaper Execution" begin
        test_epochs = EegFun.create_test_epoch_data(n = 200, fs = 256, condition = 1, n_epochs = 4, n_channels = 3, seed = 42)

        # STFT fixed window
        res_stft_fixed = EegFun.tf_stft(test_epochs; window_length = 0.2, frequencies = 5:5:25, time_steps = 0.05)
        @test res_stft_fixed isa EegFun.TimeFreqData
        @test size(res_stft_fixed.data_power, 2) == 5  # time, freq, 3 channels
        @test all(isnan.(res_stft_fixed.data_phase[!, :Ch1]))  # phase NaN when return_phase=false

        # STFT adaptive window + return_trials + return_phase
        res_stft_trials = EegFun.tf_stft(
            test_epochs;
            cycles = 5,
            frequencies = 5:5:25,
            time_steps = 0.05,
            return_trials = true,
            return_phase = true,
            pad = :both,
        )
        @test res_stft_trials isa EegFun.TimeFreqEpochData
        @test length(res_stft_trials.data_power) == 4
        @test !any(isnan.(res_stft_trials.data_phase[1][!, :Ch1]))

        # Multitaper standard
        res_mt = EegFun.tf_multitaper(test_epochs; cycles = 5, frequencies = 5:5:25, time_steps = 0.05)
        @test res_mt isa EegFun.TimeFreqData
        @test size(res_mt.data_power, 2) == 5
        @test all(isnan.(res_mt.data_phase[!, :Ch1]))

        # Multitaper return_trials + return_phase
        res_mt_trials = EegFun.tf_multitaper(
            test_epochs;
            cycles = 5,
            frequencies = 5:5:25,
            time_steps = 0.05,
            return_trials = true,
            return_phase = true,
            pad = :pre,
        )
        @test res_mt_trials isa EegFun.TimeFreqEpochData
        @test length(res_mt_trials.data_power) == 4
        @test !all(isnan.(res_mt_trials.data_phase[1][!, :Ch1]))
        @test any(!isnan, res_mt_trials.data_phase[1][!, :Ch1])

        # Validation errors
        @test_throws EegFun.EegFunError EegFun.tf_stft(test_epochs)  # neither window_length nor cycles
        @test_throws EegFun.EegFunError EegFun.tf_stft(test_epochs; window_length = 0.2, cycles = 5)  # both
        @test_throws EegFun.EegFunError EegFun.tf_stft(test_epochs; window_length = -0.1)  # negative window
        @test_throws EegFun.EegFunError EegFun.tf_multitaper(test_epochs; cycles = -5)  # negative cycles
        @test_throws EegFun.EegFunError EegFun.tf_multitaper(test_epochs; cycles = 5, frequency_smoothing = -0.4)
    end

end

@testset "TF decompositions are thread-safe (all channels vs. one at a time)" begin
    # Many channels so several threads process channels concurrently with their own buffers;
    # each channel must match the result of decomposing that channel alone.
    ep = EegFun.create_test_epoch_data(n = 500, fs = 250, n_epochs = 20, n_channels = 16, seed = 11)
    freqs = 4:2:30
    methods = (
        (d; kw...) -> EegFun.tf_morlet(d; frequencies = freqs, cycles = 5, kw...),
        (d; kw...) -> EegFun.tf_stft(d; frequencies = freqs, window_length = 0.4, time_steps = 0.02, kw...),
        (d; kw...) -> EegFun.tf_multitaper(d; frequencies = freqs, cycles = 4, time_steps = 0.02, kw...),
    )
    for tf_fn in methods, return_trials in (false, true)
        full = tf_fn(ep; return_trials = return_trials, return_phase = true)
        for ch in EegFun.channel_labels(ep)
            single = tf_fn(ep; channel_selection = EegFun.channels(ch), return_trials = return_trials, return_phase = true)
            full_power = return_trials ? full.data_power[end][!, ch] : full.data_power[!, ch]
            single_power = return_trials ? single.data_power[end][!, ch] : single.data_power[!, ch]
            @test isequal(isnan.(full_power), isnan.(single_power))
            @test filter(!isnan, full_power) ≈ filter(!isnan, single_power) rtol = 1e-9
        end
    end
end
