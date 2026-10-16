using Test
using EegFun
using DataFrames
using Statistics

@testset "TF Utils" begin

    @testset "_tf_build_row_index" begin
        df = DataFrame(freq = [1.0, 2.0, 1.0, 2.0], time = [0.0, 0.0, 0.1, 0.1], Ch1 = [10.0, 20.0, 30.0, 40.0])
        idx = EegFun._tf_build_row_index(df)
        @test length(idx) == 4
        @test idx[(1.0, 0.0)] == 1
        @test idx[(2.0, 0.0)] == 2
        @test idx[(1.0, 0.1)] == 3
        @test idx[(2.0, 0.1)] == 4
    end

    @testset "_tf_df_to_matrix" begin
        # Create a simple TF DataFrame: 3 freqs × 4 times
        freqs = [5.0, 10.0, 15.0]
        times = [0.0, 0.1, 0.2, 0.3]
        n_f = length(freqs)
        n_t = length(times)

        freq_col = repeat(freqs, inner = n_t)
        time_col = repeat(times, outer = n_f)
        ch1_vals = collect(1.0:(n_f*n_t))

        df = DataFrame(freq = freq_col, time = time_col, Ch1 = ch1_vals)

        mat = EegFun._tf_df_to_matrix(df, :Ch1, freqs, times)

        @test size(mat) == (3, 4)

        # Verify all values roundtrip correctly
        row_index = EegFun._tf_build_row_index(df)
        for (fi, f) in enumerate(freqs)
            for (ti, t) in enumerate(times)
                row = row_index[(round(f, digits = 6), round(t, digits = 6))]
                @test mat[fi, ti] == df[row, :Ch1]
            end
        end

        # Test with subset of frequencies
        mat_sub = EegFun._tf_df_to_matrix(df, :Ch1, [5.0, 15.0], times)
        @test size(mat_sub) == (2, 4)
        @test mat_sub[1, :] == mat[1, :]  # 5 Hz row
        @test mat_sub[2, :] == mat[3, :]  # 15 Hz row

        # Test with missing freq/time → NaN
        mat_miss = EegFun._tf_df_to_matrix(df, :Ch1, [5.0, 99.0], times)
        @test size(mat_miss) == (2, 4)
        @test !any(isnan, mat_miss[1, :])  # 5 Hz exists
        @test all(isnan, mat_miss[2, :])   # 99 Hz doesn't exist
    end

    @testset "_tf_matrix_to_df!" begin
        freqs = [5.0, 10.0]
        times = [0.0, 0.1, 0.2]
        n_f = length(freqs)
        n_t = length(times)

        freq_col = repeat(freqs, inner = n_t)
        time_col = repeat(times, outer = n_f)

        df = DataFrame(freq = freq_col, time = time_col, Ch1 = zeros(n_f * n_t))

        mat = [1.0 2.0 3.0; 4.0 5.0 6.0]  # [2 freqs × 3 times]
        EegFun._tf_matrix_to_df!(df, :Ch1, mat, freqs, times)

        # Read back and verify
        mat_back = EegFun._tf_df_to_matrix(df, :Ch1, freqs, times)
        @test mat_back == mat
    end

    @testset "Roundtrip with float tolerance" begin
        # Test that slight float variations still match
        freqs = [1.0 / 3.0, 2.0 / 3.0]  # imprecise floats
        times = [0.0, 0.1]

        freq_col = repeat(freqs, inner = 2)
        time_col = repeat(times, outer = 2)
        df = DataFrame(freq = freq_col, time = time_col, Ch1 = [10.0, 20.0, 30.0, 40.0])

        mat = EegFun._tf_df_to_matrix(df, :Ch1, freqs, times)
        @test !any(isnan, mat)  # All values should be found
        @test mat[1, 1] == 10.0
        @test mat[1, 2] == 20.0
        @test mat[2, 1] == 30.0
        @test mat[2, 2] == 40.0
    end

    @testset "Row index reuse across channels" begin
        freqs = [1.0, 2.0]
        times = [0.0, 0.1]
        n = 4

        df = DataFrame(
            freq = repeat(freqs, inner = 2),
            time = repeat(times, outer = 2),
            Ch1 = [1.0, 2.0, 3.0, 4.0],
            Ch2 = [10.0, 20.0, 30.0, 40.0],
        )

        row_index = EegFun._tf_build_row_index(df)

        mat1 = EegFun._tf_df_to_matrix(df, :Ch1, freqs, times, row_index)
        mat2 = EegFun._tf_df_to_matrix(df, :Ch2, freqs, times, row_index)

        @test mat1[1, 1] == 1.0
        @test mat2[1, 1] == 10.0
        @test mat1[2, 2] == 4.0
        @test mat2[2, 2] == 40.0
    end

    @testset "_filter_edges! 3D layout is [trial, freq, time]" begin
        n_trials, n_f, n_t = 4, 3, 10
        time_indices = collect(1:n_t)
        n_samples = 10
        windows = [4, 6, 8]  # samples, per frequency

        power2d = ones(n_f, n_t)
        EegFun._filter_edges!(power2d, nothing, n_f, time_indices, windows, n_samples)
        @test any(isnan, power2d)
        @test !all(isnan, power2d)

        power3d = ones(n_trials, n_f, n_t)
        conv3d = ones(ComplexF64, n_trials, n_f, n_t)
        EegFun._filter_edges!(power3d, conv3d, n_f, time_indices, windows, n_samples)

        # Every trial must get exactly the same edge mask as the averaged (2D) case
        for trial = 1:n_trials
            @test isequal(power3d[trial, :, :], power2d)
            @test isequal(isnan.(conv3d[trial, :, :]), isnan.(power2d))
        end

        # Without eegconv
        power3d_only = ones(n_trials, n_f, n_t)
        EegFun._filter_edges!(power3d_only, nothing, n_f, time_indices, windows, n_samples)
        @test isequal(power3d_only, power3d)
    end

    @testset "return_trials with filter_edges (default) for all TF methods" begin
        fixture() = EegFun.create_test_epoch_data(n = 256, fs = 256, n_epochs = 6, n_channels = 3, seed = 1)
        freqs = 2:4:30
        runners = [
            "tf_morlet" => (; kw...) -> EegFun.tf_morlet(fixture(); frequencies = freqs, cycles = 5, kw...),
            "tf_stft" => (; kw...) -> EegFun.tf_stft(fixture(); frequencies = freqs, window_length = 0.3, kw...),
            "tf_multitaper" => (; kw...) -> EegFun.tf_multitaper(fixture(); frequencies = freqs, cycles = 5, kw...),
        ]

        for (name, run) in runners
            @testset "$name" begin
                avg = run()
                trials = run(return_trials = true, return_phase = true)

                @test trials.data_power isa Vector{DataFrame}
                @test length(trials.data_power) == 6

                for ch in (:Ch1, :Ch2, :Ch3)
                    avg_power = avg.data_power[!, ch]
                    @test any(isnan, avg_power)  # edges are actually exercised

                    for (df_pow, df_phase) in zip(trials.data_power, trials.data_phase)
                        # Same edge mask in every trial as in the averaged result
                        @test isequal(isnan.(df_pow[!, ch]), isnan.(avg_power))
                        @test isequal(isnan.(df_phase[!, ch]), isnan.(avg_power))
                    end

                    # Away from the edges, the mean over trials equals the averaged power
                    valid = .!isnan.(avg_power)
                    trial_mean = vec(mean(reduce(hcat, [df[!, ch] for df in trials.data_power]), dims = 2))
                    @test isapprox(trial_mean[valid], avg_power[valid]; rtol = 1e-8)
                end
            end
        end
    end
end
