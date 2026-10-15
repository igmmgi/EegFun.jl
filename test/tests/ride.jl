using Test
using EegFun
using Statistics
using DataFrames

@testset "RIDE (Residue Iteration DEcomposition)" begin

    # 1. Shifting Utilities
    @testset "Signal and Matrix Shifting" begin
        v = [1.0, 2.0, 3.0, 4.0, 5.0]

        # In-place shift_signal!
        dest = zeros(5)
        EegFun.shift_signal!(dest, v, 2)
        @test dest == [0.0, 0.0, 1.0, 2.0, 3.0]

        EegFun.shift_signal!(dest, v, -2)
        @test dest == [3.0, 4.0, 5.0, 0.0, 0.0]

        EegFun.shift_signal!(dest, v, 0)
        @test dest == v

        EegFun.shift_signal!(dest, v, 5)
        @test dest == zeros(5)

        EegFun.shift_signal!(dest, v, -10)
        @test dest == zeros(5)

        # Out-of-place shift_signal
        @test EegFun.shift_signal(v, 1) == [0.0, 1.0, 2.0, 3.0, 4.0]
        @test EegFun.shift_signal(v, -1) == [2.0, 3.0, 4.0, 5.0, 0.0]

        # Matrix shifting
        M = [1.0 10.0; 2.0 20.0; 3.0 30.0; 4.0 40.0]
        M_shift = EegFun.shift_matrix(M, 1)
        @test M_shift[:, 1] == [0.0, 1.0, 2.0, 3.0]
        @test M_shift[:, 2] == [0.0, 10.0, 20.0, 30.0]

        M_shift_neg = EegFun.shift_matrix(M, -1)
        @test M_shift_neg[:, 1] == [2.0, 3.0, 4.0, 0.0]
        @test M_shift_neg[:, 2] == [20.0, 30.0, 40.0, 0.0]
    end

    # 2. Woody Filter & Peak Finding
    @testset "Woody Filter & Peak Detection" begin
        # Peak finding
        sig = [0.0, 1.0, 0.5, 2.0, 0.0]
        p_idx, found = EegFun.find_nearest_peak(sig, 2)
        @test found && p_idx == 2
        p_idx2, found2 = EegFun.find_nearest_peak(sig, 5)
        @test found2 && p_idx2 == 4

        # Flat signal has no local peaks
        flat = [1.0, 1.0, 1.0, 1.0]
        _, found_flat = EegFun.find_nearest_peak(flat, 2)
        @test !found_flat

        # Woody filter on synthetic single-channel tensor
        n_samples = 60
        n_trials = 8
        n_channels = 1
        data = zeros(Float64, n_samples, n_channels, n_trials)
        wave = sin.(range(0, 2π, length = 20))
        for t = 1:n_trials
            # Shift wave by t-4
            s_idx = 20 + (t - 4)
            data[s_idx:(s_idx+19), 1, t] .= wave
        end

        lats = EegFun.woody_filter(data, 10)
        @test length(lats) == n_trials
        # Latencies should correlate with trial shift
        @test cor(Float64.(lats), Float64.(1:n_trials)) > 0.9
    end

    # 3. Iteration and Decomposition
    @testset "RIDE Iteration & Decomposition" begin
        n_samples = 80
        n_channels = 2
        n_trials = 6

        data = zeros(Float64, n_samples, n_channels, n_trials)
        s_wave = sin.(range(0, 2π, length = 25))
        for t = 1:n_trials
            data[10:34, :, t] .+= s_wave
        end

        # Pre-defined S component
        comp_s = EegFun.RideComponent(:S, zeros(Int, n_trials))
        @test comp_s.name == :S
        @test comp_s.latencies == zeros(Int, n_trials)
        @test comp_s.window_radius == 0

        # Estimated C component
        comp_c = EegFun.RideComponent(:C, 10)
        @test comp_c.name == :C
        @test comp_c.latencies === nothing
        @test comp_c.window_radius == 10

        # Run decomposition
        res = EegFun.ride_decompose(data, [comp_s]; outer_iters = 2, inner_iters = 10, method = :mean)
        @test res isa EegFun.RideResult
        @test res.components == [:S]
        @test size(res.templates) == (n_samples, n_channels, 1)
        @test size(res.residue) == size(data)

        # Median method
        res_med = EegFun.ride_decompose(data, [comp_s]; outer_iters = 2, inner_iters = 10, method = :median)
        @test res_med isa EegFun.RideResult

        # Error handling
        @test_throws EegFun.EegFunError EegFun.ride_decompose(data, [EegFun.RideComponent(:S, zeros(Int, n_trials + 2))])
    end

    # 4. extract_ride (EpochData interface)
    @testset "extract_ride" begin
        n_samples = 50
        n_trials = 5
        n_channels = 3
        epochs = EegFun.create_test_epoch_data(n = n_samples, fs = 100, n_epochs = n_trials, n_channels = n_channels, seed = 1)

        comps = [EegFun.RideComponent(:S, zeros(Int, n_trials))]
        ride_dict = EegFun.extract_ride(epochs, comps; outer_iters = 1, inner_iters = 5)

        @test haskey(ride_dict, :S)
        @test ride_dict[:S] isa EegFun.ErpData
        @test length(EegFun.channel_labels(ride_dict[:S])) == n_channels
        @test nrow(ride_dict[:S].data) == n_samples
        @test occursin("RIDE_S", ride_dict[:S].condition_name)

        # Empty epochs
        empty_epochs = copy(epochs)
        empty_epochs.data = DataFrame[]
        @test_throws EegFun.EegFunError EegFun.extract_ride(empty_epochs, comps)
    end
end
