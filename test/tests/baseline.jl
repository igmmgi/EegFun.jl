using Test
using DataFrames
using Statistics
using EegFun

@testset "baseline" begin

    # Baseline over first sample for Ch1, Ch2, and Ch3
    dat = EegFun.create_test_continuous_data(n = 6)
    EegFun.baseline!(dat, (0.0, 0.0))

    @test isapprox(mean(dat.data.Ch1[1]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch2[1]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch3[1]), 0.0; atol = 1e-9)

    dat = EegFun.create_test_continuous_data(n = 6)
    dat = EegFun.baseline(dat, (0.0, 0.0))

    @test isapprox(mean(dat.data.Ch1[1]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch2[1]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch3[1]), 0.0; atol = 1e-9)

    # Baseline over entire range for all channels
    dat = EegFun.create_test_continuous_data(n = 6)
    EegFun.baseline!(dat)

    @test isapprox(mean(dat.data.Ch1), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch2), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch3), 0.0; atol = 1e-9)

    dat = EegFun.create_test_continuous_data(n = 6)
    dat = EegFun.baseline(dat)

    @test isapprox(mean(dat.data.Ch1), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch2), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch3), 0.0; atol = 1e-9)


    # EpochData: each epoch baselined independently
    dat = EegFun.create_test_epoch_data(n = 3)
    EegFun.baseline!(dat, (0.0, 0.0))

    @test isapprox(mean(dat.data[1].Ch1[1]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data[2].Ch2[1]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data[3].Ch3[1]), 0.0; atol = 1e-9)

    dat = EegFun.create_test_epoch_data(n = 3)
    EegFun.baseline!(dat, (0.0, 0.002))

    @test isapprox(mean(dat.data[1].Ch1[1:3]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data[2].Ch2[1:3]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data[3].Ch3[1:3]), 0.0; atol = 1e-9)


    # Tuple converted to indices correctly
    dat = EegFun.create_test_continuous_data(n = 6)
    EegFun.baseline!(dat, (0.003, 0.003))

    @test isapprox(mean(dat.data.Ch1[4]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch2[4]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch3[4]), 0.0; atol = 1e-9)

    dat = EegFun.create_test_continuous_data(n = 6)
    dat = EegFun.baseline(dat, (0.003, 0.003))

    @test isapprox(mean(dat.data.Ch1[4]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch2[4]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch3[4]), 0.0; atol = 1e-9)


    dat = EegFun.create_test_epoch_data(n = 6)
    EegFun.baseline!(dat, (0.003, 0.003))

    @test isapprox(mean(dat.data[1].Ch1[4]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data[2].Ch2[4]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data[3].Ch3[4]), 0.0; atol = 1e-9)

    dat = EegFun.create_test_epoch_data(n = 6)
    dat = EegFun.baseline(dat, (0.003, 0.003))

    @test isapprox(mean(dat.data[1].Ch1[4]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data[2].Ch2[4]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data[3].Ch3[4]), 0.0; atol = 1e-9)

    # Test channel selection
    dat = EegFun.create_test_continuous_data(n = 10)
    original_ch1 = copy(dat.data.Ch1)
    original_ch2 = copy(dat.data.Ch2)
    original_ch3 = copy(dat.data.Ch3)

    # Baseline only Ch1 and Ch2
    EegFun.baseline!(dat, (0.0, 0.004), channel_selection = EegFun.channels([:Ch1, :Ch2]))

    # Ch1 and Ch2 should be baselined
    @test isapprox(mean(dat.data.Ch1[1:5]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch2[1:5]), 0.0; atol = 1e-9)

    # Ch3 should not be baselined (should be different from original)
    baseline_mean_ch1 = mean(original_ch1[1:5])
    baseline_mean_ch2 = mean(original_ch2[1:5])
    @test isapprox(dat.data.Ch1, original_ch1 .- baseline_mean_ch1; atol = 1e-9)
    @test isapprox(dat.data.Ch2, original_ch2 .- baseline_mean_ch2; atol = 1e-9)
    @test dat.data.Ch3 == original_ch3  # Ch3 unchanged

    # Test with non-mutating version and channel selection
    dat = EegFun.create_test_continuous_data(n = 10)
    original_ch1 = copy(dat.data.Ch1)
    baseline_mean_ch1 = mean(original_ch1[1:5])

    dat_baselined = EegFun.baseline(dat, (0.0, 0.004), channel_selection = EegFun.channels([:Ch1]))

    @test isapprox(mean(dat_baselined.data.Ch1[1:5]), 0.0; atol = 1e-9)
    @test isapprox(dat_baselined.data.Ch1, original_ch1 .- baseline_mean_ch1; atol = 1e-9)
    @test dat_baselined !== dat  # Should be different object

    # Test empty channel selection (should warn and return early)
    dat = EegFun.create_test_continuous_data(n = 10)
    original_ch1 = copy(dat.data.Ch1)

    EegFun.baseline!(dat, (0.0, 0.004), channel_selection = EegFun.channels(Symbol[]))

    # Data should be unchanged
    @test dat.data.Ch1 == original_ch1

    # Test with EpochData and channel selection
    epochs = EegFun.create_test_epoch_data(n = 10, n_epochs = 3)
    original_epoch1_ch1 = copy(epochs.data[1].Ch1)
    baseline_mean_epoch1 = mean(original_epoch1_ch1[1:5])

    EegFun.baseline!(epochs, (0.0, 0.004), channel_selection = EegFun.channels([:Ch1]))

    @test isapprox(mean(epochs.data[1].Ch1[1:5]), 0.0; atol = 1e-9)
    @test isapprox(epochs.data[1].Ch1, original_epoch1_ch1 .- baseline_mean_epoch1; atol = 1e-9)

    # Test with ErpData
    erp = EegFun.create_test_erp_data(participant = 1, condition = 1, n_channels = 3)
    original_erp_ch1 = copy(erp.data.Ch1)
    baseline_mean_erp = mean(original_erp_ch1[1:10])

    EegFun.baseline!(erp, (-0.5, -0.491))

    @test isapprox(mean(erp.data.Ch1[1:10]), 0.0; atol = 1e-9)
    @test isapprox(erp.data.Ch1, original_erp_ch1 .- baseline_mean_erp; atol = 1e-9)

    # Test with ErpData and channel selection
    erp = EegFun.create_test_erp_data(participant = 1, condition = 1, n_channels = 3)
    original_erp_ch1 = copy(erp.data.Ch1)
    original_erp_ch2 = copy(erp.data.Ch2)
    baseline_mean_erp = mean(original_erp_ch1[1:10])

    EegFun.baseline!(erp, (-0.5, -0.491), channel_selection = EegFun.channels([:Ch1]))

    @test isapprox(mean(erp.data.Ch1[1:10]), 0.0; atol = 1e-9)
    @test isapprox(erp.data.Ch1, original_erp_ch1 .- baseline_mean_erp; atol = 1e-9)
    @test erp.data.Ch2 == original_erp_ch2  # Ch2 unchanged

    # Test baseline! without interval (uses entire range)
    dat = EegFun.create_test_continuous_data(n = 20)
    original_ch1 = copy(dat.data.Ch1)
    baseline_mean_ch1 = mean(original_ch1)

    EegFun.baseline!(dat)

    @test isapprox(mean(dat.data.Ch1), 0.0; atol = 1e-9)
    @test isapprox(dat.data.Ch1, original_ch1 .- baseline_mean_ch1; atol = 1e-9)

    # Test baseline (non-mutating) without interval
    dat = EegFun.create_test_continuous_data(n = 20)
    original_ch1 = copy(dat.data.Ch1)
    baseline_mean_ch1 = mean(original_ch1)

    dat_baselined = EegFun.baseline(dat)

    @test isapprox(mean(dat_baselined.data.Ch1), 0.0; atol = 1e-9)
    @test isapprox(dat_baselined.data.Ch1, original_ch1 .- baseline_mean_ch1; atol = 1e-9)
    @test dat_baselined !== dat

    # Test with larger baseline interval
    dat = EegFun.create_test_continuous_data(n = 100)
    original_ch1 = copy(dat.data.Ch1)
    baseline_mean_ch1 = mean(original_ch1[10:50])

    EegFun.baseline!(dat, (0.009, 0.049))

    @test isapprox(mean(dat.data.Ch1[10:50]), 0.0; atol = 1e-9)
    @test isapprox(dat.data.Ch1, original_ch1 .- baseline_mean_ch1; atol = 1e-9)

    # Test tuple with larger range
    dat = EegFun.create_test_continuous_data(n = 100, fs = 1000)
    original_ch1 = copy(dat.data.Ch1)
    # Find indices for time range 0.01 to 0.05 seconds
    time_idx_start = findfirst(x -> x >= 0.01, dat.data.time)
    time_idx_stop = findlast(x -> x <= 0.05, dat.data.time)
    baseline_mean_ch1 = mean(original_ch1[time_idx_start:time_idx_stop])

    EegFun.baseline!(dat, (0.01, 0.05))

    @test isapprox(mean(dat.data.Ch1[time_idx_start:time_idx_stop]), 0.0; atol = 1e-9)
    @test isapprox(dat.data.Ch1, original_ch1 .- baseline_mean_ch1; atol = 1e-9)

    # Test that baseline correction is applied to entire signal, not just baseline interval
    dat = EegFun.create_test_continuous_data(n = 100)
    original_ch1 = copy(dat.data.Ch1)
    baseline_mean_ch1 = mean(original_ch1[1:10])
    original_mean_ch1 = mean(original_ch1)

    EegFun.baseline!(dat, (0.0, 0.009))

    # Mean of baseline interval should be 0
    @test isapprox(mean(dat.data.Ch1[1:10]), 0.0; atol = 1e-9)
    # Entire signal should be shifted by baseline mean
    @test isapprox(dat.data.Ch1, original_ch1 .- baseline_mean_ch1; atol = 1e-9)
    # Mean of entire signal should be original_mean - baseline_mean
    @test isapprox(mean(dat.data.Ch1), original_mean_ch1 - baseline_mean_ch1; atol = 1e-9)

    # Test with EpochData - each epoch baselined independently
    epochs = EegFun.create_test_epoch_data(n = 10, n_epochs = 3)
    original_epoch_ch1_1 = copy(epochs.data[1].Ch1)
    original_epoch_ch1_2 = copy(epochs.data[2].Ch1)
    original_epoch_ch1_3 = copy(epochs.data[3].Ch1)

    baseline_mean_epoch_1 = mean(original_epoch_ch1_1[1:5])
    baseline_mean_epoch_2 = mean(original_epoch_ch1_2[1:5])
    baseline_mean_epoch_3 = mean(original_epoch_ch1_3[1:5])

    EegFun.baseline!(epochs, (0.0, 0.004))

    # Each epoch should be baselined independently
    @test isapprox(mean(epochs.data[1].Ch1[1:5]), 0.0; atol = 1e-9)
    @test isapprox(mean(epochs.data[2].Ch1[1:5]), 0.0; atol = 1e-9)
    @test isapprox(mean(epochs.data[3].Ch1[1:5]), 0.0; atol = 1e-9)

    @test isapprox(epochs.data[1].Ch1, original_epoch_ch1_1 .- baseline_mean_epoch_1; atol = 1e-9)
    @test isapprox(epochs.data[2].Ch1, original_epoch_ch1_2 .- baseline_mean_epoch_2; atol = 1e-9)
    @test isapprox(epochs.data[3].Ch1, original_epoch_ch1_3 .- baseline_mean_epoch_3; atol = 1e-9)

    # Test error handling - invalid interval (start > stop)
    # validate_baseline_interval throws ErrorException for invalid intervals
    dat = EegFun.create_test_continuous_data(n = 10)
    @test_throws Exception EegFun.baseline!(dat, (0.005, 0.001))

    # Test error handling - invalid interval (negative start time)
    dat = EegFun.create_test_continuous_data(n = 10)
    @test_throws Exception EegFun.baseline!(dat, (-1.0, 0.005))

    # Test error handling - invalid interval (stop time past data range)
    dat = EegFun.create_test_continuous_data(n = 10)
    @test_throws Exception EegFun.baseline!(dat, (0.0, 100.0))

    # Test error handling - invalid tuple (completely outside time range)
    # When tuple is outside range, _find_idx_start_end returns nothing
    dat = EegFun.create_test_continuous_data(n = 10, fs = 1000)
    @test_throws Exception EegFun.baseline!(dat, (100.0, 200.0))

    # Test with single sample baseline interval
    dat = EegFun.create_test_continuous_data(n = 20)
    original_ch1 = copy(dat.data.Ch1)
    baseline_value = original_ch1[5]  # Single sample value

    EegFun.baseline!(dat, (0.004, 0.004))

    @test isapprox(dat.data.Ch1[5], 0.0; atol = 1e-9)
    @test isapprox(dat.data.Ch1, original_ch1 .- baseline_value; atol = 1e-9)

    # Test that metadata columns are not affected
    dat = EegFun.create_test_continuous_data(n = 20)
    original_time = copy(dat.data.time)
    original_sample = copy(dat.data.sample)
    original_triggers = copy(dat.data.triggers)

    EegFun.baseline!(dat, (0.0, 0.009))

    @test dat.data.time == original_time
    @test dat.data.sample == original_sample
    @test dat.data.triggers == original_triggers

    # Test baseline! without interval on EpochData (uses entire time range)
    # This test catches the bug where dat.data.time was accessed directly on Vector{DataFrame}
    epochs_no_interval = EegFun.create_test_epoch_data(n = 20, n_epochs = 3)
    original_epoch1_ch1_no_int = copy(epochs_no_interval.data[1].Ch1)
    original_epoch2_ch1_no_int = copy(epochs_no_interval.data[2].Ch1)
    baseline_mean_epoch1_no_int = mean(original_epoch1_ch1_no_int)
    baseline_mean_epoch2_no_int = mean(original_epoch2_ch1_no_int)

    EegFun.baseline!(epochs_no_interval)

    # Each epoch should be baselined to its own mean
    @test isapprox(mean(epochs_no_interval.data[1].Ch1), 0.0; atol = 1e-9)
    @test isapprox(mean(epochs_no_interval.data[2].Ch1), 0.0; atol = 1e-9)
    @test isapprox(epochs_no_interval.data[1].Ch1, original_epoch1_ch1_no_int .- baseline_mean_epoch1_no_int; atol = 1e-9)
    @test isapprox(epochs_no_interval.data[2].Ch1, original_epoch2_ch1_no_int .- baseline_mean_epoch2_no_int; atol = 1e-9)

    # Test non-mutating version without interval on EpochData
    epochs = EegFun.create_test_epoch_data(n = 20, n_epochs = 2)
    original_epoch_ch1 = copy(epochs.data[1].Ch1)
    baseline_mean_epoch = mean(original_epoch_ch1)

    epochs_baselined = EegFun.baseline(epochs)

    @test isapprox(mean(epochs_baselined.data[1].Ch1), 0.0; atol = 1e-9)
    @test isapprox(epochs_baselined.data[1].Ch1, original_epoch_ch1 .- baseline_mean_epoch; atol = 1e-9)
    @test epochs_baselined !== epochs  # Should be different object

end
