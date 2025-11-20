using Test
using DataFrames
using Statistics
using eegfun

@testset "baseline" begin

    # 1) Baseline over first sample for Ch1, Ch2, and Ch3
    dat = create_test_data(n = 6)
    eegfun.baseline!(dat, eegfun.IntervalIndex(start = 1, stop = 1))

    @test isapprox(mean(dat.data.Ch1[1]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch2[1]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch3[1]), 0.0; atol = 1e-9)

    dat = create_test_data(n = 6)
    dat = eegfun.baseline(dat, eegfun.IntervalIndex(start = 1, stop = 1))

    @test isapprox(mean(dat.data.Ch1[1]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch2[1]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch3[1]), 0.0; atol = 1e-9)

    # 2) Baseline over entire range for all channels
    dat = create_test_data(n = 6)
    eegfun.baseline!(dat)

    @test isapprox(mean(dat.data.Ch1), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch2), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch3), 0.0; atol = 1e-9)

    dat = create_test_data(n = 6)
    dat = eegfun.baseline(dat)

    @test isapprox(mean(dat.data.Ch1), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch2), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch3), 0.0; atol = 1e-9)


    # 3) EpochData: each epoch baselined independently
    dat = create_test_epoch_data(n = 3)
    eegfun.baseline!(dat, eegfun.IntervalIndex(start = 1, stop = 1))

    @test isapprox(mean(dat.data[1].Ch1[1]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data[2].Ch2[1]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data[3].Ch3[1]), 0.0; atol = 1e-9)

    dat = create_test_epoch_data(n = 3)
    eegfun.baseline!(dat, eegfun.IntervalIndex(start = 1, stop = 3))

    @test isapprox(mean(dat.data[1].Ch1[1:3]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data[2].Ch2[1:3]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data[3].Ch3[1:3]), 0.0; atol = 1e-9)

    # 4) IntervalTime converted to indices correctly
    dat = create_test_data(n = 6)
    eegfun.baseline!(dat, eegfun.IntervalTime(start = 0.0, stop = 0.0))

    @test isapprox(mean(dat.data.Ch1[1]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch2[1]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch3[1]), 0.0; atol = 1e-9)

    # 4) IntervalTime converted to indices correctly
    dat = create_test_data(n = 6)
    dat = eegfun.baseline(dat, eegfun.IntervalTime(start = 0.0, stop = 0.0))

    @test isapprox(mean(dat.data.Ch1[1]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch2[1]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch3[1]), 0.0; atol = 1e-9)


    # 4) IntervalTime converted to indices correctly
    dat = create_test_data(n = 6)
    eegfun.baseline!(dat, eegfun.IntervalTime(start = 0.003, stop = 0.003))

    @test isapprox(mean(dat.data.Ch1[4]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch2[4]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch3[4]), 0.0; atol = 1e-9)

    dat = create_test_data(n = 6)
    dat = eegfun.baseline(dat, eegfun.IntervalTime(start = 0.003, stop = 0.003))

    @test isapprox(mean(dat.data.Ch1[4]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch2[4]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch3[4]), 0.0; atol = 1e-9)


    dat = create_test_epoch_data(n = 6)
    eegfun.baseline!(dat, eegfun.IntervalTime(start = 0.003, stop = 0.003))

    @test isapprox(mean(dat.data[1].Ch1[4]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data[2].Ch2[4]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data[3].Ch3[4]), 0.0; atol = 1e-9)

    dat = create_test_epoch_data(n = 6)
    dat = eegfun.baseline(dat, eegfun.IntervalTime(start = 0.003, stop = 0.003))

    @test isapprox(mean(dat.data[1].Ch1[4]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data[2].Ch2[4]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data[3].Ch3[4]), 0.0; atol = 1e-9)

    # Test channel selection
    dat = create_test_data(n = 10)
    original_ch1 = copy(dat.data.Ch1)
    original_ch2 = copy(dat.data.Ch2)
    original_ch3 = copy(dat.data.Ch3)
    
    # Baseline only Ch1 and Ch2
    eegfun.baseline!(dat, eegfun.IntervalIndex(start = 1, stop = 5), channel_selection = eegfun.channels([:Ch1, :Ch2]))
    
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
    dat2 = create_test_data(n = 10)
    original_ch1_2 = copy(dat2.data.Ch1)
    baseline_mean_ch1_2 = mean(original_ch1_2[1:5])
    
    dat2_baselined = eegfun.baseline(dat2, eegfun.IntervalIndex(start = 1, stop = 5), channel_selection = eegfun.channels([:Ch1]))
    
    @test isapprox(mean(dat2_baselined.data.Ch1[1:5]), 0.0; atol = 1e-9)
    @test isapprox(dat2_baselined.data.Ch1, original_ch1_2 .- baseline_mean_ch1_2; atol = 1e-9)
    @test dat2_baselined !== dat2  # Should be different object
    
    # Test empty channel selection (should warn and return early)
    dat3 = create_test_data(n = 10)
    original_ch1_3 = copy(dat3.data.Ch1)
    
    eegfun.baseline!(dat3, eegfun.IntervalIndex(start = 1, stop = 5), channel_selection = eegfun.channels(Symbol[]))
    
    # Data should be unchanged
    @test dat3.data.Ch1 == original_ch1_3
    
    # Test with EpochData and channel selection
    epochs = create_test_epoch_data(n = 10, n_epochs = 3)
    original_epoch1_ch1 = copy(epochs.data[1].Ch1)
    baseline_mean_epoch1 = mean(original_epoch1_ch1[1:5])
    
    eegfun.baseline!(epochs, eegfun.IntervalIndex(start = 1, stop = 5), channel_selection = eegfun.channels([:Ch1]))
    
    @test isapprox(mean(epochs.data[1].Ch1[1:5]), 0.0; atol = 1e-9)
    @test isapprox(epochs.data[1].Ch1, original_epoch1_ch1 .- baseline_mean_epoch1; atol = 1e-9)
    
    # Test with ErpData
    erp = create_test_erp_data(1, 1, n_channels = 3)
    original_erp_ch1 = copy(erp.data.Ch1)
    baseline_mean_erp = mean(original_erp_ch1[1:10])
    
    eegfun.baseline!(erp, eegfun.IntervalIndex(start = 1, stop = 10))
    
    @test isapprox(mean(erp.data.Ch1[1:10]), 0.0; atol = 1e-9)
    @test isapprox(erp.data.Ch1, original_erp_ch1 .- baseline_mean_erp; atol = 1e-9)
    
    # Test with ErpData and channel selection
    erp2 = create_test_erp_data(1, 1, n_channels = 3)
    original_erp2_ch1 = copy(erp2.data.Ch1)
    original_erp2_ch2 = copy(erp2.data.Ch2)
    baseline_mean_erp2 = mean(original_erp2_ch1[1:10])
    
    eegfun.baseline!(erp2, eegfun.IntervalIndex(start = 1, stop = 10), channel_selection = eegfun.channels([:Ch1]))
    
    @test isapprox(mean(erp2.data.Ch1[1:10]), 0.0; atol = 1e-9)
    @test isapprox(erp2.data.Ch1, original_erp2_ch1 .- baseline_mean_erp2; atol = 1e-9)
    @test erp2.data.Ch2 == original_erp2_ch2  # Ch2 unchanged
    
    # Test baseline! without interval (uses entire range)
    dat4 = create_test_data(n = 20)
    original_ch1_4 = copy(dat4.data.Ch1)
    baseline_mean_ch1_4 = mean(original_ch1_4)
    
    eegfun.baseline!(dat4)
    
    @test isapprox(mean(dat4.data.Ch1), 0.0; atol = 1e-9)
    @test isapprox(dat4.data.Ch1, original_ch1_4 .- baseline_mean_ch1_4; atol = 1e-9)
    
    # Test baseline (non-mutating) without interval
    dat5 = create_test_data(n = 20)
    original_ch1_5 = copy(dat5.data.Ch1)
    baseline_mean_ch1_5 = mean(original_ch1_5)
    
    dat5_baselined = eegfun.baseline(dat5)
    
    @test isapprox(mean(dat5_baselined.data.Ch1), 0.0; atol = 1e-9)
    @test isapprox(dat5_baselined.data.Ch1, original_ch1_5 .- baseline_mean_ch1_5; atol = 1e-9)
    @test dat5_baselined !== dat5
    
    # Test with larger baseline interval
    dat6 = create_test_data(n = 100)
    original_ch1_6 = copy(dat6.data.Ch1)
    baseline_mean_ch1_6 = mean(original_ch1_6[10:50])
    
    eegfun.baseline!(dat6, eegfun.IntervalIndex(start = 10, stop = 50))
    
    @test isapprox(mean(dat6.data.Ch1[10:50]), 0.0; atol = 1e-9)
    @test isapprox(dat6.data.Ch1, original_ch1_6 .- baseline_mean_ch1_6; atol = 1e-9)
    
    # Test IntervalTime with larger range
    dat7 = create_test_data(n = 100, fs = 1000)
    original_ch1_7 = copy(dat7.data.Ch1)
    # Find indices for time range 0.01 to 0.05 seconds
    time_idx_start = findfirst(x -> x >= 0.01, dat7.data.time)
    time_idx_stop = findlast(x -> x <= 0.05, dat7.data.time)
    baseline_mean_ch1_7 = mean(original_ch1_7[time_idx_start:time_idx_stop])
    
    eegfun.baseline!(dat7, eegfun.IntervalTime(start = 0.01, stop = 0.05))
    
    @test isapprox(mean(dat7.data.Ch1[time_idx_start:time_idx_stop]), 0.0; atol = 1e-9)
    @test isapprox(dat7.data.Ch1, original_ch1_7 .- baseline_mean_ch1_7; atol = 1e-9)
    
    # Test that baseline correction is applied to entire signal, not just baseline interval
    dat8 = create_test_data(n = 100)
    original_ch1_8 = copy(dat8.data.Ch1)
    baseline_mean_ch1_8 = mean(original_ch1_8[1:10])
    original_mean_ch1_8 = mean(original_ch1_8)
    
    eegfun.baseline!(dat8, eegfun.IntervalIndex(start = 1, stop = 10))
    
    # Mean of baseline interval should be 0
    @test isapprox(mean(dat8.data.Ch1[1:10]), 0.0; atol = 1e-9)
    # Entire signal should be shifted by baseline mean
    @test isapprox(dat8.data.Ch1, original_ch1_8 .- baseline_mean_ch1_8; atol = 1e-9)
    # Mean of entire signal should be original_mean - baseline_mean
    @test isapprox(mean(dat8.data.Ch1), original_mean_ch1_8 - baseline_mean_ch1_8; atol = 1e-9)
    
    # Test with EpochData - each epoch baselined independently
    epochs2 = create_test_epoch_data(n = 10, n_epochs = 3)
    original_epoch1_ch1_2 = copy(epochs2.data[1].Ch1)
    original_epoch2_ch1_2 = copy(epochs2.data[2].Ch1)
    original_epoch3_ch1_2 = copy(epochs2.data[3].Ch1)
    
    baseline_mean_epoch1_2 = mean(original_epoch1_ch1_2[1:5])
    baseline_mean_epoch2_2 = mean(original_epoch2_ch1_2[1:5])
    baseline_mean_epoch3_2 = mean(original_epoch3_ch1_2[1:5])
    
    eegfun.baseline!(epochs2, eegfun.IntervalIndex(start = 1, stop = 5))
    
    # Each epoch should be baselined independently
    @test isapprox(mean(epochs2.data[1].Ch1[1:5]), 0.0; atol = 1e-9)
    @test isapprox(mean(epochs2.data[2].Ch1[1:5]), 0.0; atol = 1e-9)
    @test isapprox(mean(epochs2.data[3].Ch1[1:5]), 0.0; atol = 1e-9)
    
    @test isapprox(epochs2.data[1].Ch1, original_epoch1_ch1_2 .- baseline_mean_epoch1_2; atol = 1e-9)
    @test isapprox(epochs2.data[2].Ch1, original_epoch2_ch1_2 .- baseline_mean_epoch2_2; atol = 1e-9)
    @test isapprox(epochs2.data[3].Ch1, original_epoch3_ch1_2 .- baseline_mean_epoch3_2; atol = 1e-9)
    
    # Test error handling - invalid interval (start > stop)
    # validate_baseline_interval throws ErrorException for invalid intervals
    dat9 = create_test_data(n = 10)
    @test_throws ErrorException eegfun.baseline!(dat9, eegfun.IntervalIndex(start = 5, stop = 1))
    
    # Test error handling - invalid interval (start out of range)
    dat10 = create_test_data(n = 10)
    @test_throws MethodError eegfun.baseline!(dat10, eegfun.IntervalIndex(start = 0, stop = 5))
    
    # Test error handling - invalid interval (stop out of range)
    dat11 = create_test_data(n = 10)
    @test_throws MethodError eegfun.baseline!(dat11, eegfun.IntervalIndex(start = 1, stop = 100))
    
    # Test error handling - invalid IntervalTime (outside time range)
    # When IntervalTime is outside range, find_idx_start_end may return nothing, causing issues
    dat12 = create_test_data(n = 10, fs = 1000)
    @test_throws MethodError eegfun.baseline!(dat12, eegfun.IntervalTime(start = 100.0, stop = 200.0))
    
    # Test with single sample baseline interval
    dat13 = create_test_data(n = 20)
    original_ch1_13 = copy(dat13.data.Ch1)
    baseline_value_ch1_13 = original_ch1_13[5]  # Single sample value
    
    eegfun.baseline!(dat13, eegfun.IntervalIndex(start = 5, stop = 5))
    
    @test isapprox(dat13.data.Ch1[5], 0.0; atol = 1e-9)
    @test isapprox(dat13.data.Ch1, original_ch1_13 .- baseline_value_ch1_13; atol = 1e-9)
    
    # Test that metadata columns are not affected
    dat14 = create_test_data(n = 20)
    original_time = copy(dat14.data.time)
    original_sample = copy(dat14.data.sample)
    original_triggers = copy(dat14.data.triggers)
    
    eegfun.baseline!(dat14, eegfun.IntervalIndex(start = 1, stop = 10))
    
    @test dat14.data.time == original_time
    @test dat14.data.sample == original_sample
    @test dat14.data.triggers == original_triggers

end
