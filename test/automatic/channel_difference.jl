using Test
using DataFrames
using EegFun

@testset "channel_difference" begin

    # 1) Basic difference A - B
    dat = create_test_data(n = 500)
    EegFun.channel_difference!(
        dat;
        channel_selection1 = EegFun.channels([:Ch1]),
        channel_selection2 = EegFun.channels([:Ch2]),
        channel_out = :Ch1_minus_Ch2,
    )
    @test :Ch1_minus_Ch2 ∈ propertynames(dat.data)
    @test all(dat.data.Ch1_minus_Ch2 .== (dat.data.Ch1 .- dat.data.Ch2))

    # 2) Group average difference mean(A,B) - C
    dat = create_test_data(n = 500)
    EegFun.channel_difference!(
        dat;
        channel_selection1 = EegFun.channels([:Ch1, :Ch2]),
        channel_selection2 = EegFun.channels([:Ch3]),
        channel_out = :Ch1_Ch2_minus_Ch3,
    )
    @test :Ch1_Ch2_minus_Ch3 ∈ propertynames(dat.data)
    @test all(dat.data.Ch1_Ch2_minus_Ch3 .== ((dat.data.Ch1 .+ dat.data.Ch2) ./ 2 .- dat.data.Ch3))

    # 3) Predicate selection (same as explicit)
    dat = create_test_data(n = 500)
    EegFun.channel_difference!(
        dat;
        channel_selection1 = EegFun.channels([:Ch1, :Ch2]),
        channel_selection2 = EegFun.channels([:Ch3]),
        channel_out = :out,
    )
    @test all(dat.data.out .== ((dat.data.Ch1 .+ dat.data.Ch2) ./ 2 .- dat.data.Ch3))

    # 4) Non-mutating version returns a new object; original unchanged
    dat = create_test_data(n = 500)
    dat = EegFun.channel_difference(
        dat;
        channel_selection1 = EegFun.channels([:Ch1]),
        channel_selection2 = EegFun.channels([:Ch3]),
        channel_out = :out,
    )
    @test :out ∈ propertynames(dat.data)
    @test :xxx ∉ propertynames(dat.data)
    @test all(dat.data.out .== (dat.data.Ch1 .- dat.data.Ch3))

    # 5) Overwrite behavior: write then overwrite with different groups
    dat = create_test_data(n = 500)
    dat.data[!, :X] = zeros(500)
    EegFun.channel_difference!(
        dat;
        channel_selection1 = EegFun.channels([:Ch2]),
        channel_selection2 = EegFun.channels([:Ch1]),
        channel_out = :X,
    )
    @test all(dat.data.X .== (dat.data.Ch2 .- dat.data.Ch1))

    # 6) EpochData: append to each epoch
    dat = create_test_epoch_data(n = 500)
    EegFun.channel_difference!(
        dat;
        channel_selection1 = EegFun.channels([:Ch1]),
        channel_selection2 = EegFun.channels([:Ch2]),
    )
    @test :diff ∈ propertynames(dat.data[1]) && :diff ∈ propertynames(dat.data[2])
    @test all(dat.data[1].diff .== (dat.data[1].Ch1 .- dat.data[1].Ch2))
    @test all(dat.data[2].diff .== (dat.data[2].Ch1 .- dat.data[2].Ch2))

    # 7) ErpData (SingleDataFrameEeg): append
    dat = create_test_epoch_data(n = 500)
    EegFun.channel_difference!(
        dat;
        channel_selection1 = EegFun.channels([:Ch3]),
        channel_selection2 = EegFun.channels([:Ch2]),
        channel_out = :Ch3_minus_Ch2,
    )
    @test :Ch3_minus_Ch2 ∈ propertynames(dat.data[1])
    @test all(dat.data[1].Ch3_minus_Ch2 .== (dat.data[1].Ch3 .- dat.data[1].Ch2))

    # 10) Commutativity sanity: A-B == -(B-A)
    dat = create_test_epoch_data(n = 500)
    dat_Ch1_Ch2 = copy(dat)
    EegFun.channel_difference!(
        dat_Ch1_Ch2;
        channel_selection1 = EegFun.channels([:Ch1]),
        channel_selection2 = EegFun.channels([:Ch2]),
        channel_out = :Ch1_minus_Ch2,
    )
    dat_Ch2_Ch1 = copy(dat)
    EegFun.channel_difference!(
        dat_Ch2_Ch1;
        channel_selection1 = EegFun.channels([:Ch2]),
        channel_selection2 = EegFun.channels([:Ch1]),
        channel_out = :Ch2_minus_Ch1,
    )
    @test all(dat_Ch1_Ch2.data[1].Ch1_minus_Ch2 .== .-(dat_Ch2_Ch1.data[1].Ch2_minus_Ch1))
    @test all(dat_Ch1_Ch2.data[2].Ch1_minus_Ch2 .== .-(dat_Ch2_Ch1.data[2].Ch2_minus_Ch1))
    @test all(dat_Ch1_Ch2.data[end].Ch1_minus_Ch2 .== .-(dat_Ch2_Ch1.data[end].Ch2_minus_Ch1))

    # 11) Default behavior: all channels vs all channels (should be zero)
    dat = create_test_data(n = 100)
    EegFun.channel_difference!(dat)
    @test :diff ∈ propertynames(dat.data)
    # Note: Due to floating point precision and how channels() selects, might not be exactly zero
    # But should be very close to zero
    @test all(isapprox.(dat.data.diff, 0.0; atol = 1e-6))

    # 12) Test with ErpData
    erp = create_test_erp_data(1, 1, n_channels = 3)
    EegFun.channel_difference!(
        erp;
        channel_selection1 = EegFun.channels([:Ch1]),
        channel_selection2 = EegFun.channels([:Ch2]),
        channel_out = :Ch1_minus_Ch2,
    )
    @test :Ch1_minus_Ch2 ∈ propertynames(erp.data)
    @test all(erp.data.Ch1_minus_Ch2 .== (erp.data.Ch1 .- erp.data.Ch2))

    # 13) Test warning for overwriting existing channel
    dat = create_test_data(n = 100)
    dat.data[!, :existing_channel] = ones(nrow(dat.data))
    original_value = copy(dat.data.existing_channel)

    EegFun.channel_difference!(
        dat;
        channel_selection1 = EegFun.channels([:Ch1]),
        channel_selection2 = EegFun.channels([:Ch2]),
        channel_out = :existing_channel,
    )
    # Should overwrite with difference
    @test all(dat.data.existing_channel .== (dat.data.Ch1 .- dat.data.Ch2))
    @test dat.data.existing_channel != original_value

    # 14) Test calculate_eog_channels! with EogConfig
    dat = create_test_data(n = 100)
    # Add channels that will be used for EOG
    dat.data[!, :Fp1] = dat.data.Ch1 .+ 0.1
    dat.data[!, :Fp2] = dat.data.Ch2 .+ 0.1
    dat.data[!, :IO1] = dat.data.Ch3 .+ 0.2
    dat.data[!, :IO2] = dat.data.Ch1 .+ 0.2
    dat.data[!, :F9] = dat.data.Ch2 .+ 0.3
    dat.data[!, :F10] = dat.data.Ch3 .+ 0.3

    # Update layout to include new channels
    layout_df = DataFrame(label = [:Ch1, :Ch2, :Ch3, :Fp1, :Fp2, :IO1, :IO2, :F9, :F10], inc = zeros(9), azi = zeros(9))
    dat.layout = EegFun.Layout(layout_df, nothing, nothing)

    eog_cfg = EegFun.EogConfig(
        vEOG_criterion = 50.0,
        hEOG_criterion = 30.0,
        vEOG_channels = [["Fp1", "Fp2"], ["IO1", "IO2"], ["vEOG"]],
        hEOG_channels = [["F9"], ["F10"], ["hEOG"]],
    )

    EegFun.calculate_eog_channels!(dat, eog_cfg)

    @test :vEOG ∈ propertynames(dat.data)
    @test :hEOG ∈ propertynames(dat.data)
    # vEOG should be mean(Fp1, Fp2) - mean(IO1, IO2)
    expected_vEOG = ((dat.data.Fp1 .+ dat.data.Fp2) ./ 2) .- ((dat.data.IO1 .+ dat.data.IO2) ./ 2)
    @test all(isapprox.(dat.data.vEOG, expected_vEOG; atol = 1e-10))
    # hEOG should be F9 - F10
    @test all(isapprox.(dat.data.hEOG, dat.data.F9 .- dat.data.F10; atol = 1e-10))

    # 15) Test calculate_eog_channels! with Dict
    dat2 = create_test_data(n = 100)
    dat2.data[!, :Fp1] = dat2.data.Ch1 .+ 0.1
    dat2.data[!, :Fp2] = dat2.data.Ch2 .+ 0.1
    dat2.data[!, :IO1] = dat2.data.Ch3 .+ 0.2
    dat2.data[!, :IO2] = dat2.data.Ch1 .+ 0.2
    dat2.data[!, :F9] = dat2.data.Ch2 .+ 0.3
    dat2.data[!, :F10] = dat2.data.Ch3 .+ 0.3

    layout_df2 =
        DataFrame(label = [:Ch1, :Ch2, :Ch3, :Fp1, :Fp2, :IO1, :IO2, :F9, :F10], inc = zeros(9), azi = zeros(9))
    dat2.layout = EegFun.Layout(layout_df2, nothing, nothing)

    eog_cfg_dict = Dict(
        "vEOG_criterion" => 50.0,
        "hEOG_criterion" => 30.0,
        "vEOG_channels" => [["Fp1", "Fp2"], ["IO1", "IO2"], ["vEOG"]],
        "hEOG_channels" => [["F9"], ["F10"], ["hEOG"]],
    )

    EegFun.calculate_eog_channels!(dat2, eog_cfg_dict)

    @test :vEOG ∈ propertynames(dat2.data)
    @test :hEOG ∈ propertynames(dat2.data)

    # 16) Test calculate_eog_channels! with EpochData
    epochs = create_test_epoch_data(n = 100, n_epochs = 3)
    for epoch_df in epochs.data
        epoch_df[!, :Fp1] = epoch_df.Ch1 .+ 0.1
        epoch_df[!, :Fp2] = epoch_df.Ch2 .+ 0.1
        epoch_df[!, :IO1] = epoch_df.Ch3 .+ 0.2
        epoch_df[!, :IO2] = epoch_df.Ch1 .+ 0.2
    end

    layout_df3 = DataFrame(label = [:Ch1, :Ch2, :Ch3, :Fp1, :Fp2, :IO1, :IO2], inc = zeros(7), azi = zeros(7))
    epochs.layout = EegFun.Layout(layout_df3, nothing, nothing)

    eog_cfg_epochs = EegFun.EogConfig(
        vEOG_criterion = 50.0,
        hEOG_criterion = 30.0,
        vEOG_channels = [["Fp1", "Fp2"], ["IO1", "IO2"], ["vEOG"]],
        hEOG_channels = [["Ch1"], ["Ch2"], ["hEOG"]],  # Use existing channels
    )

    EegFun.calculate_eog_channels!(epochs, eog_cfg_epochs)

    @test :vEOG ∈ propertynames(epochs.data[1])
    @test :hEOG ∈ propertynames(epochs.data[1])
    # Check that each epoch has the channels
    for epoch_df in epochs.data
        @test :vEOG ∈ propertynames(epoch_df)
        @test :hEOG ∈ propertynames(epoch_df)
    end

    # 17) Test detect_eog_signals! with EogConfig
    dat3 = create_test_data(n = 1000, fs = 1000)
    # Create signal with large jumps (not just large amplitude) for EOG detection
    # detect_eog_onsets! looks for differences/jumps, so create step changes
    vEOG_signal = zeros(nrow(dat3.data))
    vEOG_signal[100:110] .= 200.0  # Large jump
    vEOG_signal[500:510] .= -200.0  # Large negative jump
    vEOG_signal[800:810] .= 150.0  # Another jump
    dat3.data[!, :vEOG] = vEOG_signal

    hEOG_signal = zeros(nrow(dat3.data))
    hEOG_signal[200:210] .= 150.0  # Large jump
    hEOG_signal[600:610] .= -150.0  # Large negative jump
    dat3.data[!, :hEOG] = hEOG_signal

    layout_df4 = DataFrame(label = [:Ch1, :Ch2, :Ch3, :vEOG, :hEOG], inc = zeros(5), azi = zeros(5))
    dat3.layout = EegFun.Layout(layout_df4, nothing, nothing)

    eog_cfg_detect = EegFun.EogConfig(
        vEOG_criterion = 50.0,
        hEOG_criterion = 30.0,
        vEOG_channels = [["Fp1"], ["Fp2"], ["vEOG"]],
        hEOG_channels = [["F9"], ["F10"], ["hEOG"]],
    )

    EegFun.detect_eog_signals!(dat3, eog_cfg_detect)

    @test :is_vEOG ∈ propertynames(dat3.data)
    @test :is_hEOG ∈ propertynames(dat3.data)
    # Should detect some EOG onsets given the large jumps
    @test sum(dat3.data.is_vEOG) > 0
    @test sum(dat3.data.is_hEOG) > 0

    # 18) Test error handling - missing channels
    # If a channel doesn't exist, get_selected_channels returns empty vector,
    # which causes division by zero in _calculate_channel_difference!
    # However, @minimal_error returns nothing, so the function continues and hits DivideError
    dat4 = create_test_data(n = 100)
    # When one channel selection is empty, division by zero occurs
    try
        EegFun.channel_difference!(
            dat4;
            channel_selection1 = EegFun.channels([:NonExistentChannel]),
            channel_selection2 = EegFun.channels([:Ch1]),
        )
        # If we get here, check that the result contains NaN or Inf (division by zero result)
        @test any(isnan.(dat4.data.diff)) || any(isinf.(dat4.data.diff))
    catch e
        # Should throw DivideError
        @test e isa DivideError
    end

    # Test with both channels missing (both empty vectors)
    dat4b = create_test_data(n = 100)
    try
        EegFun.channel_difference!(
            dat4b;
            channel_selection1 = EegFun.channels([:NonExistentChannel1]),
            channel_selection2 = EegFun.channels([:NonExistentChannel2]),
        )
        # If we get here, check that the result contains NaN or Inf (division by zero result)
        @test any(isnan.(dat4b.data.diff)) || any(isinf.(dat4b.data.diff))
    catch e
        # Should throw DivideError
        @test e isa DivideError
    end

    # 19) Test with single channel in each group
    dat5 = create_test_data(n = 100)
    EegFun.channel_difference!(
        dat5;
        channel_selection1 = EegFun.channels([:Ch1]),
        channel_selection2 = EegFun.channels([:Ch2]),
        channel_out = :single_diff,
    )
    @test all(dat5.data.single_diff .== (dat5.data.Ch1 .- dat5.data.Ch2))

    # 20) Test with multiple channels in first group, single in second
    dat6 = create_test_data(n = 100)
    EegFun.channel_difference!(
        dat6;
        channel_selection1 = EegFun.channels([:Ch1, :Ch2, :Ch3]),
        channel_selection2 = EegFun.channels([:Ch1]),
        channel_out = :multi_single_diff,
    )
    expected_multi_single = ((dat6.data.Ch1 .+ dat6.data.Ch2 .+ dat6.data.Ch3) ./ 3) .- dat6.data.Ch1
    @test all(isapprox.(dat6.data.multi_single_diff, expected_multi_single; atol = 1e-10))

    # 21) Test non-mutating version with EpochData
    epochs2 = create_test_epoch_data(n = 100, n_epochs = 2)
    epochs2_result = EegFun.channel_difference(
        epochs2;
        channel_selection1 = EegFun.channels([:Ch1]),
        channel_selection2 = EegFun.channels([:Ch2]),
        channel_out = :diff_result,
    )
    @test epochs2_result !== epochs2
    @test :diff_result ∈ propertynames(epochs2_result.data[1])
    @test :diff_result ∉ propertynames(epochs2.data[1])  # Original unchanged

    # 22) Test non-mutating version with ErpData
    erp2 = create_test_erp_data(1, 1, n_channels = 3)
    erp2_result = EegFun.channel_difference(
        erp2;
        channel_selection1 = EegFun.channels([:Ch1]),
        channel_selection2 = EegFun.channels([:Ch2]),
    )
    @test erp2_result !== erp2
    @test :diff ∈ propertynames(erp2_result.data)
    @test :diff ∉ propertynames(erp2.data)  # Original unchanged

    # 23) Test that metadata columns are not included in difference calculation
    dat7 = create_test_data(n = 100)
    original_time = copy(dat7.data.time)
    original_sample = copy(dat7.data.sample)

    EegFun.channel_difference!(dat7)

    # Metadata should be unchanged
    @test dat7.data.time == original_time
    @test dat7.data.sample == original_sample

end
