using Test
using DataFrames
using EegFun

@testset "channel_difference" begin

    # Basic difference A - B
    dat = EegFun.create_test_continuous_data(n = 500)
    EegFun.channel_difference!(
        dat;
        channel_selection1 = EegFun.channels([:Ch1]),
        channel_selection2 = EegFun.channels([:Ch2]),
        channel_out = :Ch1_minus_Ch2,
    )
    @test :Ch1_minus_Ch2 ∈ propertynames(dat.data)
    @test all(dat.data.Ch1_minus_Ch2 .== (dat.data.Ch1 .- dat.data.Ch2))

    # Group average difference mean(A,B) - C
    dat = EegFun.create_test_continuous_data(n = 500)
    EegFun.channel_difference!(
        dat;
        channel_selection1 = EegFun.channels([:Ch1, :Ch2]),
        channel_selection2 = EegFun.channels([:Ch3]),
        channel_out = :Ch1_Ch2_minus_Ch3,
    )
    @test :Ch1_Ch2_minus_Ch3 ∈ propertynames(dat.data)
    @test all(dat.data.Ch1_Ch2_minus_Ch3 .== ((dat.data.Ch1 .+ dat.data.Ch2) ./ 2 .- dat.data.Ch3))

    # Non-mutating version returns a new object; original unchanged
    dat = EegFun.create_test_continuous_data(n = 500)
    dat = EegFun.channel_difference(
        dat;
        channel_selection1 = EegFun.channels([:Ch1]),
        channel_selection2 = EegFun.channels([:Ch3]),
        channel_out = :out,
    )
    @test :out ∈ propertynames(dat.data)
    @test :xxx ∉ propertynames(dat.data)
    @test all(dat.data.out .== (dat.data.Ch1 .- dat.data.Ch3))

    # Overwrite behavior: write then overwrite with different groups
    dat = EegFun.create_test_continuous_data(n = 500)
    dat.data[!, :X] = zeros(500)
    EegFun.channel_difference!(
        dat;
        channel_selection1 = EegFun.channels([:Ch2]),
        channel_selection2 = EegFun.channels([:Ch1]),
        channel_out = :X,
    )
    @test all(dat.data.X .== (dat.data.Ch2 .- dat.data.Ch1))

    # EpochData: append to each epoch
    dat = EegFun.create_test_epoch_data(n = 500)
    EegFun.channel_difference!(dat; channel_selection1 = EegFun.channels([:Ch1]), channel_selection2 = EegFun.channels([:Ch2]))
    @test :diff ∈ propertynames(dat.data[1]) && :diff ∈ propertynames(dat.data[2])
    @test all(dat.data[1].diff .== (dat.data[1].Ch1 .- dat.data[1].Ch2))
    @test all(dat.data[2].diff .== (dat.data[2].Ch1 .- dat.data[2].Ch2))

    # ErpData (SingleDataFrameEeg): append
    dat = EegFun.create_test_epoch_data(n = 500)
    EegFun.channel_difference!(
        dat;
        channel_selection1 = EegFun.channels([:Ch3]),
        channel_selection2 = EegFun.channels([:Ch2]),
        channel_out = :Ch3_minus_Ch2,
    )
    @test :Ch3_minus_Ch2 ∈ propertynames(dat.data[1])
    @test all(dat.data[1].Ch3_minus_Ch2 .== (dat.data[1].Ch3 .- dat.data[1].Ch2))

    # Commutativity sanity: A-B == -(B-A)
    dat = EegFun.create_test_epoch_data(n = 500)
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

    # Default behavior: all channels vs all channels (should be zero)
    dat = EegFun.create_test_continuous_data(n = 100)
    EegFun.channel_difference!(dat)
    @test :diff ∈ propertynames(dat.data)
    # should be very close to zero
    @test all(isapprox.(dat.data.diff, 0.0; atol = 1e-6))

    # Test with ErpData
    erp = EegFun.create_test_erp_data(participant = 1, condition = 1, n_channels = 3)
    EegFun.channel_difference!(
        erp;
        channel_selection1 = EegFun.channels([:Ch1]),
        channel_selection2 = EegFun.channels([:Ch2]),
        channel_out = :Ch1_minus_Ch2,
    )
    @test :Ch1_minus_Ch2 ∈ propertynames(erp.data)
    @test all(erp.data.Ch1_minus_Ch2 .== (erp.data.Ch1 .- erp.data.Ch2))

    # Test warning for overwriting existing channel
    dat = EegFun.create_test_continuous_data(n = 100)
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

    # Test calculate_eog_channels! with EogConfig
    dat = EegFun.create_test_continuous_data(n = 100)
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

    # Test calculate_eog_channels! with Dict
    dat = EegFun.create_test_continuous_data(n = 100)
    dat.data[!, :Fp1] = dat.data.Ch1 .+ 0.1
    dat.data[!, :Fp2] = dat.data.Ch2 .+ 0.1
    dat.data[!, :IO1] = dat.data.Ch3 .+ 0.2
    dat.data[!, :IO2] = dat.data.Ch1 .+ 0.2
    dat.data[!, :F9] = dat.data.Ch2 .+ 0.3
    dat.data[!, :F10] = dat.data.Ch3 .+ 0.3

    layout_df = DataFrame(label = [:Ch1, :Ch2, :Ch3, :Fp1, :Fp2, :IO1, :IO2, :F9, :F10], inc = zeros(9), azi = zeros(9))
    dat.layout = EegFun.Layout(layout_df, nothing, nothing)

    eog_cfg_dict = Dict(
        "vEOG_criterion" => 50.0,
        "hEOG_criterion" => 30.0,
        "vEOG_channels" => [["Fp1", "Fp2"], ["IO1", "IO2"], ["vEOG"]],
        "hEOG_channels" => [["F9"], ["F10"], ["hEOG"]],
    )

    EegFun.calculate_eog_channels!(dat, eog_cfg_dict)

    @test :vEOG ∈ propertynames(dat.data)
    @test :hEOG ∈ propertynames(dat.data)

    # Test calculate_eog_channels! with EpochData
    epochs = EegFun.create_test_epoch_data(n = 100, n_epochs = 3)
    for epoch_df in epochs.data
        epoch_df[!, :Fp1] = epoch_df.Ch1 .+ 0.1
        epoch_df[!, :Fp2] = epoch_df.Ch2 .+ 0.1
        epoch_df[!, :IO1] = epoch_df.Ch3 .+ 0.2
        epoch_df[!, :IO2] = epoch_df.Ch1 .+ 0.2
    end

    layout_df = DataFrame(label = [:Ch1, :Ch2, :Ch3, :Fp1, :Fp2, :IO1, :IO2], inc = zeros(7), azi = zeros(7))
    epochs.layout = EegFun.Layout(layout_df, nothing, nothing)

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

    # Test detect_eog_signals! with EogConfig
    dat = EegFun.create_test_continuous_data(n = 1000, fs = 1000)
    # Create signal with large jumps (not just large amplitude) for EOG detection
    # detect_eog_onsets! looks for differences/jumps, so create step changes
    vEOG_signal = zeros(nrow(dat.data))
    vEOG_signal[100:110] .= 200.0  # Large jump
    vEOG_signal[500:510] .= -200.0  # Large negative jump
    vEOG_signal[800:810] .= 150.0  # Another jump
    dat.data[!, :vEOG] = vEOG_signal

    hEOG_signal = zeros(nrow(dat.data))
    hEOG_signal[200:210] .= 150.0  # Large jump
    hEOG_signal[600:610] .= -150.0  # Large negative jump
    dat.data[!, :hEOG] = hEOG_signal

    layout_df = DataFrame(label = [:Ch1, :Ch2, :Ch3, :vEOG, :hEOG], inc = zeros(5), azi = zeros(5))
    dat.layout = EegFun.Layout(layout_df, nothing, nothing)

    eog_cfg_detect = EegFun.EogConfig(
        vEOG_criterion = 50.0,
        hEOG_criterion = 30.0,
        vEOG_channels = [["Fp1"], ["Fp2"], ["vEOG"]],
        hEOG_channels = [["F9"], ["F10"], ["hEOG"]],
    )

    EegFun.detect_eog_signals!(dat, eog_cfg_detect)

    @test :is_vEOG ∈ propertynames(dat.data)
    @test :is_hEOG ∈ propertynames(dat.data)
    # Should detect some EOG onsets given the large jumps
    @test sum(dat.data.is_vEOG) > 0
    @test sum(dat.data.is_hEOG) > 0

    # Test with multiple channels in first group, single in second
    dat = EegFun.create_test_continuous_data(n = 100)
    EegFun.channel_difference!(
        dat;
        channel_selection1 = EegFun.channels([:Ch1, :Ch2, :Ch3]),
        channel_selection2 = EegFun.channels([:Ch1]),
        channel_out = :multi_single_diff,
    )
    expected_multi_single = ((dat.data.Ch1 .+ dat.data.Ch2 .+ dat.data.Ch3) ./ 3) .- dat.data.Ch1
    @test all(isapprox.(dat.data.multi_single_diff, expected_multi_single; atol = 1e-10))

    # Non-mutating version with EpochData
    epochs = EegFun.create_test_epoch_data(n = 100, n_epochs = 2)
    epochs_result = EegFun.channel_difference(
        epochs;
        channel_selection1 = EegFun.channels([:Ch1]),
        channel_selection2 = EegFun.channels([:Ch2]),
        channel_out = :diff_result,
    )
    @test epochs_result !== epochs
    @test :diff_result ∈ propertynames(epochs_result.data[1])
    @test :diff_result ∉ propertynames(epochs.data[1])  # Original unchanged

    # Non-mutating version with ErpData
    erp = EegFun.create_test_erp_data(participant = 1, condition = 1, n_channels = 3)
    erp_result = EegFun.channel_difference(erp; channel_selection1 = EegFun.channels([:Ch1]), channel_selection2 = EegFun.channels([:Ch2]))
    @test erp_result !== erp
    @test :diff ∈ propertynames(erp_result.data)
    @test :diff ∉ propertynames(erp.data)  # Original unchanged

    # Test that metadata columns are not included in difference calculation
    dat = EegFun.create_test_continuous_data(n = 100)
    original_time = copy(dat.data.time)
    original_sample = copy(dat.data.sample)

    EegFun.channel_difference!(dat)

    # Metadata should be unchanged
    @test dat.data.time == original_time
    @test dat.data.sample == original_sample

end

@testset "TF Channel Difference" begin

    @testset "Mutating: basic difference" begin
        tf = EegFun.create_test_tf_data(n_channels = 3, power_offset = 5.0, noise = 0.0)

        EegFun.channel_difference!(
            tf,
            channel_selection1 = EegFun.channels([:Ch1]),
            channel_selection2 = EegFun.channels([:Ch2]),
            channel_out = :diff_12,
        )

        # Should have :diff_12 in both power and phase
        @test hasproperty(tf.data_power, :diff_12)
        @test hasproperty(tf.data_phase, :diff_12)

        # With same offset and no noise, difference should be ~0
        @test all(abs.(tf.data_power[!, :diff_12]) .< 1e-10)
    end

    @testset "Mutating: group difference" begin
        tf = EegFun.create_test_tf_data(n_channels = 4, power_offset = 3.0, noise = 0.0)

        EegFun.channel_difference!(
            tf,
            channel_selection1 = EegFun.channels([:Ch1, :Ch2]),
            channel_selection2 = EegFun.channels([:Ch3, :Ch4]),
            channel_out = :laterality,
        )

        @test hasproperty(tf.data_power, :laterality)
        @test hasproperty(tf.data_phase, :laterality)
    end

    @testset "Non-mutating: returns copy" begin
        tf = EegFun.create_test_tf_data(n_channels = 3)

        result = EegFun.channel_difference(
            tf,
            channel_selection1 = EegFun.channels([:Ch1]),
            channel_selection2 = EegFun.channels([:Ch2]),
            channel_out = :test_diff,
        )

        # Original should not be modified
        @test !hasproperty(tf.data_power, :test_diff)

        # Result should have the diff column
        @test hasproperty(result.data_power, :test_diff)
        @test hasproperty(result.data_phase, :test_diff)
    end

    @testset "Non-mutating: Vector{TimeFreqData}" begin
        tfs = [EegFun.create_test_tf_data(condition = 1), EegFun.create_test_tf_data(condition = 2)]

        results = EegFun.channel_difference(
            tfs,
            channel_selection1 = EegFun.channels([:Ch1]),
            channel_selection2 = EegFun.channels([:Ch2]),
            channel_out = :diff,
        )

        @test length(results) == 2
        @test all(hasproperty(r.data_power, :diff) for r in results)
        @test all(hasproperty(r.data_phase, :diff) for r in results)

        # Originals untouched
        @test !hasproperty(tfs[1].data_power, :diff)
    end

    @testset "Both DataFrames stay in sync" begin
        tf = EegFun.create_test_tf_data(n_channels = 3)

        EegFun.channel_difference!(
            tf,
            channel_selection1 = EegFun.channels([:Ch1]),
            channel_selection2 = EegFun.channels([:Ch2]),
            channel_out = :sync_test,
        )

        # Same number of rows
        @test nrow(tf.data_power) == nrow(tf.data_phase)

        # Same columns
        @test propertynames(tf.data_power) == propertynames(tf.data_phase)
    end
end
