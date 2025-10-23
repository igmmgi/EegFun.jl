using Test
using DataFrames
using eegfun

@testset "channel_average" begin

    dat = create_test_data(n = 500)

    # 1) Append averaged columns only (Symbols input)
    eegfun.channel_average!(dat, channel_selections = [eegfun.channels([:Ch1, :Ch2])])

    @test :Ch1_Ch2 ∈ propertynames(dat.data)
    @test all(dat.data.Ch1_Ch2 .== (dat.data.Ch1 .+ dat.data.Ch2) ./ 2)

    # 2) Reduce to only averages and create averaged layout
    dat = create_test_data(n = 500)
    dat = eegfun.channel_average(dat, channel_selections = [eegfun.channels([:Ch1, :Ch2])]; reduce = true)
    @test all(propertynames(dat.data) .== [:time, :sample, :triggers, :Ch1_Ch2])
    @test size(dat.layout.data, 1) == 1
    @test :inc ∈ propertynames(dat.layout.data)
    @test :azi ∈ propertynames(dat.layout.data)

    # 3) Append averaged channels to layout
    dat = create_test_data(n = 500)
    dat = eegfun.channel_average(dat, channel_selections = [eegfun.channels([:Ch2, :Ch3])])
    @test :Ch2_Ch3 ∈ propertynames(dat.data)
    @test any(dat.layout.data.label .== :Ch2_Ch3)

    # 4) Mixed Symbol input
    dat = create_test_data(n = 500)
    eegfun.channel_average!(dat, channel_selections = [eegfun.channels([:Ch1, :Ch3])])
    @test :Ch1_Ch3 ∈ propertynames(dat.data)
    @test :Ch1_Ch2 ∉ propertynames(dat.data)

    # 5) Auto-label :avg for all channels
    dat = create_test_data(n = 500)
    dat = eegfun.channel_average(dat, channel_selections = [eegfun.channels([:Ch1, :Ch2, :Ch3])]; reduce = true)
    @test :avg ∈ propertynames(dat.data)

    # 6) Custom output_labels applied and length mismatch errors
    dat = create_test_data(n = 500)
    eegfun.channel_average!(dat, channel_selections = [eegfun.channels([:Ch1, :Ch2])]; output_labels = [:output_label])
    @test :output_label ∈ propertynames(dat.data)
    @test_throws Any eegfun.channel_average!(
        dat,
        channel_selections = [eegfun.channels([:A, :B]), eegfun.channels([:A, :C])];
        output_labels = [:x],
    )

    # 7) Duplicate labels in layout (should accumulate)
    # First add B_C, then add again to verify rows accumulate
    dat = create_test_data(n = 500)
    eegfun.channel_average!(dat, channel_selections = [eegfun.channels([:Ch2, :Ch3])])
    n1 = sum(dat.layout.data.label .== :Ch2_Ch3)
    eegfun.channel_average!(dat, channel_selections = [eegfun.channels([:Ch2, :Ch3])])
    n2 = sum(dat.layout.data.label .== :Ch2_Ch3)
    @test n1 == 1 && n2 == 2

    # 9) EpochData append and reduce
    # Create simple EpochData (2 epochs) from dat
    dat = create_test_epoch_data(n = 500)
    dat = eegfun.channel_average(dat, channel_selections = [eegfun.channels([:Ch1, :Ch2])])

    @test :Ch1_Ch2 ∈ propertynames(dat.data[1]) && :Ch1_Ch2 ∈ propertynames(dat.data[2])
    @test :Ch1_Ch2 ∈ propertynames(dat.data[3]) && :Ch1_Ch2 ∈ propertynames(dat.data[4])

    dat = create_test_epoch_data(n = 500)
    dat = eegfun.channel_average(dat, channel_selections = [eegfun.channels([:Ch1, :Ch2])]; reduce = true)

    # Expect meta columns (leading) + A_B only
    cols = propertynames(dat.data[1])
    @test cols[end] == :Ch1_Ch2
    @test :Ch1 ∉ cols && :Ch2 ∉ cols && :Ch3 ∉ cols

    # 10) ErpData reduce path
    dat = create_test_epoch_data(n = 500)
    dat = eegfun.channel_average(dat, channel_selections = [eegfun.channels([:Ch1, :Ch2])]; reduce = true)
    @test all(propertynames(dat.data[1]) .== [:time, :sample, :condition, :condition_name, :epoch, :Ch1_Ch2])
    @test all(propertynames(dat.data[end]) .== [:time, :sample, :condition, :condition_name, :epoch, :Ch1_Ch2])

    # 11) Test default behavior (average all channels)
    dat = create_test_epoch_data(n = 500)
    eegfun.channel_average!(dat)  # Should use default channels() from second function
    @test :avg ∈ propertynames(dat.data[1])
    @test all(dat.data[1].avg .== (dat.data[1].Ch1 .+ dat.data[1].Ch2 .+ dat.data[1].Ch3) ./ 3)
    @test all(dat.data[end].avg .== (dat.data[end].Ch1 .+ dat.data[end].Ch2 .+ dat.data[end].Ch3) ./ 3)

    # 12) Test single channel selection with custom label
    dat = create_test_epoch_data(n = 500)
    eegfun.channel_average!(dat, channel_selections = [eegfun.channels([:Ch1, :Ch2])], output_labels = [:custom])
    @test :custom ∈ propertynames(dat.data[1])
    @test :custom ∈ propertynames(dat.data[end])
    @test all(dat.data[1].custom .== (dat.data[1].Ch1 .+ dat.data[1].Ch2) ./ 2)
    @test all(dat.data[end].custom .== (dat.data[end].Ch1 .+ dat.data[end].Ch2) ./ 2)

end
