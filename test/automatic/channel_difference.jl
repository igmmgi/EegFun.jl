using Test
using DataFrames
using eegfun

@testset "channel_difference" begin

    # 1) Basic difference A - B
    dat = create_test_data(n = 500)
    eegfun.channel_difference!(
        dat;
        channel_selection1 = eegfun.channels([:Ch1]),
        channel_selection2 = eegfun.channels([:Ch2]),
        channel_out = :Ch1_minus_Ch2,
    )
    @test :Ch1_minus_Ch2 ∈ propertynames(dat.data)
    @test all(dat.data.Ch1_minus_Ch2 .== (dat.data.Ch1 .- dat.data.Ch2))

    # 2) Group average difference mean(A,B) - C
    dat = create_test_data(n = 500)
    eegfun.channel_difference!(
        dat;
        channel_selection1 = eegfun.channels([:Ch1, :Ch2]),
        channel_selection2 = eegfun.channels([:Ch3]),
        channel_out = :Ch1_Ch2_minus_Ch3,
    )
    @test :Ch1_Ch2_minus_Ch3 ∈ propertynames(dat.data)
    @test all(dat.data.Ch1_Ch2_minus_Ch3 .== ((dat.data.Ch1 .+ dat.data.Ch2) ./ 2 .- dat.data.Ch3))

    # 3) Predicate selection (same as explicit)
    dat = create_test_data(n = 500)
    eegfun.channel_difference!(
        dat;
        channel_selection1 = eegfun.channels([:Ch1, :Ch2]),
        channel_selection2 = eegfun.channels([:Ch3]),
        channel_out = :out,
    )
    @test all(dat.data.out .== ((dat.data.Ch1 .+ dat.data.Ch2) ./ 2 .- dat.data.Ch3))

    # 4) Non-mutating version returns a new object; original unchanged
    dat = create_test_data(n = 500)
    dat = eegfun.channel_difference(
        dat;
        channel_selection1 = eegfun.channels([:Ch1]),
        channel_selection2 = eegfun.channels([:Ch3]),
        channel_out = :out,
    )
    @test :out ∈ propertynames(dat.data)
    @test :xxx ∉ propertynames(dat.data)
    @test all(dat.data.out .== (dat.data.Ch1 .- dat.data.Ch3))

    # 5) Overwrite behavior: write then overwrite with different groups
    dat = create_test_data(n = 500)
    dat.data[!, :X] = zeros(500)
    eegfun.channel_difference!(
        dat;
        channel_selection1 = eegfun.channels([:Ch2]),
        channel_selection2 = eegfun.channels([:Ch1]),
        channel_out = :X,
    )
    @test all(dat.data.X .== (dat.data.Ch2 .- dat.data.Ch1))

    # 6) EpochData: append to each epoch
    dat = create_test_epoch_data(n = 500)
    eegfun.channel_difference!(
        dat;
        channel_selection1 = eegfun.channels([:Ch1]),
        channel_selection2 = eegfun.channels([:Ch2]),
    )
    @test :diff ∈ propertynames(dat.data[1]) && :diff ∈ propertynames(dat.data[2])
    @test all(dat.data[1].diff .== (dat.data[1].Ch1 .- dat.data[1].Ch2))
    @test all(dat.data[2].diff .== (dat.data[2].Ch1 .- dat.data[2].Ch2))

    # 7) ErpData (SingleDataFrameEeg): append
    dat = create_test_epoch_data(n = 500)
    eegfun.channel_difference!(
        dat;
        channel_selection1 = eegfun.channels([:Ch3]),
        channel_selection2 = eegfun.channels([:Ch2]),
        channel_out = :Ch3_minus_Ch2,
    )
    @test :Ch3_minus_Ch2 ∈ propertynames(dat.data[1])
    @test all(dat.data[1].Ch3_minus_Ch2 .== (dat.data[1].Ch3 .- dat.data[1].Ch2))

    # 10) Commutativity sanity: A-B == -(B-A)
    dat = create_test_epoch_data(n = 500)
    dat_Ch1_Ch2 = copy(dat)
    eegfun.channel_difference!(
        dat_Ch1_Ch2;
        channel_selection1 = eegfun.channels([:Ch1]),
        channel_selection2 = eegfun.channels([:Ch2]),
        channel_out = :Ch1_minus_Ch2,
    )
    dat_Ch2_Ch1 = copy(dat)
    eegfun.channel_difference!(
        dat_Ch2_Ch1;
        channel_selection1 = eegfun.channels([:Ch2]),
        channel_selection2 = eegfun.channels([:Ch1]),
        channel_out = :Ch2_minus_Ch1,
    )
    @test all(dat_Ch1_Ch2.data[1].Ch1_minus_Ch2 .== .-(dat_Ch2_Ch1.data[1].Ch2_minus_Ch1))
    @test all(dat_Ch1_Ch2.data[2].Ch1_minus_Ch2 .== .-(dat_Ch2_Ch1.data[2].Ch2_minus_Ch1))
    @test all(dat_Ch1_Ch2.data[end].Ch1_minus_Ch2 .== .-(dat_Ch2_Ch1.data[end].Ch2_minus_Ch1))

end
