using Test
using DataFrames
using eegfun
using Statistics


@testset "correlation" begin

    @testset "correlation_matrix" begin

        dat = create_test_data(n=1000, n_channels=100)

        # Test with all channels
        cm = eegfun.correlation_matrix(dat)

        @test cm isa DataFrame
        @test :row in propertynames(cm)
        @test size(cm, 1) == 100  # 4 channels
        @test size(cm, 2) == 101  # 4 channels + row column

        # Test diagonal elements are 1.0 (correlation of each channel with itself)
        for i = 1:100
            channel_name = cm[i, :row]
            @test isapprox(cm[i, channel_name], 1.0, atol = 1e-10)
        end

        # Test symmetry
        @test isapprox(cm[1, :Ch2], cm[2, :Ch1], atol = 1e-10)
        @test isapprox(cm[1, :Ch3], cm[3, :Ch1], atol = 1e-10)

        # Test with channel selection
        cm = eegfun.correlation_matrix(dat, channel_selection = eegfun.channels([:Ch1, :Ch2]))
        @test size(cm, 1) == 2
        @test size(cm, 2) == 3  # 2 channels + row column

        # Test empty channel selection
        @test_throws ErrorException eegfun.correlation_matrix(dat, channel_selection = eegfun.channels(Symbol[]))
    end

    @testset "channel_joint_probability" begin
        dat = create_test_data(n=1000, n_channels=100)

        # Test with default parameters
        jp = eegfun.channel_joint_probability(dat)
        @test jp isa DataFrame
        @test :channel in propertynames(jp)
        @test :jp in propertynames(jp)
        @test size(jp, 1) == 100  # 4 channels

        # Test jp values are finite (can be negative due to normalization)
        for i = 1:100
            @test isfinite(jp[i, :jp])
        end

        # Test with custom parameters
        jp = eegfun.channel_joint_probability(dat, threshold = 0.3, normalize = 2, discret = 500)
        @test jp isa DataFrame
        @test size(jp, 1) == 100

        # Test with channel selection
        jp = eegfun.channel_joint_probability(dat, channel_selection = eegfun.channels([:Ch1, :Ch2]))
        @test size(jp, 1) == 2

        # Test empty channel selection
        @test_throws ErrorException eegfun.channel_joint_probability(dat, channel_selection = eegfun.channels(Symbol[]))
    end

    @testset "_correlation_matrix" begin
        dat = create_test_data(n=1000, n_channels=100)

        # Test internal function
        cm = eegfun._correlation_matrix(dat.data, collect(1:size(dat.data, 1)), [:Ch1, :Ch2, :Ch3])
        @test cm isa DataFrame
        @test :row in propertynames(cm)
        @test size(cm, 1) == 3
        @test size(cm, 2) == 4  # 3 channels + row column
    end

    @testset "_channel_joint_probability" begin
        dat = create_test_data(n=1000, n_channels=100)

        # Test internal function
        jp = eegfun._channel_joint_probability(
            dat.data,
            collect(1:size(dat.data, 1)),
            [:Ch1, :Ch2, :Ch3],
            threshold = 0.5,
            normval = 1,
            discret = 1000,
        )
        @test jp isa DataFrame
        @test :channel in propertynames(jp)
        @test :jp in propertynames(jp)
        @test size(jp, 1) == 3
    end

    @testset "_joint_probability" begin
        # Test with simple data (channels × samples)
        signal = randn(3, 100)  # 3 channels × 100 samples
        jp = eegfun._joint_probability(signal, 0.5, 1, 100)

        @test jp isa Tuple
        @test length(jp) == 2
        jp, rej = jp
        @test jp isa Vector{Float64}
        @test rej isa BitVector
        @test length(jp) == 3
        @test length(rej) == 3
    end

    @testset "compute_probability!" begin
        # Test probability computation
        data = randn(1000)
        proba_map = zeros(Float64, 1000)  # Must match data length

        result = eegfun.compute_probability!(proba_map, data, 50)

        @test result isa Vector{Float64}
        @test length(result) == 1000
        @test all(0 <= x <= 1 for x in result)  # All values between 0 and 1
    end

    @testset "_trim_extremes" begin
        # Test extreme value trimming
        x = [1.0, 2.0, 3.0, 100.0, 4.0, 5.0, 6.0, 200.0, 7.0, 8.0, 9.0]
        trimmed = eegfun._trim_extremes(x)

        @test trimmed isa SubArray{Float64,1,Vector{Float64},Tuple{UnitRange{Int64}},true}
        @test length(trimmed) < length(x)  # Should be shorter after trimming
        @test maximum(trimmed) < maximum(x)  # Should trim extreme values
        @test minimum(trimmed) >= minimum(x)  # Should trim extreme values (can be equal if no extreme low values)

        # Test with data that has both high and low extremes
        x2 = [1.0, 2.0, 3.0, 100.0, 4.0, 5.0, 6.0, 200.0, 7.0, 8.0, -50.0]
        trimmed2 = eegfun._trim_extremes(x2)
        @test maximum(trimmed2) < maximum(x2)  # Should trim high extreme values
        @test minimum(trimmed2) > minimum(x2)  # Should trim low extreme values
        @test length(trimmed2) < length(x2)  # Should be shorter after trimming
    end

end
