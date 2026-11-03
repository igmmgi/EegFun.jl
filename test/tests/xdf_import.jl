using Test
using EegFun
using DataFrames

@testset "XDF Import" begin

    # ────────────────────────────────────────────────────────────
    # 1. read_xdf function exists and has correct signature
    # ────────────────────────────────────────────────────────────
    @testset "read_xdf API" begin
        @test hasmethod(EegFun.read_xdf, Tuple{String})
        @test_throws Exception EegFun.read_xdf("/nonexistent/file.xdf")
    end

    # ────────────────────────────────────────────────────────────
    # 2. File extension validation via read_raw_data dispatcher
    # ────────────────────────────────────────────────────────────
    @testset "read_raw_data dispatches XDF and .xdf.gz" begin
        @test_throws Exception EegFun.read_raw_data("/nonexistent/file.xdf")
        @test_throws Exception EegFun.read_raw_data("/nonexistent/file.xdf.gz")
    end

    # ────────────────────────────────────────────────────────────
    # 3. Real XDF file import, channel metadata, scaling, and effective_srate fallback
    # ────────────────────────────────────────────────────────────
    @testset "Real test1.xdf loading and unit scaling" begin
        xdf_file = joinpath(dirname(dirname(@__DIR__)), "resources", "data", "xdf", "test1.xdf")
        @test isfile(xdf_file)

        raw = EegFun.read_xdf(xdf_file)
        @test raw isa EegFun.ExtensibleDataFormat.XdfData

        dat = EegFun.create_eegfun_data(raw)
        @test dat isa EegFun.ContinuousData
        @test EegFun.sample_rate(dat) == 250
        @test EegFun.n_channels(dat) == 16
        @test :Fp1 in EegFun.channel_labels(dat)
        @test :Oz in EegFun.channel_labels(dat)
        @test count(!=(0), dat.data.trigger) > 0

        # Also test create_eegfun_data(::XdfData, ::Layout) & copy isolation
        lay = copy(dat.layout)
        dat_with_lay = EegFun.create_eegfun_data(raw, lay)
        @test EegFun.sample_rate(dat_with_lay) == 250
        dat_with_lay.layout.data[1, :label] = :Mutated
        @test lay.data[1, :label] == :Fp1

        # Test XDF unit scale helper & effective_srate fallback helper
        @test EegFun._unit_to_microvolt_scale("V") == 1e6
        @test EegFun._unit_to_microvolt_scale("volts") == 1e6
        @test EegFun._unit_to_microvolt_scale("mV") == 1e3
        @test EegFun._unit_to_microvolt_scale("millivolts") == 1e3
        @test EegFun._unit_to_microvolt_scale("microvolts") == 1.0
        @test EegFun._unit_to_microvolt_scale("uV") == 1.0

        mock_stream = (header = (nominal_srate = 0.0, effective_srate = 512.3),)
        @test EegFun._xdf_srate(mock_stream) == 512.3
    end

end # @testset "XDF Import"
