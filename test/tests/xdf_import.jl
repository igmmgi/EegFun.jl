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
        @test EegFun._unit_to_microvolt_scale("nV") == 1e-3

        # BrainVision V/mV unit scaling in create_eegfun_data
        bv_mod = EegFun.BrainVisionDataFormat
        bv_hdr = bv_mod.BrainVisionHeader(
            "test.eeg",
            "test.vmrk",
            "BINARY",
            "MULTIPLEXED",
            "IEEE_FLOAT_32",
            3,
            4000.0,
            250.0,
            ["ChV", "ChmV", "ChuV"],
            ["REF", "REF", "REF"],
            [1.0, 1.0, 1.0],
            ["V", "mV", "uV"],
            4,
            1,
            0,
            nothing,
        )
        bv_raw = bv_mod.BrainVisionData("subj1.vhdr", bv_hdr, ones(Float64, 4, 3), bv_mod.BrainVisionMarker[])
        bv_dat = EegFun.create_eegfun_data(bv_raw)
        @test all(==(1e6), bv_dat.data.ChV)
        @test all(==(1e3), bv_dat.data.ChmV)
        @test all(==(1.0), bv_dat.data.ChuV)

        # XDF categorical string markers receive distinct codes and out-of-range markers are ignored
        xdf_mod = EegFun.ExtensibleDataFormat
        eeg_hdr = xdf_mod.XdfStreamHeader("EEG", "EEG", 2, 100.0, 100.0, "float32", Dict{String,Any}(), "")
        eeg_str = xdf_mod.XdfStream(1, eeg_hdr, collect(0.0:0.01:0.19), zeros(Float32, 20, 2), "")
        mrk_hdr = xdf_mod.XdfStreamHeader("Markers", "Markers", 1, 0.0, 0.0, "string", Dict{String,Any}(), "")
        mrk_ts = [-5.0, 0.03, 0.08, 0.13, 0.17, 10.0]
        mrk_vals = reshape(["early_oob", "target", "standard", "deviant", "target", "late_oob"], :, 1)
        mrk_str = xdf_mod.XdfStream(2, mrk_hdr, mrk_ts, mrk_vals, "")
        mock_xdf = xdf_mod.XdfData("mock.xdf", Dict{String,Any}(), Dict{Int,xdf_mod.XdfStream}(1 => eeg_str, 2 => mrk_str))
        xdf_dat = EegFun.create_eegfun_data(mock_xdf)
        @test xdf_dat.data.trigger[1] == 0 && xdf_dat.data.trigger[end] == 0
        t_codes = [xdf_dat.data.trigger[4], xdf_dat.data.trigger[9], xdf_dat.data.trigger[14]]
        @test all(>(0), t_codes) && length(unique(t_codes)) == 3
        @test xdf_dat.data.trigger[18] == xdf_dat.data.trigger[4]
        @test EegFun.trigger_count(mock_xdf).data == EegFun.trigger_count(xdf_dat).data
    end

end # @testset "XDF Import"
