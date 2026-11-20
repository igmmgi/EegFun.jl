# Tests for FieldTrip `.mat` file import.

using Test
using EegFun

@testset "FieldTrip Import" begin
    # Test file path
    test_file = joinpath(pkgdir(EegFun), "resources", "data", "fieldtrip", "continuous.mat")
    @test isfile(test_file)

    layout_file = joinpath(pkgdir(EegFun), "resources", "layouts", "biosemi", "biosemi64.csv")
    @test isfile(layout_file)
    layout = EegFun.read_layout(layout_file)

    @testset "Basic Loading Continuous" begin
        dat = EegFun.read_fieldtrip(test_file, layout)

        @test dat isa EegFun.ContinuousData
        @test size(dat.data, 1) > 0
        @test :time in propertynames(dat.data)
        @test :sample in propertynames(dat.data)

        ch_labels = EegFun.channel_labels(dat)
        @test length(ch_labels) > 0
        @test :Fp1 in ch_labels || :FP1 in ch_labels

        # Verify layout copy isolation
        dat.layout.data[1, :label] = :Mutated
        @test layout.data[1, :label] != :Mutated
    end

    @testset "Basic Loading Epochs" begin
        epoch_file = joinpath(pkgdir(EegFun), "resources", "data", "fieldtrip", "epochs.mat")
        @test isfile(epoch_file)

        dat = EegFun.read_fieldtrip(epoch_file, layout)
        @test dat isa EegFun.EpochData
        @test length(dat.data) > 0
        @test :time in propertynames(dat.data[1])
        @test :sample in propertynames(dat.data[1])
        @test :epoch in propertynames(dat.data[1])
        @test dat.data[1].epoch[1] == 1
    end

    @testset "cfg.event Dict extraction and epoched Status channel stripping" begin
        mktempdir() do tmpdir
            lay2 = EegFun.create_test_layout(n_channels = 2)
            cont_mat = joinpath(tmpdir, "ft_events.mat")
            ft_cont = Dict{String,Any}(
                "label" => Any["Ch1", "Ch2", "Status"],
                "fsample" => 100.0,
                "trial" => Any[ones(Float64, 3, 20)],
                "time" => Any[collect(0.0:0.01:0.19)],
                "cfg" => Dict{String,Any}(
                    "event" => Dict{String,Any}(
                        "sample" => Any[2.0, 10.0],
                        "value" => Any[11.0, "S 22"],
                        "type" => Any["Stimulus", "Response"],
                    ),
                ),
            )
            EegFun.MAT.matwrite(cont_mat, Dict("data" => ft_cont))
            cdat = EegFun.read_fieldtrip(cont_mat, lay2)
            @test cdat isa EegFun.ContinuousData
            @test EegFun.channel_labels(cdat) == [:Ch1, :Ch2]
            @test cdat.data.trigger[2] == 11
            @test cdat.data.trigger_info[2] == "Stimulus"
            @test cdat.data.trigger[10] == 22
            @test cdat.data.trigger_info[10] == "Response"

            ep_mat = joinpath(tmpdir, "ft_ep_status.mat")
            ft_ep = Dict{String,Any}(
                "label" => Any["Ch1", "Ch2", "Status"],
                "fsample" => 100.0,
                "trial" => Any[ones(Float64, 3, 10), fill(2.0, 3, 10)],
                "time" => Any[collect(-0.02:0.01:0.07), collect(-0.02:0.01:0.07)],
                "trialinfo" => [5.0, 6.0],
            )
            EegFun.MAT.matwrite(ep_mat, Dict("data" => ft_ep))
            edat = EegFun.read_fieldtrip(ep_mat, lay2)
            @test edat isa EegFun.EpochData
            @test EegFun.channel_labels(edat) == [:Ch1, :Ch2]
            @test :Status ∉ propertynames(edat.data[1])
        end
    end
end
