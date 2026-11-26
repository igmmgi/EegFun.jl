# Tests for average, mastoid, and custom-channel EEG rereferencing.

using Test
using DataFrames
using JLD2
using Statistics
using EegFun

@testset "rereference" begin

    dat = EegFun.create_test_continuous_data(n_channels = 2)
    # 1) Average reference (:avg) subtracts mean of all channels from each selected channel
    EegFun.rereference!(dat, :avg)
    # After average reference, the mean of all channels should be approximately zero
    mean_of_channels = mean([dat.data.Ch1, dat.data.Ch2])
    @test all(abs.(mean_of_channels) .< 1e-10)
    @test EegFun.reference(dat.analysis_info) == :avg

    # 2) Single channel reference uses Ch1 as reference
    dat = EegFun.create_test_continuous_data(n_channels = 2)
    dat_ref = EegFun.rereference(dat, :Ch1)
    # After single channel reference, the reference channel should be zero
    @test all(dat_ref.data.Ch1 .== 0.0)
    # Other channels should be original minus reference
    @test all(dat_ref.data.Ch2 .== dat.data.Ch2 .- dat.data.Ch1)
    @test EegFun.reference(dat_ref.analysis_info) == :Ch1

    # 4) EpochData: per-epoch reference computed independently
    dat = EegFun.create_test_epoch_data(n_epochs = 2, n_channels = 2, condition = 1)
    EegFun.rereference!(dat, :avg)
    # After average reference, the mean of all channels should be approximately zero for each epoch
    for epoch in dat.data
        mean_of_channels = mean([epoch.Ch1, epoch.Ch2])
        @test all(abs.(mean_of_channels) .< 1e-10)
    end

    # 5) Non-mutating version returns a new object; original unchanged
    dat = EegFun.create_test_epoch_data(n_epochs = 2, n_channels = 2, condition = 1)
    dat_ref = EegFun.rereference(dat, :Ch1)
    @test :Ch1 ∈ propertynames(dat_ref.data[1]) && :Ch1 ∈ propertynames(dat.data[1])
    @test !all(dat_ref.data[1].Ch1 .== dat.data[1].Ch1)


    # 7) Channel included in both reference and selection becomes zero
    dat = EegFun.create_test_epoch_data(n_epochs = 2, n_channels = 2, condition = 1)
    EegFun.rereference!(dat, [:Ch1]; channel_selection = EegFun.channels([:Ch1]))
    @test all(dat.data[1].Ch1 .== 0.0)

    # 9) Missing reference channel should throw
    dat = EegFun.create_test_epoch_data(n_epochs = 2, n_channels = 2, condition = 1)
    @test_throws Any EegFun.rereference!(copy(dat), [:Z]; channel_selection = EegFun.channels([:A]))

    # 10) Ergonomic channel selection and reference inputs
    dat_cont = EegFun.create_test_continuous_data(n_channels = 3)
    res_sym = EegFun.rereference(dat_cont, :Ch1; channel_selection = :Ch2)
    @test all(res_sym.data.Ch2 .== dat_cont.data.Ch2 .- dat_cont.data.Ch1)
    @test res_sym.data.Ch3 == dat_cont.data.Ch3

    res_str = EegFun.rereference(dat_cont, "Ch1"; channel_selection = ["Ch2", "Ch3"])
    @test all(res_str.data.Ch2 .== dat_cont.data.Ch2 .- dat_cont.data.Ch1)
    @test all(res_str.data.Ch3 .== dat_cont.data.Ch3 .- dat_cont.data.Ch1)

    # 11) Vector dispatch for rereference and rereference!
    erp1 = EegFun.create_test_erp_data(n_channels = 3, condition = 1)
    erp2 = EegFun.create_test_erp_data(n_channels = 3, condition = 2)
    vec_res = EegFun.rereference([erp1, erp2], :Ch1)
    @test length(vec_res) == 2
    @test all(vec_res[1].data.Ch1 .== 0.0)
    @test all(vec_res[2].data.Ch1 .== 0.0)

    EegFun.rereference!([erp1, erp2], :Ch1)
    @test all(erp1.data.Ch1 .== 0.0)
    @test all(erp2.data.Ch1 .== 0.0)
end



@testset "Batch Rereference" begin
    isdefined(@__MODULE__, :run_batch_jld2_contract_tests!) ||
        Base.include(@__MODULE__, joinpath(pkgdir(EegFun), "test", "tests", "test_helpers.jl"))

    run_batch_jld2_contract_tests!(
        (pat; kwargs...) -> EegFun.rereference(pat; reference_selection = :avg, kwargs...);
        make_data = p -> [EegFun.create_test_erp_data(participant = p, condition = c) for c in (1, 2)],
        pattern = "erps_unrejected",
        n_participants = 3,
        log_file = "rereference.log",
        log_snippets = ["Batch rereferencing started", "Reference settings: avg", "Batch operation complete"],
        expected_auto_dir = "rereferenced_erps_unrejected_avg",
    )

    mktempdir() do test_dir
        orig_erps = [EegFun.create_test_erp_data(participant = 1, condition = c) for c in (1, 2)]
        jldsave(joinpath(test_dir, "1_erps_unrejected.jld2"); data = orig_erps)

        # Reference modes (:avg, :Ch2, [:Ch1], [:Ch1, :Ch2]) & average-ref zero-mean + metadata check
        for (tag, ref_sel) in (("avg", :avg), ("mast", :Ch2), ("single", [:Ch1]), ("multi", [:Ch1, :Ch2]))
            odir = joinpath(test_dir, "ref_$tag")
            res = EegFun.rereference("erps_unrejected"; input_dir = test_dir, reference_selection = ref_sel, output_dir = odir)
            @test res.success == 1 && res.errors == 0
            if ref_sel === :avg
                reref_erp = load(joinpath(odir, "1_erps_unrejected.jld2"), "data")[1]
                @test reref_erp.sample_rate == orig_erps[1].sample_rate && reref_erp.n_epochs == orig_erps[1].n_epochs
                @test reref_erp.layout.data.label == orig_erps[1].layout.data.label
                @test reref_erp.analysis_info.reference == :avg && orig_erps[1].analysis_info.reference == :none
                @test reref_erp.data.time == orig_erps[1].data.time
                @test maximum(abs, (reref_erp.data.Ch1 .+ reref_erp.data.Ch2 .+ reref_erp.data.Ch3) ./ 3) < 1e-10
            end
        end

        # EpochData batch support
        jldsave(
            joinpath(test_dir, "1_epochs_unrejected.jld2");
            data = EegFun.create_test_epoch_data_vector(conditions = 1:2, n_channels = 3),
        )
        ep_out = joinpath(test_dir, "reref_epochs")
        @test EegFun.rereference("epochs_unrejected"; input_dir = test_dir, reference_selection = :avg, output_dir = ep_out).success == 1
        ep_loaded = load(joinpath(ep_out, "1_epochs_unrejected.jld2"), "data")
        @test length(ep_loaded) == 2 && ep_loaded[1] isa EegFun.EpochData

        # Invalid reference channel & unrecognized/empty JLD2 edge cases
        res_bad = EegFun.rereference(
            "erps_unrejected";
            input_dir = test_dir,
            reference_selection = [:InvalidChannel],
            output_dir = joinpath(test_dir, "bad_ref"),
        )
        @test res_bad.errors > 0

        edge_dir = mkpath(joinpath(test_dir, "edge"))
        jldsave(joinpath(edge_dir, "1_erps_unrejected.jld2"); data = EegFun.ErpData[])
        jldsave(joinpath(edge_dir, "2_erps_unrejected.jld2"); other_data = "test")
        res_edge = EegFun.rereference("erps_unrejected"; input_dir = edge_dir, output_dir = joinpath(edge_dir, "out"))
        @test res_edge.success == 1 && res_edge.errors == 1
    end

    @testset "In-memory :avg rereference respects channel_selection" begin
        erp_sel = EegFun.create_test_erp_data(n_channels = 3)
        erp_sel.data.Ch1 .= 10.0
        erp_sel.data.Ch2 .= 20.0
        erp_sel.data.Ch3 .= 1000.0
        EegFun.rereference!(erp_sel, :avg; channel_selection = EegFun.channels([:Ch1, :Ch2]))
        @test all(erp_sel.data.Ch1 .≈ -5.0)
        @test all(erp_sel.data.Ch2 .≈ 5.0)
        @test all(erp_sel.data.Ch3 .≈ 1000.0)
    end
end
