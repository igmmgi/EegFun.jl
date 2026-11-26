# Tests for pre/post edge-mirroring padding and unmirroring.

using Test
using DataFrames
using EegFun

@testset "Data Mirroring" begin
    empty_lay = EegFun.Layout(DataFrame(), nothing, nothing, nothing)
    make_ep(dfs, cond = 1, fs = 100) = EegFun.EpochData("test_data", cond, "condition_$cond", dfs, empty_lay, fs, EegFun.AnalysisInfo())
    make_erp(df, fs = 100) = EegFun.ErpData("test_data", 1, "condition_1", df, empty_lay, fs, EegFun.AnalysisInfo(), 10)

    t_base = collect(-0.2:0.1:0.2)
    ns = length(t_base)

    @testset "Mutating & non-mutating roundtrip ($kind, $side)" for kind in (:epoch, :erp), side in (:pre, :post, :both)
        df1 = DataFrame(time = copy(t_base), Cz = collect(1.0:ns), Pz = collect(ns:-1.0:1.0))
        dat =
            kind === :epoch ? make_ep([df1, DataFrame(time = copy(t_base), Cz = collect(1.0:ns) .* 2, Pz = collect(ns:-1.0:1.0) .* 2)]) :
            make_erp(df1)
        get_df(x) = kind === :epoch ? x.data[1] : x.data

        # Non-mutating mirror & unmirror
        mir = EegFun.mirror(dat, side)
        @test nrow(get_df(dat)) == ns
        @test nrow(get_df(mir)) > (side === :both ? 2 * ns : ns)
        td = diff(get_df(mir).time)
        @test all(td .≈ td[1])

        unmir = EegFun.unmirror(mir, side)
        @test nrow(get_df(unmir)) == ns
        @test get_df(unmir).time ≈ t_base
        @test get_df(unmir).Cz ≈ collect(1.0:ns)

        # Mutating mirror! & unmirror! (double roundtrip consistency)
        for _ = 1:2
            EegFun.mirror!(dat, side)
            @test nrow(get_df(dat)) > (side === :both ? 2 * ns : ns)
            EegFun.unmirror!(dat, side)
            @test nrow(get_df(dat)) == ns
            @test get_df(dat).time ≈ t_base
            @test get_df(dat).Cz ≈ collect(1.0:ns)
            @test get_df(dat).Pz ≈ collect(ns:-1.0:1.0)
        end
    end

    @testset "Invalid side, multi-epoch & metadata preservation" begin
        ep_inv = make_ep([DataFrame(time = copy(t_base), Cz = collect(1.0:ns))])
        @test_throws Exception EegFun.mirror!(ep_inv, :invalid)
        @test_throws Exception EegFun.unmirror!(ep_inv, :middle)

        ep_multi = make_ep([
            DataFrame(time = copy(t_base), Cz = collect(1.0:ns)),
            DataFrame(time = copy(t_base), Cz = collect(1.0:ns) .* 2),
            DataFrame(time = copy(t_base), Cz = collect(ns:-1.0:1.0)),
        ])
        orig_czs = [copy(df.Cz) for df in ep_multi.data]
        EegFun.mirror!(ep_multi, :both)
        EegFun.unmirror!(ep_multi, :both)
        @test [df.Cz for df in ep_multi.data] ≈ orig_czs

        ep_meta = make_ep(
            [
                DataFrame(
                    time = copy(t_base),
                    Cz = collect(1.0:ns),
                    Pz = collect(ns:-1.0:1.0),
                    trial = fill(1, ns),
                    condition = fill(2, ns),
                    response = fill("left", ns),
                ),
            ],
            2,
        )
        EegFun.mirror!(ep_meta, :both)
        @test all(ep_meta.data[1].trial .== 1) && ep_meta.condition == 2 && all(ep_meta.data[1].response .== "left")
        EegFun.unmirror!(ep_meta, :both)
        @test all(ep_meta.data[1].trial .== 1) && all(ep_meta.data[1].condition .== 2) && all(ep_meta.data[1].response .== "left")
    end

    @testset "Mirroring pattern validation" begin
        epochs = make_ep([DataFrame(time = [0.0, 0.1, 0.2, 0.3, 0.4], Cz = [1.0, 2.0, 3.0, 4.0, 5.0])], 1, 10)
        for (side, expected) in (
            (:pre, [5.0, 4.0, 3.0, 2.0, 1.0, 2.0, 3.0, 4.0, 5.0]),
            (:post, [1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0]),
            (:both, [5.0, 4.0, 3.0, 2.0, 1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0]),
        )
            EegFun.mirror!(epochs, side)
            @test epochs.data[1].Cz ≈ expected
            if side === :both
                @test sum(diff(epochs.data[1].Cz) .== 0) == 0
                orig = epochs.data[1].Cz[5:9]
                @test epochs.data[1].Cz[1:4] ≈ reverse(orig[2:end])
                @test epochs.data[1].Cz[10:end] ≈ reverse(orig[1:(end-1)])
            end
            EegFun.unmirror!(epochs, side)
        end
    end
end
