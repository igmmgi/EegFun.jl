# Tests for EEG DataFrame accessors, channel/metadata column queries, and data helpers.

using Test
using Statistics
using EegFun
using DataFrames
using OrderedCollections

@testset "Data Utilities" begin
    make_lay(chs::Vector{Symbol}) =
        EegFun.Layout(DataFrame(label = chs, inc = zeros(length(chs)), azi = zeros(length(chs))), nothing, nothing, nothing)
    base_lay2 = make_lay([:Fz, :Cz])
    base_df3 = DataFrame(time = [0.1, 0.2, 0.3], sample = [1, 2, 3], Fz = [1.0, 2.0, 3.0], Cz = [4.0, 5.0, 6.0])
    make_cont3(; info = EegFun.AnalysisInfo()) = EegFun.ContinuousData("test_data", copy(base_df3), copy(base_lay2), 250, info)
    make_ep2() = EegFun.EpochData(
        "test_data",
        1,
        "condition_1",
        [
            DataFrame(time = [0.1, 0.2], sample = [1, 2], Fz = [1.0, 2.0], Cz = [3.0, 4.0]),
            DataFrame(time = [0.1, 0.2], sample = [1, 2], Fz = [5.0, 6.0], Cz = [7.0, 8.0]),
        ],
        copy(base_lay2),
        250,
        EegFun.AnalysisInfo(),
    )

    @testset "Column identification, data access & MultiDataFrameEeg" begin
        df_ext = DataFrame(
            time = [0.1, 0.2, 0.3],
            sample = [1, 2, 3],
            Fz = [1.0, 2.0, 3.0],
            Cz = [4.0, 5.0, 6.0],
            Pz = [7.0, 8.0, 9.0],
            vEOG = [0.1, 0.2, 0.3],
            hEOG = [0.4, 0.5, 0.6],
        )
        cd_ext = EegFun.ContinuousData("test_data", df_ext, make_lay([:Fz, :Cz, :Pz]), 250, EegFun.AnalysisInfo())
        @test EegFun._get_cols_by_group(cd_ext, :channels) == [:Fz, :Cz, :Pz]
        @test EegFun._get_cols_by_group(cd_ext, :metadata) == [:time, :sample]
        @test EegFun._get_cols_by_group(cd_ext, :extra) == [:vEOG, :hEOG]
        @test_throws Exception EegFun._get_cols_by_group(cd_ext, :invalid)

        empty_cd = EegFun.ContinuousData("test_data", DataFrame(time = [0.1], sample = [1]), EegFun.Layout(), 250, EegFun.AnalysisInfo())
        @test (
            EegFun._get_cols_by_group(empty_cd, :channels),
            EegFun._get_cols_by_group(empty_cd, :metadata),
            EegFun._get_cols_by_group(empty_cd, :extra),
        ) == (Symbol[], [:time, :sample], Symbol[])

        cd = make_cont3()
        @test EegFun.all_data(cd) == base_df3
        @test EegFun.all_labels(cd) == [:time, :sample, :Fz, :Cz] == EegFun.all_labels(base_df3)
        @test EegFun.meta_labels(cd) == [:time, :sample] && EegFun.meta_data(cd) == base_df3[:, [:time, :sample]]
        @test EegFun.channel_labels(cd) == [:Fz, :Cz] && EegFun.channel_data(cd) == base_df3[:, [:Fz, :Cz]]
        @test isempty(EegFun.extra_labels(cd)) && isempty(EegFun.extra_data(cd))

        ep = make_ep2()
        all_df = EegFun.all_data(ep)
        @test size(all_df) == (4, 4)
        @test EegFun.all_labels(ep) == [:time, :sample, :Fz, :Cz] == EegFun.all_labels(ep, 1)
        @test EegFun.meta_data(ep) == all_df[:, [:time, :sample]] && EegFun.meta_data(ep, 1) == ep.data[1][:, [:time, :sample]]
        @test EegFun.channel_data(ep) == all_df[:, [:Fz, :Cz]] && EegFun.channel_data(ep, 1) == ep.data[1][:, [:Fz, :Cz]]
        @test isempty(EegFun.extra_data(ep))
        @test EegFun.n_samples(ep) == 2 == EegFun.n_samples(ep, 1) && EegFun.n_epochs(ep) == 2
        @test EegFun.duration(ep) == 0.1 == EegFun.duration(ep, 1)
    end

    @testset "Convenience, viewing, DataFrame conversion & math/rename helpers" begin
        info = EegFun.AnalysisInfo(reference = :avg, hp_filter = 0.1, lp_filter = 30.0)
        cd = make_cont3(info = info)
        @test EegFun.sample_rate(cd) == 250 && EegFun.sample_rate(base_df3) == 10
        @test EegFun.reference(cd) == :avg == EegFun.reference(info) && EegFun.filter_info(info) == [0.1, 30.0]
        @test EegFun.n_samples(cd) == 3 == EegFun.n_samples(base_df3) && EegFun.n_channels(cd) == 2 && EegFun.n_layout(base_lay2) == 2
        @test EegFun.n_epochs(cd) == 1 && isapprox(EegFun.duration(cd), 0.2; atol = 1e-10)
        @test EegFun.has_channels(cd, [:Fz, :Cz]) && !EegFun.has_channels(cd, [:Fz, :Pz])
        cd2 = EegFun.ContinuousData(
            "test_data",
            DataFrame(time = [0.1, 0.2], sample = [1, 2], Fz = [1.0, 2.0], Pz = [7.0, 8.0]),
            make_lay([:Fz, :Pz]),
            250,
            info,
        )
        @test EegFun.common_channels(cd, cd2) == [:Fz]
        @test EegFun.n_epochs(
            EegFun.ErpData("test_data", 1, "condition_1", copy(base_df3), copy(base_lay2), 250, EegFun.AnalysisInfo(), 10),
        ) == 10

        df5 = DataFrame(time = collect(0.1:0.1:0.5), sample = 1:5, Fz = collect(1.0:5.0), Cz = collect(6.0:10.0))
        cd5 = EegFun.ContinuousData("test_data", df5, copy(base_lay2), 250, EegFun.AnalysisInfo())
        @test_nowarn EegFun.viewer(cd5)
        @test_nowarn EegFun.viewer(df5)
        @test EegFun.head(cd5; n = 3) == df5[1:3, :] && EegFun.tail(cd5; n = 3) == df5[3:5, :]
        @test EegFun.head(cd5; n = 10) == df5 == EegFun.tail(cd5; n = 10)
        empty_cd = EegFun.ContinuousData("test_data", DataFrame(), EegFun.Layout(), 250, EegFun.AnalysisInfo())
        @test isempty(EegFun.head(empty_cd)) && isempty(EegFun.tail(empty_cd))

        ep = make_ep2()
        @test size(EegFun.to_data_frame(ep)) == (4, 4) && nrow(EegFun.to_data_frame([ep, copy(ep)])) == 8
        @test isempty(
            EegFun.to_data_frame(EegFun.EpochData("test_data", 1, "condition_1", DataFrame[], copy(base_lay2), 250, EegFun.AnalysisInfo())),
        )
        @test isempty(EegFun.to_data_frame(EegFun.EpochData[]))

        df_abc = DataFrame(A = [1.0, 2.0, 3.0], B = [4.0, 5.0, 6.0], C = [7.0, 8.0, 9.0])
        @test EegFun._colmeans(df_abc, [:A, :B]) == [2.5, 3.5, 4.5] == EegFun._colmeans(df_abc, [1, 2])
        mat = [1.0 2.0; 3.0 4.0; 5.0 6.0]
        @test EegFun._colmeans(mat) == [1.5, 3.5, 5.5] == EegFun._colmeans(mat, [1, 2])

        @test EegFun._subset_dataframe(df5, [:time, :Fz], [1, 3, 5]) == df5[[1, 3, 5], [:time, :Fz]]
        @test EegFun._subset_dataframes([df5, df5], [1, 2], [:time, :Fz], [1, 3, 5]) ==
              [df5[[1, 3, 5], [:time, :Fz]], df5[[1, 3, 5], [:time, :Fz]]]

        single_ep = EegFun.epoch_to_continuous(ep, 1)
        @test single_ep isa EegFun.ContinuousData &&
              single_ep.data == ep.data[1] &&
              EegFun.epoch_to_continuous(ep, 2) isa EegFun.ContinuousData

        cd_ren = make_cont3()
        ren_copy = EegFun.rename_channel(cd_ren, Dict(:Fz => :Fz_new, :Cz => :Cz_new))
        @test :Fz in propertynames(cd_ren.data) && :Fz_new in propertynames(ren_copy.data) && :Fz ∉ propertynames(ren_copy.data)
        EegFun.rename_channel!(cd_ren, Dict(:Fz => :Fz_new, :Cz => :Cz_new))
        @test cd_ren.layout.data.label == [:Fz_new, :Cz_new] && :Fz_new in propertynames(cd_ren.data)
        @test isnothing(EegFun.log_pretty_table(DataFrame(A = [1, 2, 3], B = [4, 5, 6])))
    end

    @testset "Channel/sample selection, predicates, subset & ylimits" begin
        df_v = hcat(base_df3, DataFrame(vEOG = [0.1, 0.2, 0.3]))
        cd_v = EegFun.ContinuousData("test_data", df_v, copy(base_lay2), 250, EegFun.AnalysisInfo())
        @test EegFun.get_selected_channels(cd_v, EegFun.channels([:Fz])) == [:time, :sample, :Fz]
        @test EegFun.get_selected_channels(cd_v, EegFun.channels([:Fz]); include_meta = false) == [:Fz]
        @test EegFun.get_selected_channels(cd_v, EegFun.channels([:vEOG]); include_extra = true) == [:time, :sample, :vEOG]
        for (pred, expected) in (
            (EegFun.channels([:Cz, :Fz]), [:Cz, :Fz]),
            (EegFun.channels_not([:Cz]), [:Fz]),
            (EegFun.channels([2, 1]), [:Cz, :Fz]),
            (EegFun.channels_not([2]), [:Fz]),
            (EegFun.channels(1:2), [:Fz, :Cz]),
            (EegFun.channels_not(1:1), [:Cz]),
            (EegFun.channels(1), [:Fz]),
            (EegFun.channels(:Fz), [:Fz]),
            (EegFun.channels(), [:Fz, :Cz]),
            (EegFun.channels([:Fz, :NonExistent]), [:Fz]),
            (EegFun.channels([:Fz, :Cz, :Fz]), [:Fz, :Cz]),
            (EegFun.channels([1, 100]), [:Fz]),
            (EegFun.channels_not([:NonExistent]), [:Fz, :Cz]),
            (EegFun.channels_not([:Fz, :Fz]), [:Cz]),
            (EegFun.channels_not(1:2), Symbol[]),
        )
            @test EegFun.get_selected_channels(cd_v, pred; include_meta = false, include_extra = false) == expected
        end
        df_bool = DataFrame(time = [0.1, 0.2, 0.3], sample = [1, 2, 3], Fz = [1.0, 2.0, 3.0], flag = [true, false, true])
        @test EegFun.get_selected_samples(df_bool, EegFun.samples(:flag)) == [1, 3]
        @test EegFun.samples()(df_bool) == [true, true, true]
        @test EegFun.samples(:flag)(df_bool) == [true, false, true] && EegFun.samples_not(:flag)(df_bool) == [false, true, false]

        chs = [:Fz, :Cz, :Pz]
        @test EegFun.channels()(chs) == [true, true, true]
        @test EegFun.channels([:Fz, :Cz])(chs) == [true, true, false] == EegFun.channels([1, 2])(chs) == EegFun.channels(1:2)(chs)
        @test EegFun.channels(:Fz)(chs) == [true, false, false]
        @test EegFun.channels_not([:Fz])(chs) == [false, true, true] == EegFun.channels_not(:Fz)(chs)
        @test EegFun.channels_not([1, 2])(chs) == [false, false, true]

        for (pos_fn, neg_fn) in ((EegFun.components, EegFun.components_not), (EegFun.epochs, EegFun.epochs_not))
            ids = collect(1:5)
            @test pos_fn()(ids) == trues(5) && pos_fn([1, 3, 5])(ids) == [true, false, true, false, true]
            @test pos_fn(1:3)(ids) == [true, true, true, false, false] && pos_fn(2)(ids) == [false, true, false, false, false]
            @test neg_fn([1, 3, 5])(ids) == [false, true, false, true, false] && neg_fn(2)(ids) == [true, false, true, true, true]
        end

        dvec = [EegFun.create_test_erp_data(condition = c) for c = 1:3]
        @test EegFun.conditions("condition_1")(dvec) == [1] && EegFun.conditions(["condition_3", "condition_1"])(dvec) == [3, 1]
        @test EegFun.conditions_not("condition_1")(dvec) == [false, true, true]
        @test EegFun.conditions_not(["condition_1", "condition_2"])(dvec) == [false, false, true]
        @test EegFun.conditions_not("condition_1", "condition_3")(dvec) == [false, true, false]

        sub_cd = EegFun.subset(cd_v; channel_selection = EegFun.channels([:Fz]), sample_selection = EegFun.samples(), include_extra = true)
        @test sub_cd isa EegFun.ContinuousData && nrow(sub_cd.data) == 3 && :Fz in propertynames(sub_cd.data)
        erp_d = EegFun.ErpData("test_data", 1, "condition_1", copy(base_df3), copy(base_lay2), 250, EegFun.AnalysisInfo(), 5)
        @test EegFun.subset(erp_d; channel_selection = EegFun.channels([:Fz])).n_epochs == 5
        ep = make_ep2()
        @test length(EegFun.subset(ep; channel_selection = EegFun.channels([:Fz]), epoch_selection = EegFun.epochs([1])).data) == 1

        ylims = EegFun._ylimits(erp_d)
        @test ylims isa Tuple{Float64,Float64} && ylims[1] < ylims[2] && EegFun._ylimits(make_cont3()) == ylims
        yl_ep = EegFun._ylimits(ep)
        @test yl_ep isa Tuple{Float64,Float64} && yl_ep[1] < yl_ep[2]
    end

    @testset "Base Array Interface and Conversions" begin
        df = DataFrame(time = [0.1, 0.2, 0.3], sample = [1, 2, 3], Fz = [1.0, 2.0, 3.0], Cz = [4.0, 5.0, 6.0])
        layout = EegFun.Layout(DataFrame(label = [:Fz, :Cz], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing, nothing)
        analysis_info = EegFun.AnalysisInfo()
        cont = EegFun.ContinuousData("test", df, layout, 250, analysis_info)

        epoch1 = DataFrame(time = [0.1, 0.2], sample = [1, 2], Fz = [1.0, 2.0], Cz = [3.0, 4.0])
        epoch2 = DataFrame(time = [0.1, 0.2], sample = [1, 2], Fz = [5.0, 6.0], Cz = [7.0, 8.0])
        ep = EegFun.EpochData("test", 1, "cond", [epoch1, epoch2], layout, 250, analysis_info)

        # ndims, isempty
        @test ndims(cont) == 2
        @test ndims(ep) == 3
        @test !isempty(cont)
        @test !isempty(ep)

        # EegData is not an AbstractArray: indexing returns subsets, so no scalar eltype / linear indexing
        @test !hasmethod(lastindex, Tuple{EegFun.ContinuousData})
        @test !hasmethod(firstindex, Tuple{EegFun.ContinuousData})
        @test eltype(ep) == DataFrame          # what iterating epoch data yields
        @test length(ep) == 2
        @test collect(ep) == [epoch1, epoch2]
        @test all(x -> x isa DataFrame, ep)

        # firstindex, lastindex (per dimension), axes
        @test firstindex(cont, 1) == 1
        @test firstindex(cont, 2) == 1
        @test lastindex(cont, 1) == 3
        @test lastindex(cont, 2) == 2
        @test size(cont[2:end, :]) == (2, 2)
        @test axes(cont, 1) == Base.OneTo(3)
        @test axes(cont, 2) == Base.OneTo(2)
        @test axes(ep, 1) == Base.OneTo(2)

        # begin indexing
        @test size(cont[begin:2, :]) == (2, 2)
        @test size(ep[begin:1, :, :]) == (1, 2, 2)

        # Shorthand channel subsetting
        @test EegFun.channel_labels(cont[:Cz]) == [:Cz]
        @test EegFun.channel_labels(cont[[:Fz, :Cz]]) == [:Fz, :Cz]
        @test EegFun.channel_labels(cont["Cz"]) == [:Cz]
        @test EegFun.channel_labels(cont[["Fz", "Cz"]]) == [:Fz, :Cz]
        @test EegFun.channel_labels(ep[:Cz]) == [:Cz]

        # CartesianIndex getindex & setindex!
        @test cont[CartesianIndex(2, 2)] isa EegFun.ContinuousData
        cont_copy = copy(cont)
        cont_copy[CartesianIndex(1, 1)] = 99.0
        @test cont_copy.data[1, :Fz] == 99.0

        # Tensor conversions & vec
        @test size(Array(cont)) == (3, 2)
        @test size(Array(ep)) == (2, 2, 2)
        # No implicit conversions (would silently drop layout/metadata); explicit constructors only
        @test_throws MethodError convert(Matrix, cont)
        @test_throws MethodError convert(Array, ep)
        @test_throws MethodError convert(DataFrame, cont)
        @test Matrix{Float64}(cont) == Matrix(cont)
        @test Array{Float64}(cont) == Array(cont)
        @test Array{Float64}(ep) == Array(ep)
        @test Array{Float64,2}(cont) == Matrix(cont)
        @test Array{Float64,3}(ep) == Array(ep)
        @test length(vec(cont)) == 6
        @test length(vec(ep)) == 8

        # Channel membership & keys
        @test :Cz in cont
        @test !(:Oz in cont)
        @test "Cz" in cont
        @test haskey(cont, :Cz)
        @test !haskey(cont, :Oz)
        @test keys(cont) == [:Fz, :Cz]
        @test :Cz in ep
        @test keys(ep) == [:Fz, :Cz]

        # Equality & isapprox
        @test cont == copy(cont)
        @test isequal(cont, copy(cont))
        @test ep == copy(ep)
        @test isequal(ep, copy(ep))
        @test cont ≈ copy(cont)
        @test ep ≈ copy(ep)

        cont_diff = copy(cont)
        cont_diff.data[1, :Fz] += 1e-8
        @test cont != cont_diff
        @test isapprox(cont, cont_diff; atol = 1e-6)
        @test !isapprox(cont, ep)

        # isapprox also requires matching sample rates
        cont_sr = copy(cont)
        cont_sr.sample_rate = 500
        @test !isapprox(cont, cont_sr)

        # hash is consistent with isequal (usable as Dict keys / in Sets)
        @test hash(cont) == hash(copy(cont))
        @test hash(ep) == hash(copy(ep))
        @test hash(cont) != hash(cont_diff)
        @test length(Set([cont, copy(cont), cont_diff])) == 2
        @test Dict(cont => 1)[copy(cont)] == 1
        @test hash(layout) == hash(copy(layout))
        @test hash(analysis_info) == hash(copy(analysis_info))

        # NaN: `==` follows IEEE (NaN != NaN), `isequal` treats NaNs as equal
        cont_nan = copy(cont)
        cont_nan.data[1, :Fz] = NaN
        @test cont_nan != copy(cont_nan)
        @test isequal(cont_nan, copy(cont_nan))
        @test hash(cont_nan) == hash(copy(cont_nan))

        # Neighbors equality by value
        nb1 = EegFun.Neighbors([:Cz], [1.0], [1.0])
        nb2 = EegFun.Neighbors([:Cz], [1.0], [1.0])
        @test nb1 == nb2
        @test isequal(nb1, nb2)
        @test hash(nb1) == hash(nb2)
        @test nb1 != EegFun.Neighbors([:Pz], [1.0], [1.0])

        # summary
        @test summary(cont) == "3 timepoints × 2 channels ContinuousData (\"test\")"
        @test summary(ep) == "2 epochs × 2 timepoints × 2 channels EpochData (\"test\")"
    end

    @testset "Polymorphic Accessors (Layout, TF, Stats)" begin
        layout_df = DataFrame(label = [:Fz, :Cz, :Pz], x = [0.0, 0.0, 0.0], y = [1.0, 0.0, -1.0], z = [0.0, 1.0, 0.0])
        layout = EegFun.Layout(layout_df, nothing, nothing, nothing)
        @test EegFun.channel_labels(layout) == [:Fz, :Cz, :Pz]
        @test EegFun.n_channels(layout) == 3

        # TF Data frequency_vector
        tf_df =
            DataFrame(freq = [4.0, 4.0, 8.0, 8.0, 12.0, 12.0], time = [0.0, 0.1, 0.0, 0.1, 0.0, 0.1], Fz = fill(1.0, 6), Cz = fill(2.0, 6))
        tf_data = EegFun.TimeFreqData("test_tf", 1, "Cond1", tf_df, copy(tf_df), layout, 100, :wavelet, nothing, EegFun.AnalysisInfo())
        @test EegFun.frequency_vector(tf_data) == [4.0, 8.0, 12.0]
        @test EegFun.channel_labels(tf_data) == [:Fz, :Cz]
        @test EegFun.n_channels(tf_data) == 2

        tf_epoch = EegFun.TimeFreqEpochData(
            "test_tf_ep",
            1,
            "Cond1",
            [tf_df, copy(tf_df)],
            [copy(tf_df), copy(tf_df)],
            layout,
            100,
            :wavelet,
            nothing,
            EegFun.AnalysisInfo(),
        )
        @test EegFun.frequency_vector(tf_epoch) == [4.0, 8.0, 12.0]

        spec_df = DataFrame(freq = [2.0, 4.0, 6.0], Fz = [1.0, 2.0, 3.0], Cz = [2.0, 3.0, 4.0])
        spec_data = EegFun.SpectrumData("test_spec", 1, "Cond1", spec_df, layout, 100, :welch, EegFun.AnalysisInfo())
        @test EegFun.frequency_vector(spec_data) == [2.0, 4.0, 6.0]

        # TFStatisticalData
        tf_analysis = EegFun.TFAnalysisData(:paired, [zeros(2, 2, 3, 2), zeros(2, 2, 3, 2)], [4.0, 8.0, 12.0], [0.0, 0.1])
        tf_stat = EegFun.TFStatisticalData([tf_data, tf_data], tf_analysis)
        @test EegFun.time_vector(tf_stat) == [0.0, 0.1]
        @test EegFun.frequency_vector(tf_stat) == [4.0, 8.0, 12.0]
        @test EegFun.channel_labels(tf_stat) == [:Fz, :Cz]
        @test EegFun.n_channels(tf_stat) == 2
        @test EegFun.n_samples(tf_stat) == 2
    end

    @testset "Zero-epoch guards, vcat validation, TF rename, circular phase, and NaN-safe ylimits" begin
        layout = EegFun.Layout(DataFrame(label = [:Fz, :Cz], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing, nothing)
        info = EegFun.AnalysisInfo()
        empty_ep = EegFun.EpochData("empty", 1, "cond_1", DataFrame[], layout, 250, info)

        # 0-epoch and empty collection guards
        @test EegFun.all_labels(empty_ep) == Symbol[]
        @test EegFun.time_vector(empty_ep) == Float64[]
        @test EegFun.n_samples(empty_ep) == 0
        @test EegFun.duration(empty_ep) == 0.0
        @test EegFun.channel_data(empty_ep, :Fz) == Float64[]
        @test isempty(EegFun.meta_data(empty_ep))
        @test size(Array(empty_ep)) == (0, 0, 2)
        @test lastindex(empty_ep) == 0
        @test EegFun.channel_labels(EegFun.EegData[]) == Symbol[]
        @test EegFun.extra_labels(EegFun.EegData[]) == Symbol[]

        # lastindex & vcat validation on EpochData
        ep1 = EegFun.EpochData(
            "f1",
            1,
            "c1",
            [DataFrame(time = [0.0, 0.1], sample = [1, 2], epoch = [1, 1], Fz = [1.0, 2.0], Cz = [3.0, 4.0])],
            layout,
            250,
            info,
        )
        @test lastindex(ep1) == 1
        ep_diff_ch = EegFun.EpochData(
            "f2",
            1,
            "c1",
            [DataFrame(time = [0.0, 0.1], sample = [1, 2], epoch = [1, 1], Fz = [1.0, 2.0])],
            EegFun.Layout(DataFrame(label = [:Fz], inc = [0.0], azi = [0.0]), nothing, nothing, nothing),
            250,
            info,
        )
        ep_diff_samp = EegFun.EpochData(
            "f3",
            1,
            "c1",
            [DataFrame(time = [0.0, 0.1, 0.2], sample = [1, 2, 3], epoch = [1, 1, 1], Fz = [1.0, 2.0, 3.0], Cz = [4.0, 5.0, 6.0])],
            layout,
            250,
            info,
        )
        @test_throws EegFun.EegFunError vcat(ep1, ep_diff_ch)
        @test_throws EegFun.EegFunError vcat(ep1, ep_diff_samp)

        # _have_same_structure on ContinuousData (no .condition field) and isapprox time vector
        c1 = EegFun.ContinuousData("c1", DataFrame(time = [0.0, 0.1], Fz = [1.0, 2.0], Cz = [3.0, 4.0]), copy(layout), 10, info)
        c2 = EegFun.ContinuousData("c2", DataFrame(time = [0.0, 0.1 + 1e-12], Fz = [5.0, 6.0], Cz = [7.0, 8.0]), copy(layout), 10, info)
        @test EegFun._have_same_structure(c1, c2)

        # rename_channel! on TimeFreqData and TimeFreqEpochData renames data_phase as well as data_power
        tf = EegFun.create_test_tf_data(n_channels = 2)
        EegFun.rename_channel!(tf, Dict(:Ch1 => :Frontal))
        @test :Frontal in propertynames(tf.data_power)
        @test :Frontal in propertynames(tf.data_phase)
        @test :Ch1 ∉ propertynames(tf.data_phase)

        tf_ep = EegFun.create_test_tf_epoch_data(n_channels = 2, n_epochs = 2)
        EegFun.rename_channel!(tf_ep, Dict(:Ch1 => :Frontal))
        @test all(df -> :Frontal in propertynames(df) && :Ch1 ∉ propertynames(df), tf_ep.data_power)
        @test all(df -> :Frontal in propertynames(df) && :Ch1 ∉ propertynames(df), tf_ep.data_phase)

        # _condition_channels excludes boolean :extra columns; circular phase difference & mean
        c_flag = copy(c1)
        c_flag.data[!, :is_extreme] = [true, false]
        @test EegFun._condition_channels(c_flag) == [:Fz, :Cz]

        tf1 = EegFun.create_test_tf_data(condition = 1, n_channels = 1)
        tf2 = EegFun.create_test_tf_data(condition = 2, n_channels = 1)
        tf1.data_phase[!, :Ch1] .= π - 0.1
        tf2.data_phase[!, :Ch1] .= -π + 0.1
        derived_diff = EegFun._derive_condition(
            tf1,
            1,
            "difference_1_2",
            (
                tf1.data_power .- tf2.data_power,
                DataFrame(time = tf1.data_phase.time, freq = tf1.data_phase.freq, Ch1 = tf1.data_phase.Ch1 .- tf2.data_phase.Ch1),
            ),
        )
        @test all(isapprox.(derived_diff.data_phase[!, :Ch1], -0.2, atol = 1e-10))

        # _find_closest_time_index empty guard & NaN-safe ylimits + include_extra
        @test_throws EegFun.EegFunError EegFun._find_closest_time_index(Float64[], 0.0)
        df_nan = DataFrame(
            time = [0.0, 0.1, 0.2],
            Fz = [NaN, -2.0, 4.0],
            Cz = [NaN, NaN, NaN],
            vEOG = [-50.0, 0.0, 50.0],
            flag = [true, false, true],
        )
        cont_nan = EegFun.ContinuousData("nan", df_nan, copy(layout), 10, info)
        ylims_no_extra = EegFun._ylimits(cont_nan; include_extra = false)
        ylims_extra = EegFun._ylimits(cont_nan; include_extra = true)
        @test ylims_extra[2] > ylims_no_extra[2]
    end
end
