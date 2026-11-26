# Tests for continuous and epoch artifact detection, flagging, counting, and rejection.

using Test
using DataFrames
using JLD2
using OrderedCollections
using Random
using EegFun

@testset "artifact_detection" begin
    function _test_detector_contract!(
        is_fn,
        is_mut_fn!,
        n_fn;
        prefix::String,
        count_col::Symbol,
        thresh::Real,
        high_thresh::Real,
        inject_epoch!::Function,
    )
        dat = EegFun.create_test_continuous_data_with_artifacts()
        orig_cols = names(dat.data)

        # Combined & separate non-mutating
        mask = is_fn(dat, thresh)
        @test mask isa AbstractVector{Bool} && length(mask) == nrow(dat.data) && any(mask)
        @test !any(is_fn(dat, high_thresh))
        @test_throws Exception is_fn(dat, -1)

        df_sep = is_fn(dat, thresh; mode = :separate)
        c1_sym, c2_sym = Symbol("$(prefix)_Ch1_$(thresh)"), Symbol("$(prefix)_Ch2_$(thresh)")
        @test df_sep isa DataFrame && c1_sym in propertynames(df_sep) && c2_sym in propertynames(df_sep)
        @test sum(df_sep[!, c1_sym]) == 0 && sum(df_sep[!, c2_sym]) > 0

        df_sub = is_fn(dat, thresh; channel_selection = EegFun.channels([:Ch1]), mode = :separate)
        @test c1_sym in propertynames(df_sub) && c2_sym ∉ propertynames(df_sub)
        @test_throws Exception is_fn(dat, thresh; channel_selection = EegFun.channels(Symbol[]))

        # Mutating combined, separate, custom channel_out, sample_selection & interval_selection
        dat_comb = copy(dat)
        is_mut_fn!(dat_comb, thresh)
        @test setdiff(names(dat_comb.data), orig_cols) == ["$(prefix)_$(thresh)"]

        dat_sep = copy(dat)
        is_mut_fn!(dat_sep, thresh; mode = :separate)
        @test Set(setdiff(names(dat_sep.data), orig_cols)) == Set([string(c1_sym), string(c2_sym)])

        dat_cust = copy(dat)
        is_mut_fn!(dat_cust, thresh; channel_out = :custom_flag)
        @test :custom_flag in propertynames(dat_cust.data) && any(dat_cust.data.custom_flag)
        @test_throws Exception is_mut_fn!(copy(dat), thresh; channel_selection = EegFun.channels(Symbol[]))

        s_mask = falses(nrow(dat.data))
        s_mask[100:200] .= true
        dat_samp = copy(dat)
        is_mut_fn!(dat_samp, thresh; sample_selection = _ -> s_mask)
        @test !any(dat_samp.data[!, Symbol("$(prefix)_$(thresh)")][500:505])

        # Counting & interval_selection
        n_all = n_fn(dat, thresh)
        n_early = n_fn(dat, thresh; interval_selection = EegFun.times(0.0, 0.2))
        @test 0 < n_early < n_all
        @test n_fn(dat, thresh; interval_selection = EegFun.times(0.3, 0.4)) == 0
        @test n_fn(dat, high_thresh) == 0
        @test !any(is_fn(dat, thresh; interval_selection = EegFun.times(0.0, 0.2))[600:end])
        dat_iv = copy(dat)
        is_mut_fn!(dat_iv, thresh; interval_selection = EegFun.times(0.0, 0.2))
        @test !any(dat_iv.data[!, Symbol("$(prefix)_$(thresh)")][600:end])

        cnt_df = n_fn(dat, thresh; mode = :separate)
        @test cnt_df isa DataFrame && cnt_df[cnt_df.channel .== :Ch2, count_col][1] > cnt_df[cnt_df.channel .== :Ch1, count_col][1]
        @test nrow(n_fn(dat, thresh; channel_selection = EegFun.channels([:Ch1]), mode = :separate)) == 1
        @test_throws Exception n_fn(dat, thresh; channel_selection = EegFun.channels(Symbol[]))

        # EpochData & Vector{EpochData}
        epochs = EegFun.create_test_epoch_data(n_epochs = 5, n_channels = 3)
        inject_epoch!(epochs.data[1], :Ch1)
        inject_epoch!(epochs.data[3], :Ch2)
        is_mut_fn!(epochs, 100; mode = :combined)
        ep_col = Symbol("$(prefix)_100")
        @test all(df -> ep_col in propertynames(df), epochs.data)
        @test any(epochs.data[1][!, ep_col]) && any(epochs.data[3][!, ep_col])

        ep_sep = EegFun.create_test_epoch_data(n_epochs = 3, n_channels = 2)
        inject_epoch!(ep_sep.data[1], :Ch1)
        is_mut_fn!(ep_sep, 100; mode = :separate)
        @test Symbol("$(prefix)_Ch1_100") in propertynames(ep_sep.data[1])

        ep_sel = EegFun.create_test_epoch_data(n_epochs = 5, n_channels = 2)
        inject_epoch!(ep_sel.data[1], :Ch1)
        is_mut_fn!(ep_sel, 100; epoch_selection = EegFun.epochs([1, 3]))
        @test ep_col in propertynames(ep_sel.data[1]) && ep_col in propertynames(ep_sel.data[3])
        @test_throws Exception is_mut_fn!(epochs, 100; channel_selection = EegFun.channels(Symbol[]))
        @test_throws Exception is_mut_fn!(epochs, 100; epoch_selection = EegFun.epochs(Int[]))

        ep_list = [EegFun.create_test_epoch_data(n_epochs = 3, n_channels = 2) for _ = 1:2]
        inject_epoch!(ep_list[1].data[1], :Ch1)
        is_mut_fn!(ep_list, 100)
        @test all(ep -> ep_col in propertynames(ep.data[1]), ep_list)
    end

    @testset "detect_eog_onsets!" begin
        dat = EegFun.create_test_continuous_data_with_artifacts()
        EegFun.detect_eog_onsets!(dat, 75, :Ch2, :is_eog_onset)
        @test sum(dat.data.is_eog_onset) == 3 && all(abs.(dat.data.Ch2[dat.data.is_eog_onset]) .>= 75)
        @test dat.data.is_eog_onset[100] && dat.data.is_eog_onset[500] && dat.data.is_eog_onset[800]

        dat.data[!, :is_eog_onset] .= false
        EegFun.detect_eog_onsets!(dat, 20, :Ch2, :is_eog_onset; step_size = 10)
        @test any(dat.data.is_eog_onset)
        @test_throws Exception EegFun.detect_eog_onsets!(dat, 20, :NonExistentChannel, :bad_out)
        @test :bad_out ∉ propertynames(dat.data)

        dat.data[!, :vEOG] = dat.data.Ch2 .+ randn(nrow(dat.data)) .* 5
        dat.data[!, :hEOG] = dat.data.Ch1 .+ randn(nrow(dat.data)) .* 5
        dat.layout = EegFun.Layout(DataFrame(label = [:Ch1, :Ch2, :vEOG, :hEOG], inc = zeros(4), azi = zeros(4)), nothing, nothing, nothing)
        EegFun.detect_eog_onsets!(
            dat,
            Dict(
                "vEOG_criterion" => 50.0,
                "hEOG_criterion" => 30.0,
                "vEOG_channels" => [["Fp1"], ["Fp2"], ["vEOG"]],
                "hEOG_channels" => [["F9"], ["F10"], ["hEOG"]],
            ),
        )
        @test :is_vEOG in propertynames(dat.data) && :is_hEOG in propertynames(dat.data)
    end

    @testset "Extreme & step value detectors" begin
        sig = [1.0, 2.0, 50.0, 3.0, -40.0, 4.0, 5.0]
        @test EegFun._is_extreme_value(sig, 10.0) == [false, false, true, false, true, false, false]
        mask_ext = Vector{Bool}(undef, length(sig))
        EegFun._is_extreme_value!(mask_ext, sig, 10.0)
        @test mask_ext == [false, false, true, false, true, false, false]
        @test EegFun._is_step_value(sig, 10.0) == [false, false, true, true, true, true, false]

        dat_art = EegFun.create_test_continuous_data_with_artifacts()
        ext20 = EegFun.is_extreme_value(dat_art, 20)
        @test all(ext20[100:110]) && all(ext20[500:505]) && all(ext20[800:802])
        @test EegFun.n_extreme_value(dat_art, 20) == 20 && EegFun.n_extreme_value(dat_art, 101) == 17

        _test_detector_contract!(
            EegFun.is_extreme_value,
            EegFun.is_extreme_value!,
            EegFun.n_extreme_value;
            prefix = "is_extreme_value",
            count_col = :n_extreme,
            thresh = 20,
            high_thresh = 1000,
            inject_epoch! = (df, ch) -> (df[!, ch] .+= 200.0),
        )
        _test_detector_contract!(
            EegFun.is_step_value,
            EegFun.is_step_value!,
            EegFun.n_step_value;
            prefix = "is_step_value",
            count_col = :n_step,
            thresh = 50,
            high_thresh = 500,
            inject_epoch! = (df, ch) -> (df[50, ch] = 500.0),
        )
    end

    @testset "Rejection, detect_bad_epochs_automatic & EpochRejectionInfo" begin
        rejs = [EegFun.Rejection(:Ch1, 1), EegFun.Rejection(:Ch2, 2), EegFun.Rejection(:Ch1, 1)]
        @test length(EegFun.unique_rejections(rejs)) == 2
        @test EegFun.unique_channels(rejs) == [:Ch1, :Ch2] && EegFun.unique_epochs(rejs) == [1, 2]

        epochs = EegFun.create_test_epoch_data(n_epochs = 15, n_channels = 3)
        epochs.data[2][!, :Ch1] .= 200.0
        epochs.data[5][!, :Ch2] .= -200.0
        epochs.data[7][!, :Ch1] .= randn(nrow(epochs.data[7])) .* 100.0

        info_abs = EegFun.detect_bad_epochs_automatic(epochs; z_criterion = 0, abs_criterion = 150.0)
        @test isnothing(info_abs.z_rejections) && !isnothing(info_abs.abs_rejections) && !isempty(info_abs.rejected)

        info_z = EegFun.detect_bad_epochs_automatic(epochs; z_criterion = 2.0, abs_criterion = 0)
        @test isnothing(info_z.abs_rejections) && !isnothing(info_z.z_rejections) && length(info_z.z_rejections.z_measures) == 6

        info_both = EegFun.detect_bad_epochs_automatic(epochs; z_criterion = 3.0, abs_criterion = 150.0, z_measures = [:variance, :max])
        @test info_both.z_rejections.z_measures == [:variance, :max] && !isnothing(info_both.abs_rejections)
        @test occursin("EpochRejectionInfo", sprint(show, info_both))
        @test EegFun.get_rejected(info_both) isa Vector{EegFun.Rejection}
        @test EegFun.unique_channels(info_both) isa Vector{Symbol} && EegFun.unique_epochs(info_both) isa Vector{Int}

        df_all = EegFun.all_data(info_both)
        @test names(df_all) == ["condition_number", "condition_name", "epoch", "channel", "reason"] && 2 in df_all.epoch
        @test nrow(EegFun.all_data([info_both, info_both])) == 2 * nrow(df_all)

        ep_list = [epochs, copy(epochs)]
        rej_list = EegFun.detect_bad_epochs_automatic(ep_list; abs_criterion = 150.0, z_criterion = 0)
        @test length(EegFun.get_rejected(rej_list)) == 2 && length(EegFun.unique_rejections(rej_list)) == 2
        @test length(EegFun.unique_channels(rej_list)) == 2 && !isempty(sprint(show, rej_list))

        @test_throws Exception EegFun.detect_bad_epochs_automatic(epochs; z_criterion = 0, abs_criterion = 0)
        @test_throws Exception EegFun.detect_bad_epochs_automatic(epochs; z_criterion = -1, abs_criterion = 100)
        @test_throws Exception EegFun.detect_bad_epochs_automatic(epochs; z_criterion = 3, abs_criterion = -1)
        @test_throws Exception EegFun.detect_bad_epochs_automatic(epochs; z_criterion = 3, z_measures = [:invalid_measure])
        @test_throws Exception EegFun.detect_bad_epochs_automatic(epochs; z_criterion = 3, channel_selection = EegFun.channels(Symbol[]))
    end

    @testset "channel_repairable! and repair_artifacts(!)" begin
        lay3 = EegFun.Layout(
            DataFrame(label = [:Ch1, :Ch2, :Ch3], inc = zeros(3), azi = zeros(3)),
            OrderedDict(
                :Ch1 => EegFun.Neighbors([:Ch2, :Ch3], [1.0, 1.0], [0.5, 0.5]),
                :Ch2 => EegFun.Neighbors([:Ch1, :Ch3], [1.0, 1.0], [0.5, 0.5]),
                :Ch3 => EegFun.Neighbors([:Ch1, :Ch2], [1.0, 1.0], [0.5, 0.5]),
            ),
            nothing,
            nothing,
        )
        epochs = EegFun.create_test_epoch_data(n_epochs = 4, n_channels = 3)
        epochs.layout = lay3
        epochs.data[2][!, :Ch1] .= 200.0
        rej = EegFun.detect_bad_epochs_automatic(epochs; abs_criterion = 150.0, z_criterion = 0)
        @test_throws EegFun.EegFunError EegFun._repair_artifacts_neighbor!(copy(epochs), rej)

        # Auto-calls channel_repairable! when repaired === nothing
        ep_auto = copy(epochs)
        rej_auto = copy(rej)
        EegFun.repair_artifacts!(ep_auto, rej_auto; method = :neighbor_interpolation)
        @test !isnothing(rej_auto.repaired) && ep_auto.data[2][1, :Ch1] != 200.0

        EegFun.channel_repairable!(rej, lay3)
        @test_throws EegFun.EegFunError EegFun.repair_artifacts!(copy(epochs), rej; method = :unknown_method)
        rep_copy = EegFun.repair_artifacts(epochs, rej)
        @test rep_copy !== epochs && rep_copy.data[2][1, :Ch1] != 200.0 && epochs.data[2][1, :Ch1] == 200.0
        @test EegFun._repair_artifacts_neighbor!(epochs, rej) === epochs

        # Vector{EpochData} & lazy neighbor computation + unrepairable skip
        ep_vec = [copy(epochs), copy(epochs)]
        ep_vec[1].data[2][!, :Ch1] .= 200.0
        rej_vec = EegFun.detect_bad_epochs_automatic(ep_vec; abs_criterion = 150.0, z_criterion = 0)
        EegFun.channel_repairable!(rej_vec, lay3)
        @test length(EegFun.repair_artifacts(ep_vec, rej_vec)) == 2

        fresh_lay = EegFun.Layout(
            DataFrame(label = [:Ch1, :Ch2, :Ch3, :Ch4, :Ch5], inc = [0.0, 15.0, 15.0, 90.0, 90.0], azi = [0.0, 0.0, 90.0, 0.0, 180.0]),
            nothing,
            nothing,
            nothing,
        )
        ep_fresh = EegFun.create_test_epoch_data(n_epochs = 2, n_channels = 5)
        ep_fresh.layout = fresh_lay
        ep_fresh.data[1][!, :Ch1] .= 200.0
        rej_fresh = EegFun.detect_bad_epochs_automatic(ep_fresh; abs_criterion = 150.0, z_criterion = 0)
        EegFun.repair_artifacts!(ep_fresh, rej_fresh; method = :neighbor_interpolation)
        @test :Ch1 in rej_fresh.repaired[1] && maximum(abs, ep_fresh.data[1].Ch1) < 50.0

        lay_none = EegFun.Layout(
            DataFrame(label = [:Ch1, :Ch2], inc = zeros(2), azi = zeros(2)),
            OrderedDict(:Ch1 => EegFun.Neighbors(Symbol[], Float64[], Float64[]), :Ch2 => EegFun.Neighbors([:Ch1], [1.0], [1.0])),
            nothing,
            nothing,
        )
        ep_none = EegFun.create_test_epoch_data(n_epochs = 3, n_channels = 2)
        ep_none.layout = lay_none
        ep_none.data[1][!, :Ch1] .= 200.0
        rej_none = EegFun.detect_bad_epochs_automatic(ep_none; abs_criterion = 150.0, z_criterion = 0)
        EegFun.channel_repairable!(rej_none, lay_none)
        @test isempty(rej_none.repaired) && !isempty(rej_none.skipped)
    end

    @testset "Flatline, peak-to-peak, bridged channels & empty-signal guards" begin
        dat = EegFun.create_test_continuous_data_with_artifacts()
        dat.data[1:Int(0.5*dat.sample_rate), :Ch1] .= 0.0
        EegFun.is_flatline!(dat, 0.1, 0.2)
        EegFun.is_flatline!(dat, 0.1, 0.2; mode = :separate, channel_selection = EegFun.channels([:Ch1, :Ch2]))
        @test sum(dat.data.is_flatline) > 0 && "is_flatline_Ch1" in names(dat.data)
        @test length(EegFun.is_flatline(dat, 0.1, 0.2)) == nrow(dat.data) && EegFun.n_flatline(dat, 0.1, 0.2) > 0

        dat.data[500, :Ch2] = 500.0
        dat.data[505, :Ch2] = -500.0
        EegFun.is_peak_to_peak!(dat, 800.0, 0.2)
        @test dat.data[500, :is_peak_to_peak] && EegFun.n_peak_to_peak(dat, 800.0, 0.2) > 0

        # EpochData flatline & peak-to-peak (:combined & :separate)
        ep_fp = EegFun.create_test_epoch_data(n = 200, n_epochs = 3, n_channels = 3)
        ep_fp.data[2][!, :Ch1] .= 5.0
        ep_fp.data[1][50, :Ch1] = 500.0
        ep_fp.data[1][51, :Ch1] = -500.0
        EegFun.is_flatline!(ep_fp, 0.01, 0.05)
        EegFun.is_flatline!(ep_fp, 0.01, 0.05; mode = :separate)
        EegFun.is_peak_to_peak!(ep_fp, 800.0, 0.05)
        EegFun.is_peak_to_peak!(ep_fp, 800.0, 0.05; mode = :separate)
        @test count(ep_fp.data[2].is_flatline) == 200 && count(ep_fp.data[2].is_flatline_Ch1) == 200
        @test count(ep_fp.data[1].is_peak_to_peak) > 0 && count(ep_fp.data[1].is_peak_to_peak_Ch1) > 0

        # Bridged channels on ContinuousData, EpochData, <2 channels & Vector{EpochData}
        dat_br = EegFun.create_test_continuous_data_with_artifacts()
        dat_br.data[!, :Ch2] = dat_br.data.Ch1 .+ 0.001 .* randn(nrow(dat_br.data))
        @test nrow(EegFun.find_bridged_channels(dat_br, 0.99)) >= 1

        ep_br = EegFun.create_test_epoch_data(n = 200, n_epochs = 4, n_channels = 3)
        for ep in ep_br
            ep[!, :Ch2] = ep.Ch1 .+ 0.001 .* randn(nrow(ep))
        end
        @test nrow(EegFun.find_bridged_channels(ep_br, 0.99)) >= 1
        res_empty_br = EegFun.find_bridged_channels(ep_br, 0.99; channel_selection = EegFun.channels([:Ch1]))
        @test nrow(res_empty_br) == 0 && names(res_empty_br) == ["channel_1", "channel_2", "correlation"]
        @test length(EegFun.find_bridged_channels([ep_br, ep_br], 0.99)) == 2

        @test isempty(EegFun._is_step_value(Float64[], 50.0))
        @test isempty(EegFun._is_peak_to_peak(Float64[], 50.0, 10))
        @test isempty(EegFun._is_flatline(Float64[], 0.1, 10))
        @test isempty(EegFun._is_extreme_value(Float64[], 50.0))
        @test EegFun._is_extreme_value(Float32[10.0, 60.0], 50.0) == [false, true]
        m32 = falses(2)
        EegFun._is_extreme_value!(m32, Float32[10.0, 60.0], 50.0)
        @test m32 == [false, true]
    end

    @testset "Spherical spline repair, subset_bad_data & MultiDataFrameEeg overloads" begin
        dat_sp = EegFun.create_test_epoch_data(n = 80, n_epochs = 3, n_channels = 4)
        dat_sp.layout.data.x3 = [0.0, 0.1, 0.0, -0.1]
        dat_sp.layout.data.y3 = [0.1, 0.0, -0.1, 0.0]
        dat_sp.layout.data.z3 = [0.0, 0.0, 0.0, 0.0]
        orig_ep1_ch2 = copy(dat_sp.data[1].Ch2)
        dat_sp.data[2][!, :Ch2] .= 300.0
        for ch in (:Ch1, :Ch2, :Ch3, :Ch4)
            dat_sp.data[3][!, ch] .= 300.0
        end
        rej_sp = EegFun.detect_bad_epochs_automatic(dat_sp; abs_criterion = 150.0, z_criterion = 0)
        dat_nonmut = EegFun.repair_artifacts(dat_sp, rej_sp; method = :spherical_spline)
        @test isnothing(rej_sp.repaired) && all(dat_nonmut.data[2].Ch2 .!= 300.0) && dat_nonmut.data[1].Ch2 ≈ orig_ep1_ch2

        @test_logs (:warn,) match_mode = :any EegFun.repair_artifacts!(dat_sp, rej_sp; method = :spherical_spline)
        @test rej_sp.repaired[2] == [:Ch2] && length(rej_sp.skipped[3]) == 4 && !haskey(rej_sp.repaired, 3)

        mktempdir() do tmpdir
            ep_a = EegFun.create_test_epoch_data_vector(conditions = 1:1, n_epochs = 2)
            for fn in ("1_epochs.jld2", "sub-1_epochs.jld2", "10_epochs.jld2", "21_epochs.jld2")
                JLD2.jldsave(joinpath(tmpdir, fn); data = ep_a)
            end
            JLD2.jldsave(
                joinpath(tmpdir, "epoch_summary.jld2");
                data = DataFrame(file = ["1", "10", "21"], condition = [1, 1, 1], percentage = [50.0, 95.0, 95.0]),
            )
            JLD2.jldsave(joinpath(tmpdir, "file_summary.jld2"); data = DataFrame(file = ["1", "10", "21"]))
            EegFun.subset_bad_data(tmpdir, 75.0; subset_directory = "bad_out")
            out_sub = joinpath(tmpdir, "bad_out")
            @test isfile(joinpath(out_sub, "1_epochs.jld2")) && isfile(joinpath(out_sub, "sub-1_epochs.jld2"))
            @test !isfile(joinpath(out_sub, "10_epochs.jld2")) && !isfile(joinpath(out_sub, "21_epochs.jld2"))
        end

        ep_det = EegFun.create_test_epoch_data(n = 100, n_epochs = 4, n_channels = 2)
        ep_det.data[1][10, :Ch1] = 500.0
        ep_det.data[1][11, :Ch1] = -500.0
        ep_det.data[2][!, :Ch2] .= 0.0
        @test any(EegFun.is_extreme_value(ep_det, 200.0)[1]) && EegFun.n_extreme_value(ep_det, 200.0) >= 2
        @test EegFun.n_extreme_value(ep_det, 200.0; mode = :separate) isa DataFrame
        @test any(EegFun.is_step_value(ep_det, 200.0)[1]) && EegFun.n_step_value(ep_det, 200.0) >= 1
        @test any(EegFun.is_peak_to_peak(ep_det, 800.0, 0.05)[1]) && EegFun.n_peak_to_peak(ep_det, 800.0, 0.05) >= 1
        @test all(EegFun.is_flatline(ep_det, 0.01, 0.05)[2]) && EegFun.n_flatline(ep_det, 0.01, 0.05) == 100
        @test length(EegFun.is_extreme_value([ep_det, ep_det], 200.0)) == 2
        @test length(EegFun.n_extreme_value([ep_det, ep_det], 200.0)) == 2

        ep_schema = copy(ep_det)
        EegFun.is_extreme_value!(ep_schema, 200.0; epoch_selection = EegFun.epochs([1]))
        EegFun.is_step_value!(ep_schema, 200.0; epoch_selection = EegFun.epochs([1]))
        EegFun.is_peak_to_peak!(ep_schema, 800.0, 0.05; epoch_selection = EegFun.epochs([1]))
        EegFun.is_flatline!(ep_schema, 0.01, 0.05; epoch_selection = EegFun.epochs([2]))
        @test nrow(EegFun.all_data(ep_schema)) == 400
        @test !any(ep_schema.data[3][!, Symbol("is_extreme_value_200.0")]) && !any(ep_schema.data[1].is_flatline)

        ep_flat_z = EegFun.create_test_epoch_data(n = 100, n_epochs = 8, n_channels = 2, seed = 77)
        ep_flat_z.data[1][!, :Ch1] .= 0.0
        ep_flat_z.data[2][!, :Ch1] .= 0.01 .* randn(MersenneTwister(88), 100)
        ep_flat_z.data[2][50, :Ch1] = 250.0
        @test any(
            r -> r.epoch == 2 && r.channel == :Ch1,
            EegFun.detect_bad_epochs_automatic(ep_flat_z; z_criterion = 2.0, abs_criterion = 0.0).rejected,
        )
    end
end
