# Tests for region-of-interest and multi-channel averaging.

using Test
using DataFrames
using EegFun
using JLD2
using Statistics

@testset "channel_average" begin

    dat = EegFun.create_test_continuous_data(n = 500)

    # Append averaged columns only (Symbols input)
    EegFun.channel_average!(dat, channel_selections = [EegFun.channels([:Ch1, :Ch2])])

    @test :Ch1_Ch2 ∈ propertynames(dat.data)
    @test all(dat.data.Ch1_Ch2 .== (dat.data.Ch1 .+ dat.data.Ch2) ./ 2)

    # Reduce to only averages and create averaged layout
    dat = EegFun.create_test_continuous_data(n = 500)
    dat = EegFun.channel_average(dat, channel_selections = [EegFun.channels([:Ch1, :Ch2])]; reduce = true)
    @test all(propertynames(dat.data) .== [:time, :sample, :trigger, :Ch1_Ch2])
    @test size(dat.layout.data, 1) == 1
    @test :inc ∈ propertynames(dat.layout.data)
    @test :azi ∈ propertynames(dat.layout.data)

    # Append averaged channels to layout
    dat = EegFun.create_test_continuous_data(n = 500)
    dat = EegFun.channel_average(dat, channel_selections = [EegFun.channels([:Ch2, :Ch3])])
    @test :Ch2_Ch3 ∈ propertynames(dat.data)
    @test any(dat.layout.data.label .== :Ch2_Ch3)

    # Mutating version
    dat = EegFun.create_test_continuous_data(n = 500)
    EegFun.channel_average!(dat, channel_selections = [EegFun.channels([:Ch1, :Ch3])])
    @test :Ch1_Ch3 ∈ propertynames(dat.data)
    @test :Ch1_Ch2 ∉ propertynames(dat.data)

    # Auto-label :avg for all channels
    dat = EegFun.create_test_continuous_data(n = 500)
    dat = EegFun.channel_average(dat, channel_selections = [EegFun.channels([:Ch1, :Ch2, :Ch3])]; reduce = true)
    @test :avg ∈ propertynames(dat.data)

    # Custom output_labels applied 
    dat = EegFun.create_test_continuous_data(n = 500)
    EegFun.channel_average!(dat, channel_selections = [EegFun.channels([:Ch1, :Ch2])]; output_labels = [:output_label])
    @test :output_label ∈ propertynames(dat.data)

    # Duplicate labels in layout (should accumulate)
    # First add B_C, then add again to verify rows accumulate
    dat = EegFun.create_test_continuous_data(n = 500)
    EegFun.channel_average!(dat, channel_selections = [EegFun.channels([:Ch2, :Ch3])])
    n1 = sum(dat.layout.data.label .== :Ch2_Ch3)
    EegFun.channel_average!(dat, channel_selections = [EegFun.channels([:Ch2, :Ch3])])
    n2 = sum(dat.layout.data.label .== :Ch2_Ch3)
    @test n1 == 1 && n2 == 2

    # EpochData append and reduce
    # Create simple EpochData (2 epochs) from dat
    dat = EegFun.create_test_epoch_data(n = 500)
    dat = EegFun.channel_average(dat, channel_selections = [EegFun.channels([:Ch1, :Ch2])])

    @test :Ch1_Ch2 ∈ propertynames(dat.data[1]) && :Ch1_Ch2 ∈ propertynames(dat.data[2])
    @test :Ch1_Ch2 ∈ propertynames(dat.data[3]) && :Ch1_Ch2 ∈ propertynames(dat.data[4])

    dat = EegFun.create_test_epoch_data(n = 500)
    dat = EegFun.channel_average(dat, channel_selections = [EegFun.channels([:Ch1, :Ch2])]; reduce = true)

    # Expect meta columns (leading) + A_B only
    cols = propertynames(dat.data[1])
    @test cols[end] == :Ch1_Ch2
    @test :Ch1 ∉ cols && :Ch2 ∉ cols && :Ch3 ∉ cols

    # ErpData reduce path
    dat = EegFun.create_test_epoch_data(n = 500)
    dat = EegFun.channel_average(dat, channel_selections = [EegFun.channels([:Ch1, :Ch2])]; reduce = true)
    # condition and condition_name are in struct, not DataFrame
    @test all(propertynames(dat.data[1]) .== [:time, :sample, :epoch, :Ch1_Ch2])
    @test all(propertynames(dat.data[end]) .== [:time, :sample, :epoch, :Ch1_Ch2])
    @test hasproperty(dat, :condition)
    @test hasproperty(dat, :condition_name)

    # Test default behavior (average all channels)
    dat = EegFun.create_test_epoch_data(n = 500)
    EegFun.channel_average!(dat)  # Should use default channels() from second function
    @test :avg ∈ propertynames(dat.data[1])
    @test all(dat.data[1].avg .== (dat.data[1].Ch1 .+ dat.data[1].Ch2 .+ dat.data[1].Ch3) ./ 3)
    @test all(dat.data[end].avg .== (dat.data[end].Ch1 .+ dat.data[end].Ch2 .+ dat.data[end].Ch3) ./ 3)

    # Test single channel selection with custom label
    dat = EegFun.create_test_epoch_data(n = 500)
    EegFun.channel_average!(dat, channel_selections = [EegFun.channels([:Ch1, :Ch2])], output_labels = [:custom])
    @test :custom ∈ propertynames(dat.data[1])
    @test :custom ∈ propertynames(dat.data[end])
    @test all(dat.data[1].custom .== (dat.data[1].Ch1 .+ dat.data[1].Ch2) ./ 2)
    @test all(dat.data[end].custom .== (dat.data[end].Ch1 .+ dat.data[end].Ch2) ./ 2)

end

@testset "Batch Average Channels" begin
    isdefined(@__MODULE__, :run_batch_jld2_contract_tests!) ||
        Base.include(@__MODULE__, joinpath(pkgdir(EegFun), "test", "tests", "test_helpers.jl"))

    run_batch_jld2_contract_tests!(
        (pat; kwargs...) -> EegFun.channel_average(pat, [EegFun.channels([:Ch1, :Ch2]), EegFun.channels([:Ch3, :Ch4])]; kwargs...);
        make_data = _ -> EegFun.create_test_batch_erp_data(n_conditions = 2, n_channels = 7),
        pattern = "erps_unrejected",
        log_file = "channel_average.log",
        log_snippets = ["Batch channel averaging started", "channel_average"],
    )

    mktempdir() do test_dir
        orig_erps = EegFun.create_test_batch_erp_data(n_conditions = 2, n_channels = 7)
        jldsave(joinpath(test_dir, "1_erps_unrejected.jld2"); data = orig_erps)
        jldsave(joinpath(test_dir, "2_erps_uncorrected.jld2"); data = orig_erps)

        # Default labels, layout append, metadata & sample_rate preservation
        out_def = joinpath(test_dir, "out_def")
        EegFun.channel_average(
            "erps_unrejected",
            [EegFun.channels([:Ch1, :Ch2]), EegFun.channels([:Ch3, :Ch4])];
            input_dir = test_dir,
            output_dir = out_def,
        )
        erps_def = load(joinpath(out_def, "1_erps_unrejected.jld2"), "data")
        @test erps_def[1].data.avg_1 ≈ (erps_def[1].data.Ch1 .+ erps_def[1].data.Ch2) ./ 2
        @test hasproperty(erps_def[1].data, :avg_2)
        @test nrow(erps_def[1].layout.data) == nrow(orig_erps[1].layout.data) + 2
        @test :avg_1 in erps_def[1].layout.data.label
        @test erps_def[1].sample_rate == orig_erps[1].sample_rate
        @test hasproperty(erps_def[1].data, :time) && hasproperty(erps_def[1].data, :sample)

        # Reduce mode, custom labels, overlapping groups, single-channel & all-channels groups
        out_red = joinpath(test_dir, "out_red")
        EegFun.channel_average(
            "erps_unrejected",
            [EegFun.channels([:Ch1, :Ch2]), EegFun.channels([:Ch3, :Ch4])];
            output_labels = [:Group1, :Group2],
            input_dir = test_dir,
            output_dir = out_red,
            reduce = true,
        )
        erps_red = load(joinpath(out_red, "1_erps_unrejected.jld2"), "data")
        @test erps_red[1].layout.data.label == [:Group1, :Group2]
        @test hasproperty(erps_red[1].data, :Group1) && !hasproperty(erps_red[1].data, :Ch1)
        @test hasproperty(erps_red[1].data, :time)

        out_ov = joinpath(test_dir, "out_ov")
        EegFun.channel_average(
            "erps_unrejected",
            [EegFun.channels([:Ch1, :Ch2, :Ch5]), EegFun.channels([:Ch5]), EegFun.channels_not([:Ch5])];
            output_labels = [:Ov1, :Single5, :NotCh5],
            input_dir = test_dir,
            output_dir = out_ov,
        )
        erps_ov = load(joinpath(out_ov, "1_erps_unrejected.jld2"), "data")
        @test erps_ov[1].data.Ov1 ≈ (erps_ov[1].data.Ch1 .+ erps_ov[1].data.Ch2 .+ erps_ov[1].data.Ch5) ./ 3
        @test erps_ov[1].data.Single5 ≈ erps_ov[1].data.Ch5
        @test hasproperty(erps_ov[1].data, :NotCh5)

        out_all = joinpath(test_dir, "out_all")
        EegFun.channel_average(
            "erps_unrejected",
            [EegFun.channels()];
            output_labels = [:global_avg],
            input_dir = test_dir,
            output_dir = out_all,
            reduce = true,
        )
        erps_all = load(joinpath(out_all, "1_erps_unrejected.jld2"), "data")
        expected_all = vec(mean(Matrix(orig_erps[1].data[!, [:Ch1, :Ch2, :Ch3, :Ch4, :Ch5, :Ch6, :Ch7]]), dims = 2))
        @test erps_all[1].data.global_avg ≈ expected_all

        # Pattern matching across multiple files
        @test EegFun.channel_average(
            "erps",
            [EegFun.channels([:Ch1, :Ch2])];
            input_dir = test_dir,
            output_dir = joinpath(test_dir, "pat_all"),
        ).success == 2

        # EpochData batch support
        ep_dir = mkpath(joinpath(test_dir, "ep_in"))
        jldsave(joinpath(ep_dir, "1_epochs.jld2"); data = [EegFun.create_test_epoch_data()])
        ep_out = joinpath(test_dir, "ep_out")
        @test EegFun.channel_average(
            "epochs",
            [EegFun.channels([:Ch1, :Ch2])];
            output_labels = [:Group1],
            input_dir = ep_dir,
            output_dir = ep_out,
        ).success == 1
        ep_loaded = load(joinpath(ep_out, "1_epochs.jld2"), "data")
        @test ep_loaded[1] isa EegFun.EpochData && hasproperty(ep_loaded[1].data[1], :Group1)

        # Error cases: mismatched labels, nonexistent channels, empty selection
        @test_throws Exception EegFun.channel_average(
            "erps_unrejected",
            [EegFun.channels([:Ch1, :Ch2]), EegFun.channels([:Ch3, :Ch4])];
            output_labels = [:OnlyOne],
            input_dir = test_dir,
        )
        res_inv = EegFun.channel_average(
            "erps_unrejected",
            [EegFun.channels([:NonExistent])];
            input_dir = test_dir,
            output_dir = joinpath(test_dir, "inv_out"),
        )
        @test res_inv.success == 0 && res_inv.errors == 1
        res_emp =
            EegFun.channel_average("erps_unrejected", [ch -> Symbol[]]; input_dir = test_dir, output_dir = joinpath(test_dir, "emp_out"))
        @test res_emp.success == 0 && res_emp.errors == 1
    end
end

@testset "TF Channel Average" begin

    @testset "Mutating: basic average all channels" begin
        tf = EegFun.create_test_tf_data(n_channels = 3, power_offset = 3.0, noise = 0.0)

        EegFun.channel_average!(tf)

        # Should have added :avg column to both power and phase
        @test hasproperty(tf.data_power, :avg)
        @test hasproperty(tf.data_phase, :avg)

        # :avg should be the mean of Ch1, Ch2, Ch3 (all 3.0 with noise=0)
        @test all(tf.data_power[!, :avg] .≈ 3.0)
    end

    @testset "Mutating: specific channel group" begin
        tf = EegFun.create_test_tf_data(n_channels = 3, power_offset = 2.0, noise = 0.0)

        EegFun.channel_average!(tf, channel_selections = [EegFun.channels([:Ch1, :Ch2])], output_labels = [:frontal])

        @test hasproperty(tf.data_power, :frontal)
        @test hasproperty(tf.data_phase, :frontal)
        @test all(tf.data_power[!, :frontal] .≈ 2.0)
    end

    @testset "Mutating: reduce mode" begin
        tf = EegFun.create_test_tf_data(n_channels = 3, power_offset = 1.0, noise = 0.0)

        EegFun.channel_average!(tf, channel_selections = [EegFun.channels([:Ch1, :Ch2])], output_labels = [:avg_group], reduce = true)

        # Should only have metadata + averaged columns
        @test hasproperty(tf.data_power, :time)
        @test hasproperty(tf.data_power, :freq)
        @test hasproperty(tf.data_power, :avg_group)
        @test !hasproperty(tf.data_power, :Ch1)
        @test !hasproperty(tf.data_power, :Ch3)

        # Same for phase
        @test hasproperty(tf.data_phase, :avg_group)
        @test !hasproperty(tf.data_phase, :Ch1)
    end

    @testset "Non-mutating: returns copy" begin
        tf = EegFun.create_test_tf_data(n_channels = 3)

        result = EegFun.channel_average(tf)

        # Original should not be modified
        @test !hasproperty(tf.data_power, :avg)

        # Result should have the averaged column
        @test hasproperty(result.data_power, :avg)
        @test hasproperty(result.data_phase, :avg)
    end

    @testset "Non-mutating: Vector{TimeFreqData}" begin
        tfs = [EegFun.create_test_tf_data(condition = 1), EegFun.create_test_tf_data(condition = 2)]

        results = EegFun.channel_average(tfs)

        @test length(results) == 2
        @test all(hasproperty(r.data_power, :avg) for r in results)
        @test all(hasproperty(r.data_phase, :avg) for r in results)

        # Originals untouched
        @test !hasproperty(tfs[1].data_power, :avg)
    end

    @testset "Both DataFrames stay in sync" begin
        tf = EegFun.create_test_tf_data(n_channels = 3)

        EegFun.channel_average!(tf, channel_selections = [EegFun.channels([:Ch1, :Ch2])], output_labels = [:test_avg])

        # Same number of rows
        @test nrow(tf.data_power) == nrow(tf.data_phase)

        # Same columns added
        @test propertynames(tf.data_power) == propertynames(tf.data_phase)
    end

    @testset "TimeFreqEpochData: power and phase both averaged" begin
        tfe = EegFun.create_test_tf_epoch_data(n_channels = 3, n_epochs = 2, power_offset = 2.0, noise = 0.0)

        EegFun.channel_average!(tfe, channel_selections = [EegFun.channels([:Ch1, :Ch2])], output_labels = [:grp])

        for i = 1:2
            @test hasproperty(tfe.data_power[i], :grp)
            @test hasproperty(tfe.data_phase[i], :grp)
            @test all(tfe.data_power[i][!, :grp] .≈ 2.0)
        end
    end

    @testset "TimeFreqEpochData: reduce mode" begin
        tfe = EegFun.create_test_tf_epoch_data(n_channels = 3, n_epochs = 2)

        EegFun.channel_average!(tfe, channel_selections = [EegFun.channels([:Ch1, :Ch2])], output_labels = [:grp], reduce = true)

        for dfs in (tfe.data_power, tfe.data_phase), df in dfs
            @test propertynames(df) == [:time, :freq, :grp]
        end
        @test EegFun.channel_labels(tfe) == [:grp]
    end

    @testset "Ergonomic channel_selection and Vector channel_average!" begin
        # Singular keyword channel_selection
        dat = EegFun.create_test_continuous_data(n = 100)
        EegFun.channel_average!(dat; channel_selection = [:Ch1, :Ch2], output_labels = "my_mean")
        @test :my_mean ∈ propertynames(dat.data)
        @test dat.data.my_mean ≈ (dat.data.Ch1 .+ dat.data.Ch2) ./ 2

        # Un-nested single group in channel_selections
        dat2 = EegFun.create_test_continuous_data(n = 100)
        EegFun.channel_average!(dat2; channel_selections = ["Ch1", "Ch2"], output_labels = :my_mean2)
        @test :my_mean2 ∈ propertynames(dat2.data)

        # In-place Vector method
        dats = [EegFun.create_test_continuous_data(n = 50), EegFun.create_test_continuous_data(n = 50)]
        EegFun.channel_average!(dats; channel_selection = [:Ch1, :Ch2], output_labels = :vec_avg)
        @test :vec_avg ∈ propertynames(dats[1].data)
        @test :vec_avg ∈ propertynames(dats[2].data)
    end

    @testset "Default extra-channel exclusion and circular TF phase mean" begin
        # 1. Default selection excludes extra channels while keeping layout channels containing "_"
        df_ext = DataFrame(
            time = [0.0, 0.01, 0.02],
            sample = [1, 2, 3],
            Ch1 = [2.0, 2.0, 2.0],
            Fp1_bip = [4.0, 4.0, 4.0],
            vEOG = [1000.0, 1000.0, 1000.0],
            interval = [true, false, true],
        )
        lay_ext = EegFun.Layout(DataFrame(label = [:Ch1, :Fp1_bip], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing, nothing)
        cont_ext = EegFun.ContinuousData("ext_test", df_ext, lay_ext, 100, EegFun.AnalysisInfo())
        avg_ext = EegFun.channel_average(cont_ext)
        @test avg_ext.data.avg ≈ [3.0, 3.0, 3.0] # mean of Ch1 (2) and Fp1_bip (4), excluding vEOG
        @test EegFun.extra_labels(avg_ext) == [:vEOG, :interval]
        @test EegFun.channel_labels(avg_ext) == [:Ch1, :Fp1_bip, :avg]

        # Symmetric equatorial ROI (T7/T8 centroid = (0, 0, 0)) stays in layout for reduce = false and reduce = true
        df_eq = DataFrame(time = [0.0, 0.01], sample = [1, 2], T7 = [2.0, 4.0], T8 = [6.0, 8.0])
        lay_eq = EegFun.Layout(DataFrame(label = [:T7, :T8], inc = [90.0, 90.0], azi = [180.0, 0.0]), nothing, nothing, nothing)
        erp_eq = EegFun.ErpData("eq", 1, "c1", df_eq, lay_eq, 100, EegFun.AnalysisInfo(), 10)
        erp_eq_red = EegFun.channel_average(erp_eq; channel_selections = [EegFun.channels([:T7, :T8])], reduce = true)
        @test EegFun.channel_labels(erp_eq_red) == [:avg]
        @test erp_eq_red.data.avg ≈ [4.0, 6.0]

        # 2. Circular mean on TimeFreqData.data_phase and TimeFreqEpochData.data_phase across ±π branch cut
        tf_circ = EegFun.create_test_tf_data(n_channels = 2, power_offset = 4.0, noise = 0.0)
        tf_circ.data_phase[!, :Ch1] .= 3.1
        tf_circ.data_phase[!, :Ch2] .= -3.1
        EegFun.channel_average!(tf_circ)
        # Arithmetic mean would give 0.0 (opposite direction!), whereas circular mean gives ±π
        @test all(abs.(abs.(tf_circ.data_phase[!, :avg]) .- π) .< 1e-6)

        tfe_circ = EegFun.create_test_tf_epoch_data(n_channels = 2, n_epochs = 2, power_offset = 4.0, noise = 0.0)
        for df_p in tfe_circ.data_phase
            df_p[!, :Ch1] .= 3.1
            df_p[!, :Ch2] .= -3.1
        end
        EegFun.channel_average!(tfe_circ)
        for df_p in tfe_circ.data_phase
            @test all(abs.(abs.(df_p[!, :avg]) .- π) .< 1e-6)
        end
    end
end
