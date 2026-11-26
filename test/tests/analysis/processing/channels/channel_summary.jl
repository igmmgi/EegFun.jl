# Tests for per-channel amplitude and variance summary tables.

using Test
using EegFun
using DataFrames
using JLD2
using Statistics
using CSV


@testset "EegFun" begin
    @testset "channel_summary" begin

        @testset "_channel_summary_impl core function" begin
            dat = EegFun.create_test_continuous_data()

            # Test basic functionality
            result = EegFun._channel_summary_impl(dat.data, collect(1:50), [:Ch1, :Ch2])

            @test result isa DataFrame
            @test nrow(result) == 2  # Two channels
            @test :channel in propertynames(result)
            @test :min in propertynames(result)
            @test :max in propertynames(result)
            @test :std in propertynames(result)
            @test :range in propertynames(result)
            @test :var in propertynames(result)
            @test :zvar in propertynames(result)

            # Test column order
            @test propertynames(result)[1] == :channel

            # Test range calculation
            @test all(result.range .≈ result.max .- result.min)

            # Test variance calculation
            selected_data = dat.data[collect(1:50), [:Ch1, :Ch2]]
            expected_var = var.(eachcol(selected_data))
            @test result.var ≈ expected_var

            # Test that channels are correctly identified
            @test Set(result.channel) == Set([:Ch1, :Ch2])
        end

        @testset "_channel_summary_impl input validation" begin
            dat = EegFun.create_test_continuous_data(n_channels = 4)

            @test_throws Exception EegFun._channel_summary_impl(dat.data, Int[], [:Ch1, :Ch2])
            @test_throws Exception EegFun._channel_summary_impl(dat.data, [1, 2, 3], Symbol[])
            @test_throws Exception EegFun._channel_summary_impl(dat.data, [1, 2, 3], [:NonExistent])
            @test_throws Exception EegFun._channel_summary_impl(dat.data, [3000], [:Ch1])  # Out of bounds
            @test_throws Exception EegFun._channel_summary_impl(dat.data, [0], [:Ch1])     # Below bounds

        end

        @testset "_channel_summary_impl edge cases" begin
            dat = EegFun.create_test_continuous_data(n_channels = 4)

            # Test single sample
            result = EegFun._channel_summary_impl(dat.data, [1], [:Ch1])
            @test nrow(result) == 1
            @test result.range[1] == 0.0  # min == max for single sample
            @test isnan(result.std[1]) || result.std[1] == 0.0  # std is NaN for single sample in Julia

            # Test single channel
            result = EegFun._channel_summary_impl(dat.data, collect(1:10), [:Ch1])
            @test nrow(result) == 1
            @test result.channel[1] == :Ch1

            # Test with constant signal (zero variance case)
            dat.data[!, :Ch4] .= 1.0
            result = EegFun._channel_summary_impl(dat.data, collect(1:50), [:Ch4])  # Ch4 is constant
            @test result.var[1] ≈ 0.0
            @test result.zvar[1] == 0.0  # Should handle zero variance case
        end

        @testset "channel_summary SingleDataFrameEeg" begin
            dat = EegFun.create_test_continuous_data(n_channels = 4)

            # Test basic functionality with defaults
            result = EegFun.channel_summary(dat)
            @test result isa DataFrame
            @test nrow(result) == 4  # Four layout channels (Ch1, Ch2, Ch3, Ch4)
            @test :channel in propertynames(result)

            # Test that all layout channels are included by default
            @test Set(result.channel) == Set([:Ch1, :Ch2, :Ch3, :Ch4])

            # Test specific channel selection
            result_subset = EegFun.channel_summary(dat, channel_selection = EegFun.channels([:Ch1, :Ch2]))
            @test nrow(result_subset) == 2
            @test Set(result_subset.channel) == Set([:Ch1, :Ch2])

            # Test sample selection (first half)
            n_samples = nrow(dat.data)
            result_samples = EegFun.channel_summary(dat, sample_selection = x -> 1:nrow(x) .<= div(nrow(x), 2))
            @test result_samples isa DataFrame
            @test nrow(result_samples) == 4

            # Test include_meta and include_extra flags
            result_meta = EegFun.channel_summary(dat, include_meta = true)
            @test :time in result_meta.channel || :trigger in result_meta.channel  # Should include meta columns

            # Test statistics are reasonable
            result = EegFun.channel_summary(dat)
            @test all(result.range .>= 0)  # Range should be non-negative
            @test all(result.var .>= 0)    # Variance should be non-negative
            @test all(result.std .>= 0)    # Standard deviation should be non-negative
        end

        @testset "channel_summary SingleDataFrameEeg input validation" begin
            # Test empty data
            empty_df = DataFrame()
            layout = EegFun.Layout(DataFrame(label = Symbol[], inc = Float64[], azi = Float64[]), nothing, nothing, nothing)
            empty_dat = EegFun.ContinuousData("test_data", empty_df, layout, 100, EegFun.AnalysisInfo())
            @test_throws Exception EegFun.channel_summary(empty_dat)
        end

        @testset "channel_summary MultiDataFrameEeg" begin
            epoch_dat = EegFun.create_test_epoch_data()

            # Test basic functionality
            result = EegFun.channel_summary(epoch_dat)
            @test result isa DataFrame
            @test :epoch in propertynames(result)
            @test :channel in propertynames(result)

            # Test that we have results for each epoch and channel
            n_epochs = EegFun.n_epochs(epoch_dat)
            n_channels = 3  # A, B, C
            @test nrow(result) == n_epochs * n_channels

            # Test epoch numbers are preserved
            original_epoch_numbers = [df.epoch[1] for df in epoch_dat.data]
            result_epoch_numbers = unique(result.epoch)
            @test Set(result_epoch_numbers) == Set(original_epoch_numbers)

            # Test channel selection works across epochs
            result_subset = EegFun.channel_summary(epoch_dat, channel_selection = EegFun.channels([:Ch1]))
            @test nrow(result_subset) == n_epochs  # One channel per epoch
            @test all(result_subset.channel .== :Ch1)

            # Test that statistics vary across epochs (should be different)
            epochs_Ch1 = result[result.channel .== :Ch1, :]
            if nrow(epochs_Ch1) > 1
                # Should have some variation across epochs (unless data is identical)
                @test length(unique(epochs_Ch1.var)) > 1 || all(epochs_Ch1.var .== 0)
            end
        end

        @testset "channel_summary MultiDataFrameEeg input validation" begin
            # Test empty epoch data
            empty_layout = EegFun.Layout(DataFrame(label = Symbol[], inc = Float64[], azi = Float64[]), nothing, nothing, nothing)
            empty_epochs = EegFun.EpochData("test_data", 1, "condition_1", DataFrame[], empty_layout, 100, EegFun.AnalysisInfo())
            @test_throws Exception EegFun.channel_summary(empty_epochs)
        end

        @testset "channel_summary statistical properties" begin
            dat = EegFun.create_test_continuous_data(n_channels = 4)
            dat.data[!, :Ch4] .= 1.0
            result = EegFun.channel_summary(dat)

            # Test that constant channel has zero variance
            constant_row = result[result.channel .== :Ch4, :]
            @test nrow(constant_row) == 1
            @test constant_row.var[1] ≈ 0.0 atol = 1e-10
            @test constant_row.std[1] ≈ 0.0 atol = 1e-10
            @test constant_row.range[1] ≈ 0.0 atol = 1e-10

            # Test that varying channels have non-zero variance
            varying_rows = result[result.channel .!= :Ch4, :]
            @test all(varying_rows.var .> 0)
            @test all(varying_rows.std .> 0)
            @test all(varying_rows.range .> 0)

            # Test z-scored variance properties
            if length(unique(result.var)) > 1  # If there's variance in variances
                @test abs(mean(result.zvar)) < 1e-10  # Should be approximately zero mean
                @test abs(std(result.zvar) - 1.0) < 1e-10  # Should have unit std
            end
        end

        @testset "channel_summary consistency" begin
            dat = EegFun.create_test_continuous_data()

            # Test that results are consistent between calls
            result1 = EegFun.channel_summary(dat)
            result2 = EegFun.channel_summary(dat)
            @test result1 == result2

            # Test that subsetting gives expected results
            full_result = EegFun.channel_summary(dat)
            subset_result = EegFun.channel_summary(dat, channel_selection = EegFun.channels([:Ch1, :Ch2]))

            # The subset should match the corresponding rows from the full result
            full_subset = full_result[in([:Ch1, :Ch2]).(full_result.channel), :]
            @test subset_result.channel == full_subset.channel
            @test subset_result.min ≈ full_subset.min
            @test subset_result.max ≈ full_subset.max
            @test subset_result.std ≈ full_subset.std
        end

        @testset "channel_summary integration with selection functions" begin
            dat = EegFun.create_test_continuous_data(n_channels = 4)

            # Test with different selection functions
            result_all = EegFun.channel_summary(dat, channel_selection = EegFun.channels())
            @test nrow(result_all) == 4

            # Test channel selection by name
            result_specific = EegFun.channel_summary(dat, channel_selection = EegFun.channels([:Ch1]))
            @test nrow(result_specific) == 1
            @test result_specific.channel[1] == :Ch1

            # Test sample selection by range
            result_half = EegFun.channel_summary(dat, sample_selection = x -> 1:nrow(x) .<= div(nrow(x), 2))
            @test result_half isa DataFrame
            @test nrow(result_half) == 4  # Same number of channels

            # MultiDataFrameEeg preserves :epoch column IDs (e.g. [2, 5] after rejection)
            ep_nonseq = EegFun.create_test_epoch_data(n = 30, n_epochs = 2, n_channels = 2)
            ep_nonseq.data[1].epoch .= 2
            ep_nonseq.data[2].epoch .= 5
            res_nonseq = EegFun.channel_summary(ep_nonseq)
            @test sort(unique(res_nonseq.epoch)) == [2, 5]

            # Default samples() hits fast path and matches explicit sample predicate identically
            res_pred = EegFun.channel_summary(ep_nonseq; sample_selection = df -> trues(nrow(df)))
            @test res_nonseq == res_pred
        end

    end # channel_summary testset
end # EegFun testset


@testset "Batch Channel Summary" begin
    mktempdir() do test_dir
        for p in (1, 2)
            jldsave(joinpath(test_dir, "$(p)_erps_unrejected.jld2"); data = EegFun.create_test_batch_erp_data(n_conditions = 2))
        end

        @testset "Basic CSV schema, logging, existing dir & overwrite" begin
            out_dir = mkpath(joinpath(test_dir, "summary_output"))
            touch(joinpath(out_dir, "dummy.txt"))
            csv_file = joinpath(out_dir, "channel_summary.csv")
            write(csv_file, "stale")

            EegFun.channel_summary("erps_unrejected"; input_dir = test_dir, output_dir = out_dir)
            @test isfile(joinpath(out_dir, "dummy.txt"))
            @test isfile(csv_file)
            log_txt = read(joinpath(out_dir, "channel_summary.log"), String)
            @test occursin("Batch channel summary started", log_txt) && occursin("erps_unrejected", log_txt)

            results = CSV.read(csv_file, DataFrame)
            @test nrow(results) == 12
            for col in (:file, :condition, :channel, :min, :max, :std, :range, :var, :zvar)
                @test hasproperty(results, col)
            end
            @test Set(results.channel) == Set(["Ch1", "Ch2", "Ch3"])
            @test Set(results.file) == Set(["1_erps_unrejected", "2_erps_unrejected"])
            for g in groupby(results, [:file, :condition])
                @test nrow(g) == 3
                @test mean(g.zvar) ≈ 0.0 atol = 1e-10
                @test std(g.zvar; corrected = false) ≈ 1.0 atol = 0.2
            end

            # Custom output_file name
            cust_dir = joinpath(test_dir, "custom_out")
            EegFun.channel_summary("erps_unrejected"; input_dir = test_dir, output_dir = cust_dir, output_file = "my_custom_summary")
            @test isfile(joinpath(cust_dir, "my_custom_summary.csv")) && isfile(joinpath(cust_dir, "my_custom_summary.log"))
        end

        @testset "Selection filters, error handling & partial failures" begin
            for (tag, kw, expected_rows) in (
                ("p1", (; participant_selection = EegFun.participants(1)), 6),
                ("p12", (; participant_selection = EegFun.participants([1, 2])), 12),
                ("c1", (; condition_selection = EegFun.conditions(1)), 6),
                ("c12", (; condition_selection = EegFun.conditions([1, 2])), 12),
                ("ch12", (; channel_selection = EegFun.channels([:Ch1, :Ch2])), 8),
                ("not3", (; channel_selection = EegFun.channels_not([:Ch3])), 8),
                (
                    "comb",
                    (;
                        participant_selection = EegFun.participants(1),
                        condition_selection = EegFun.conditions(1),
                        channel_selection = EegFun.channels([:Ch1, :Ch2]),
                    ),
                    2,
                ),
            )
                odir = joinpath(test_dir, "sel_$tag")
                EegFun.channel_summary("erps_unrejected"; input_dir = test_dir, output_dir = odir, kw...)
                @test nrow(CSV.read(joinpath(odir, "channel_summary.csv"), DataFrame)) == expected_rows
            end

            @test_throws Exception EegFun.channel_summary("erps_unrejected"; input_dir = joinpath(test_dir, "nonexistent"))
            empty_dir = mkpath(joinpath(test_dir, "empty_match"))
            @test isnothing(EegFun.channel_summary("erps_unrejected"; input_dir = empty_dir))

            oor_dir = joinpath(test_dir, "oor_out")
            EegFun.channel_summary(
                "erps_unrejected";
                input_dir = test_dir,
                output_dir = oor_dir,
                condition_selection = EegFun.conditions(5),
            )
            @test !isfile(joinpath(oor_dir, "channel_summary.csv"))

            pdir = mkpath(joinpath(test_dir, "partial_in"))
            jldsave(joinpath(pdir, "1_erps_unrejected.jld2"); data = EegFun.create_test_batch_erp_data(n_conditions = 2))
            jldsave(joinpath(pdir, "2_erps_unrejected.jld2"); data = "invalid_data")
            pout = joinpath(test_dir, "partial_out")
            EegFun.channel_summary("erps_unrejected"; input_dir = pdir, output_dir = pout)
            @test nrow(CSV.read(joinpath(pout, "channel_summary.csv"), DataFrame)) == 6
        end

        @testset "Analytical stats, non-sequential conditions & _fill_zvar!" begin
            sdir = mkpath(joinpath(test_dir, "stats_in"))
            n_s = 100
            df = DataFrame(
                time = collect(range(0.0, 1.0, length = n_s)),
                sample = 1:n_s,
                condition = fill(1, n_s),
                Ch1 = fill(5.0, n_s),
                Ch2 = collect(1.0:100.0),
                Ch3 = zeros(n_s),
            )
            lay = EegFun.Layout(DataFrame(label = [:Ch1, :Ch2, :Ch3], inc = zeros(3), azi = zeros(3)), nothing, nothing, nothing)
            jldsave(
                joinpath(sdir, "1_erps_stats.jld2");
                data = [EegFun.ErpData("test_data", 1, "condition_1", df, lay, 256.0, EegFun.AnalysisInfo(), 1)],
            )
            sout = joinpath(test_dir, "stats_out")
            EegFun.channel_summary("erps_stats"; input_dir = sdir, output_dir = sout)
            res = CSV.read(joinpath(sout, "channel_summary.csv"), DataFrame)

            ch1 = only(eachrow(res[res.channel .== "Ch1", :]))
            @test ch1.min ≈ 5.0 && ch1.max ≈ 5.0 && isapprox(ch1.std, 0.0; atol = 1e-10) && isapprox(ch1.range, 0.0; atol = 1e-10)
            ch2 = only(eachrow(res[res.channel .== "Ch2", :]))
            @test ch2.min ≈ 1.0 && ch2.max ≈ 100.0 && ch2.range ≈ 99.0 && ch2.std ≈ std(1.0:100.0) && ch2.var ≈ var(1.0:100.0)
            ch3 = only(eachrow(res[res.channel .== "Ch3", :]))
            @test ch3.min ≈ 0.0 && ch3.max ≈ 0.0 && isapprox(ch3.std, 0.0; atol = 1e-10)

            cdir = mkpath(joinpath(test_dir, "nonseq_cond"))
            jldsave(
                joinpath(cdir, "1_erps_nonseq.jld2");
                data = [
                    EegFun.create_test_erp_data(participant = 1, condition = 2, n_channels = 2),
                    EegFun.create_test_erp_data(participant = 1, condition = 5, n_channels = 2),
                ],
            )
            nout = joinpath(test_dir, "summary_nonseq")
            EegFun.channel_summary("erps_nonseq"; input_dir = cdir, output_dir = nout)
            @test sort(unique(CSV.read(joinpath(nout, "channel_summary.csv"), DataFrame).condition)) == [2, 5]

            zv = [1.0, 1.0, 1.0]
            EegFun._fill_zvar!(zv, [4.5, 4.5, 4.5])
            @test zv == [0.0, 0.0, 0.0]
        end
    end
end
