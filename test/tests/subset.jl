using Test
using EegFun
using DataFrames

@testset "Subsetting functions" begin

    @testset "subset with interval_selection" begin
        # Create test ERP data
        erp = EegFun.create_test_erp_data()

        # Test interval_selection with tuple
        @testset "Tuple interval" begin
            result = EegFun.subset(erp, interval_selection = (0.1, 0.3))
            @test nrow(result.data) < nrow(erp.data)
            @test all(0.1 .<= result.data.time .<= 0.3)
        end

        # Test interval_selection with nothing (all samples)
        @testset "Nothing interval (all samples)" begin
            result = EegFun.subset(erp, interval_selection = nothing)
            @test nrow(result.data) == nrow(erp.data)
        end

        # Test interval_selection with Tuple
        @testset "Tuple" begin
            interval = (0.2, 0.4)
            result = EegFun.subset(erp, interval_selection = interval)
            @test nrow(result.data) < nrow(erp.data)
            @test all(0.2 .<= result.data.time .<= 0.4)
        end
    end

    @testset "subset with sample_selection" begin
        erp = EegFun.create_test_erp_data()

        # Test sample_selection with predicate function
        @testset "Predicate function" begin
            result = EegFun.subset(erp, sample_selection = x -> x.time .>= 0.2)
            @test nrow(result.data) < nrow(erp.data)
            @test all(result.data.time .>= 0.2)
        end

        # Test sample_selection with samples() helper
        @testset "samples() helper with column" begin
            erp = EegFun.create_test_erp_data()
            # Use column-based filtering (samples() is for boolean predicates, not time)
            result = EegFun.subset(erp, sample_selection = x -> x.time .>= 0.1 .&& x.time .<= 0.3)
            @test nrow(result.data) < nrow(erp.data)
            @test all(0.1 .<= result.data.time .<= 0.3)
        end
    end

    @testset "Combining sample_selection and interval_selection" begin
        erp = EegFun.create_test_erp_data()

        # Both should work together with AND logic
        @testset "Both parameters together" begin
            # interval_selection: 0.1-0.5, sample_selection: >= 0.3
            # Result should be intersection: 0.3-0.5
            result = EegFun.subset(erp, interval_selection = (0.1, 0.5), sample_selection = x -> x.time .>= 0.3)
            @test nrow(result.data) < nrow(erp.data)
            @test all(result.data.time .>= 0.3)
            @test all(result.data.time .<= 0.5)
        end

        # Verify non-overlapping selections return empty
        @testset "Non-overlapping selections" begin
            result = EegFun.subset(erp, interval_selection = (0.1, 0.2), sample_selection = x -> x.time .>= 0.8)
            @test nrow(result.data) == 0
        end
    end

    @testset "subset with different data types" begin
        @testset "ErpData" begin
            erp = EegFun.create_test_erp_data()
            result = EegFun.subset(erp, interval_selection = (0.2, 0.4))
            @test result isa EegFun.ErpData
            @test nrow(result.data) < nrow(erp.data)
            @test all(0.2 .<= result.data.time .<= 0.4)
        end

        @testset "EpochData" begin
            epochs = EegFun.create_test_epoch_data()

            # Test interval_selection
            @testset "Interval selection" begin
                result = EegFun.subset(epochs, interval_selection = (0.1, 0.3))
                @test result isa EegFun.EpochData
                @test all(df -> nrow(df) < nrow(epochs.data[1]), result.data)
                @test all(df -> all(0.1 .<= df.time .<= 0.3), result.data)
            end

            # Test channel selection with epochs
            @testset "Channel and interval selection" begin
                result = EegFun.subset(epochs, channel_selection = EegFun.channels([:Ch1]), interval_selection = (0.0, 0.5))
                @test result isa EegFun.EpochData
                # Check only Ch1 remains in all epochs
                for df in result.data
                    channel_cols = names(df)[occursin.(r"^Ch[0-9]+$", names(df))]
                    @test length(channel_cols) == 1
                    @test "Ch1" in channel_cols
                end
            end

            # Test combining both selections
            @testset "Combined selections" begin
                result = EegFun.subset(epochs, interval_selection = (0.1, 0.8), sample_selection = x -> x.time .< 0.5)
                @test result isa EegFun.EpochData
                @test all(df -> all(df.time .>= 0.1) && all(df.time .< 0.5), result.data)
            end

            @testset "Epoch selection" begin
                # Index-based selection
                result_idx = EegFun.subset(epochs, epoch_selection = EegFun.epochs([1]))
                @test length(result_idx.data) == 1

                # Predicate-based selection over DataFrames
                # Suppose we want to only keep epochs that actually have rows
                result_pred = EegFun.subset(epochs, epoch_selection = EegFun.epochs(df -> nrow(df) > 0))
                @test length(result_pred.data) == length(epochs.data)

                # Predicate that removes all
                result_pred_empty = EegFun.subset(epochs, epoch_selection = EegFun.epochs(df -> nrow(df) > 100000))
                @test length(result_pred_empty.data) == 0
            end
        end

        @testset "ContinuousData" begin
            continuous = EegFun.create_test_continuous_data()

            # Test interval_selection
            @testset "Interval selection" begin
                result = EegFun.subset(continuous, interval_selection = (1.0, 2.0))
                @test result isa EegFun.ContinuousData
                @test nrow(result.data) < nrow(continuous.data)
                @test all(1.0 .<= result.data.time .<= 2.0)
            end

            # Test channel selection
            @testset "Channel selection" begin
                result = EegFun.subset(continuous, channel_selection = EegFun.channels([:Ch1, :Ch2]))
                @test result isa EegFun.ContinuousData
                channel_cols = names(result.data)[occursin.(r"^Ch[0-9]+$", names(result.data))]
                @test length(channel_cols) == 2
                @test "Ch1" in channel_cols
                @test "Ch2" in channel_cols
            end

            # Test combining both selections
            @testset "Combined selections" begin
                result = EegFun.subset(continuous, channel_selection = EegFun.channels([:Ch1]), interval_selection = (2.0, 4.0))
                @test result isa EegFun.ContinuousData
                @test all(2.0 .<= result.data.time .<= 4.0)
                channel_cols = names(result.data)[occursin.(r"^Ch[0-9]+$", names(result.data))]
                @test length(channel_cols) == 1
            end

            # Test sample_selection with continuous data
            @testset "Sample selection predicate" begin
                result = EegFun.subset(continuous, sample_selection = x -> x.time .>= 3.0)
                @test result isa EegFun.ContinuousData
                @test all(result.data.time .>= 3.0)
            end
        end

        @testset "TimeFreqData" begin
            tf = EegFun.create_test_tf_data(n_channels = 3)

            @testset "Channel selection" begin
                result = EegFun.subset(tf, channel_selection = EegFun.channels([:Ch1]))
                @test result isa EegFun.TimeFreqData
                # Both power and phase should be subset
                power_chs = names(result.data_power)[occursin.(r"^Ch[0-9]+$", names(result.data_power))]
                phase_chs = names(result.data_phase)[occursin.(r"^Ch[0-9]+$", names(result.data_phase))]
                @test length(power_chs) == 1
                @test "Ch1" in power_chs
                @test power_chs == phase_chs
            end

            @testset "Power and phase consistency" begin
                result = EegFun.subset(tf, channel_selection = EegFun.channels([:Ch1, :Ch2]))
                # Same number of rows in both
                @test nrow(result.data_power) == nrow(result.data_phase)
                # Same columns in both
                @test names(result.data_power) == names(result.data_phase)
            end

            @testset "Defaults (all data)" begin
                result = EegFun.subset(tf)
                @test nrow(result.data_power) == nrow(tf.data_power)
                @test nrow(result.data_phase) == nrow(tf.data_phase)
            end
        end

        @testset "Vector{TimeFreqData}" begin
            tfs = [
                EegFun.create_test_tf_data(condition = 1),
                EegFun.create_test_tf_data(condition = 2),
                EegFun.create_test_tf_data(condition = 3),
            ]

            @testset "Condition selection" begin
                result = EegFun.subset(tfs, condition_selection = EegFun.conditions([1, 2]))
                @test length(result) == 2
                @test all(r -> r isa EegFun.TimeFreqData, result)
            end

            @testset "Channel selection on vector" begin
                result = EegFun.subset(tfs, channel_selection = EegFun.channels([:Ch1]))
                @test length(result) == 3
                for r in result
                    power_chs = names(r.data_power)[occursin.(r"^Ch[0-9]+$", names(r.data_power))]
                    @test length(power_chs) == 1
                end
            end
        end
    end

    @testset "subset with channel_selection" begin
        erp = EegFun.create_test_erp_data()

        @testset "Channel and interval selection" begin
            result = EegFun.subset(erp, channel_selection = EegFun.channels([:Ch1, :Ch2]), interval_selection = (0.1, 0.3))

            # Check channels
            channel_cols = names(result.data)[occursin.(r"^Ch[0-9]+$", names(result.data))]
            @test length(channel_cols) == 2
            @test "Ch1" in channel_cols
            @test "Ch2" in channel_cols

            # Check time
            @test all(0.1 .<= result.data.time .<= 0.3)
        end
    end

    @testset "Edge cases" begin
        erp = EegFun.create_test_erp_data()

        @testset "Empty interval" begin
            # Interval outside data range
            result = EegFun.subset(erp, interval_selection = (10.0, 20.0))
            @test nrow(result.data) == 0
        end

        @testset "Single sample interval" begin
            time_points = erp.data.time
            single_time = time_points[1]
            result = EegFun.subset(erp, interval_selection = (single_time, single_time))
            @test nrow(result.data) >= 1  # At least the matching point
        end

        @testset "Defaults (all data)" begin
            result = EegFun.subset(erp)
            @test nrow(result.data) == nrow(erp.data)
        end
    end

    @testset "times() helper function" begin
        @testset "times() returns nothing" begin
            @test isnothing(EegFun.times())
        end

        @testset "times(start, stop) returns tuple" begin
            @test EegFun.times(0.1, 0.5) == (0.1, 0.5)
        end

        @testset "times(tuple) returns tuple" begin
            @test EegFun.times((0.2, 0.6)) == (0.2, 0.6)
        end



        @testset "times(Tuple) passes through" begin
            interval = (0.3, 0.7)
            @test EegFun.times(interval) === interval
        end
    end
end
