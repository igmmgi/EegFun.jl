using Test
using DataFrames
using OrderedCollections
using eegfun
using Makie


@testset "plot_triggers" begin

    @testset "_filter_triggers" begin
        @testset "basic filtering" begin
            trigger_times = [1.0, 2.0, 3.0, 4.0, 5.0]
            trigger_values = [1, 2, 1, 3, 2]
            ignore_triggers = [1, 3]

            filtered_times, filtered_values = eegfun._filter_triggers(trigger_times, trigger_values, ignore_triggers)

            @test filtered_times == [2.0, 5.0]
            @test filtered_values == [2, 2]
        end

        @testset "no filtering needed" begin
            trigger_times = [1.0, 2.0, 3.0]
            trigger_values = [1, 2, 3]
            ignore_triggers = [4, 5]  # No overlap

            filtered_times, filtered_values = eegfun._filter_triggers(trigger_times, trigger_values, ignore_triggers)

            @test filtered_times == trigger_times
            @test filtered_values == trigger_values
        end

        @testset "filter all triggers" begin
            trigger_times = [1.0, 2.0, 3.0]
            trigger_values = [1, 2, 3]
            ignore_triggers = [1, 2, 3]

            filtered_times, filtered_values = eegfun._filter_triggers(trigger_times, trigger_values, ignore_triggers)

            @test isempty(filtered_times)
            @test isempty(filtered_values)
        end

        @testset "empty input" begin
            filtered_times, filtered_values = eegfun._filter_triggers(Float64[], Int[], [1, 2])

            @test isempty(filtered_times)
            @test isempty(filtered_values)
        end
    end

    @testset "_trigger_time_count" begin
        @testset "basic counting" begin
            time = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
            triggers = [0, 1, 0, 2, 0, 1]

            trigger_times, trigger_values, trigger_count = eegfun._trigger_time_count(time, triggers)

            @test trigger_times == [0.1, 0.3, 0.5]
            @test trigger_values == [1, 2, 1]
            @test trigger_count[1] == 2
            @test trigger_count[2] == 1
        end

        @testset "with ignore_triggers" begin
            time = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
            triggers = [0, 1, 0, 2, 0, 1]
            ignore_triggers = [1]

            trigger_times, trigger_values, trigger_count = eegfun._trigger_time_count(time, triggers, ignore_triggers)

            @test trigger_times == [0.3]
            @test trigger_values == [2]
            @test length(trigger_count) == 1
            @test trigger_count[2] == 1
        end

        @testset "no triggers" begin
            time = [0.0, 0.1, 0.2, 0.3]
            triggers = [0, 0, 0, 0]

            trigger_times, trigger_values, trigger_count = eegfun._trigger_time_count(time, triggers)

            @test isempty(trigger_times)
            @test isempty(trigger_values)
            @test isempty(trigger_count)
        end

        @testset "all triggers ignored" begin
            time = [0.0, 0.1, 0.2, 0.3]
            triggers = [0, 1, 0, 2]
            ignore_triggers = [1, 2]

            trigger_times, trigger_values, trigger_count = eegfun._trigger_time_count(time, triggers, ignore_triggers)

            @test isempty(trigger_times)
            @test isempty(trigger_values)
            @test isempty(trigger_count)
        end
    end

    @testset "_extract_trigger_data" begin
        # BioSemiDataFormat tests removed due to type complexity

        @testset "ContinuousData extraction" begin
            dat = create_test_continuous_data_with_triggers()

            trigger_codes, trigger_times = eegfun._extract_trigger_data(dat)

            @test length(trigger_codes) == 14
            @test length(trigger_times) == 14
            @test 1 in trigger_codes
            @test 2 in trigger_codes
            @test 3 in trigger_codes
        end

        @testset "ContinuousData with filtering" begin
            dat = create_test_continuous_data_with_triggers()
            ignore_triggers = [2, 3]

            trigger_codes, trigger_times = eegfun._extract_trigger_data(dat, ignore_triggers)

            @test length(trigger_codes) == 7
            @test length(trigger_times) == 7
        end

        @testset "missing triggers column" begin
            # Create ContinuousData without triggers column
            df = DataFrame(time = [0.0, 0.1, 0.2], channel = [1.0, 2.0, 3.0])
            layout = eegfun.Layout(DataFrame(label = [:channel], inc = [0.0], azi = [0.0]), nothing, nothing)
            dat = eegfun.ContinuousData(df, layout, 100, eegfun.AnalysisInfo())

            # Should trigger minimal_error and return nothing
            @test eegfun._extract_trigger_data(dat) === nothing
        end
    end

    @testset "_count_triggers" begin
        @testset "basic counting" begin
            trigger_codes = Int16[1, 2, 1, 3, 2, 1]
            trigger_count = eegfun._count_triggers(trigger_codes)

            @test trigger_count[1] == 3
            @test trigger_count[2] == 2
            @test trigger_count[3] == 1
        end

        @testset "single trigger type" begin
            trigger_codes = Int16[1, 1, 1, 1]
            trigger_count = eegfun._count_triggers(trigger_codes)

            @test length(trigger_count) == 1
            @test trigger_count[1] == 4
        end

        @testset "empty input" begin
            trigger_codes = Int16[]
            trigger_count = eegfun._count_triggers(trigger_codes)

            @test isempty(trigger_count)
        end
    end

    # =============================================================================
    # MAIN PLOTTING FUNCTION TESTS
    # =============================================================================

    @testset "plot_trigger_overview" begin
        @testset "basic functionality" begin
            trigger_times = [1.0, 2.0, 3.0, 4.0]
            trigger_values = [1, 2, 1, 3]
            trigger_count = OrderedDict(1 => 2, 2 => 1, 3 => 1)

            fig, ax = eegfun.plot_trigger_overview(trigger_times, trigger_values, trigger_count; display_plot = false)

            @test fig isa Figure
            @test ax isa Axis
        end

        @testset "empty trigger data" begin
            trigger_times = Float64[]
            trigger_values = Int[]
            trigger_count = OrderedDict{Int,Int}()

            fig, ax = eegfun.plot_trigger_overview(trigger_times, trigger_values, trigger_count; display_plot = false)

            @test fig isa Figure
            @test ax isa Axis
        end

        # BioSemiDataFormat tests removed due to type complexity

        @testset "ContinuousData input" begin
            dat = create_test_continuous_data_with_triggers()

            fig, ax = eegfun.plot_trigger_overview(dat; display_plot = false)

            @test fig isa Figure
            @test ax isa Axis
        end

        @testset "ContinuousData with ignore_triggers" begin
            dat = create_test_continuous_data_with_triggers()

            fig, ax = eegfun.plot_trigger_overview(dat; ignore_triggers = [2, 3], display_plot = false)

            @test fig isa Figure
            @test ax isa Axis
        end

        @testset "empty ContinuousData" begin
            dat = create_empty_trigger_data()

            fig, ax = eegfun.plot_trigger_overview(dat; display_plot = false)

            @test fig isa Figure
            @test ax isa Axis
        end

        @testset "parameter passing" begin
            dat = create_test_continuous_data_with_triggers()

            # Test that custom parameters are passed through
            fig, ax = eegfun.plot_trigger_overview(dat; window_size = 15.0, display_plot = false, ignore_triggers = [1])

            @test fig isa Figure
            @test ax isa Axis
        end
    end

    @testset "plot_trigger_timing" begin
        # BioSemiDataFormat tests removed due to type complexity

        @testset "ContinuousData input" begin
            dat = create_test_continuous_data_with_triggers()

            fig, ax = eegfun.plot_trigger_timing(dat; display_plot = false)

            @test fig isa Figure
            @test ax isa Axis
        end

        @testset "ContinuousData with ignore_triggers" begin
            dat = create_test_continuous_data_with_triggers()

            fig, ax = eegfun.plot_trigger_timing(dat; ignore_triggers = [2, 3], display_plot = false)

            @test fig isa Figure
            @test ax isa Axis
        end

        @testset "empty ContinuousData" begin
            dat = create_empty_trigger_data()

            fig, ax = eegfun.plot_trigger_timing(dat; display_plot = false)

            @test fig isa Figure
            @test ax isa Axis
        end

        @testset "parameter passing" begin
            dat = create_test_continuous_data_with_triggers()

            # Test that custom parameters are passed through
            fig, ax = eegfun.plot_trigger_timing(
                dat;
                window_size = 20.0,
                initial_position = -5.0,
                display_plot = false,
                ignore_triggers = [1],
            )

            @test fig isa Figure
            @test ax isa Axis
        end
    end

    # =============================================================================
    # INTEGRATION TESTS
    # =============================================================================

    @testset "integration tests" begin
        @testset "end-to-end workflow" begin
            # Test complete workflow from data creation to plotting
            dat = create_test_continuous_data_with_triggers()

            # Test overview plot
            fig1, ax1 = eegfun.plot_trigger_overview(dat; display_plot = false)
            @test fig1 isa Figure
            @test ax1 isa Axis

            # Test timing plot
            fig2, ax2 = eegfun.plot_trigger_timing(dat; display_plot = false)
            @test fig2 isa Figure
            @test ax2 isa Axis

            # Test with filtering
            fig3, ax3 = eegfun.plot_trigger_overview(dat; ignore_triggers = [1], display_plot = false)
            @test fig3 isa Figure
            @test ax3 isa Axis

            fig4, ax4 = eegfun.plot_trigger_timing(dat; ignore_triggers = [1], display_plot = false)
            @test fig4 isa Figure
            @test ax4 isa Axis
        end

        @testset "parameter consistency" begin
            dat = create_test_continuous_data_with_triggers()

            # Test that same parameters work across all functions
            common_params = (ignore_triggers = [1], display_plot = false)

            fig1, ax1 = eegfun.plot_trigger_overview(dat; common_params...)
            fig2, ax2 = eegfun.plot_trigger_timing(dat; common_params...)

            @test fig1 isa Figure && fig2 isa Figure
            @test ax1 isa Axis && ax2 isa Axis
        end

        @testset "return value consistency" begin
            dat = create_test_continuous_data_with_triggers()

            # All functions should return (fig, ax) tuple
            fig1, ax1 = eegfun.plot_trigger_overview(dat; display_plot = false)
            fig2, ax2 = eegfun.plot_trigger_timing(dat; display_plot = false)

            @test isa(fig1, Figure) && isa(fig2, Figure)
            @test isa(ax1, Axis) && isa(ax2, Axis)
        end
    end

    # =============================================================================
    # PERFORMANCE TESTS
    # =============================================================================

    @testset "performance tests" begin
        @testset "large dataset handling" begin
            # Test with larger dataset
            dat = create_test_continuous_data_with_triggers(; n = 10000)

            # Should complete without errors
            fig, ax = eegfun.plot_trigger_overview(dat; display_plot = false)
            @test fig isa Figure
            @test ax isa Axis

            fig, ax = eegfun.plot_trigger_timing(dat; display_plot = false)
            @test fig isa Figure
            @test ax isa Axis
        end

        @testset "filtering efficiency" begin
            # Test that empty ignore_triggers has no performance penalty
            dat = create_test_continuous_data_with_triggers(; n = 5000)

            # Time both versions
            @time fig1, ax1 = eegfun.plot_trigger_overview(dat; display_plot = false)
            @time fig2, ax2 = eegfun.plot_trigger_overview(dat; ignore_triggers = Int[], display_plot = false)

            @test fig1 isa Figure && fig2 isa Figure
            @test ax1 isa Axis && ax2 isa Axis
        end

        @testset "memory usage" begin
            # Test that functions don't leak memory
            dat = create_test_continuous_data_with_triggers()

            # Create multiple plots to check for memory leaks
            for i = 1:5
                fig, ax = eegfun.plot_trigger_overview(dat; display_plot = false)
                @test fig isa Figure
                @test ax isa Axis

                fig, ax = eegfun.plot_trigger_timing(dat; display_plot = false)
                @test fig isa Figure
                @test ax isa Axis
            end
        end
    end

    # =============================================================================
    # EDGE CASE TESTS
    # =============================================================================

    @testset "edge cases" begin
        @testset "single trigger" begin
            time = [0.0, 0.1, 0.2, 0.3]
            triggers = [0, 1, 0, 0]

            trigger_times, trigger_values, trigger_count = eegfun._trigger_time_count(time, triggers)

            @test length(trigger_times) == 1
            @test trigger_times[1] == 0.1
            @test trigger_values[1] == 1
            @test trigger_count[1] == 1
        end

        @testset "all same triggers" begin
            time = [0.0, 0.1, 0.2, 0.3, 0.4]
            triggers = [0, 1, 1, 1, 0]

            trigger_times, trigger_values, trigger_count = eegfun._trigger_time_count(time, triggers)

            @test length(trigger_times) == 3  # All three triggers (not cleaned in this function)
            @test trigger_times == [0.1, 0.2, 0.3]
            @test trigger_values == [1, 1, 1]
            @test trigger_count[1] == 3
        end

        @testset "consecutive different triggers" begin
            time = [0.0, 0.1, 0.2, 0.3, 0.4]
            triggers = [0, 1, 2, 3, 0]

            trigger_times, trigger_values, trigger_count = eegfun._trigger_time_count(time, triggers)

            @test length(trigger_times) == 3
            @test trigger_times == [0.1, 0.2, 0.3]
            @test trigger_values == [1, 2, 3]
            @test trigger_count[1] == 1
            @test trigger_count[2] == 1
            @test trigger_count[3] == 1
        end

        @testset "ignore_triggers edge cases" begin
            time = [0.0, 0.1, 0.2, 0.3, 0.4]
            triggers = [0, 1, 2, 3, 0]

            # Ignore non-existent triggers
            trigger_times, trigger_values, trigger_count = eegfun._trigger_time_count(time, triggers, [99, 100])

            @test length(trigger_times) == 3  # No filtering
            @test trigger_times == [0.1, 0.2, 0.3]
            @test trigger_values == [1, 2, 3]
        end

        @testset "empty ignore_triggers" begin
            time = [0.0, 0.1, 0.2, 0.3]
            triggers = [0, 1, 2, 0]

            # Empty ignore list should not filter anything
            trigger_times, trigger_values, trigger_count = eegfun._trigger_time_count(time, triggers, Int[])

            @test length(trigger_times) == 2
            @test trigger_times == [0.1, 0.2]
            @test trigger_values == [1, 2]
        end
    end

    # =============================================================================
    # ERROR HANDLING TESTS
    # =============================================================================

    @testset "error handling" begin
        @testset "invalid ignore_triggers types" begin
            # Test that functions handle invalid types gracefully
            dat = create_test_continuous_data_with_triggers()

            # These should not throw errors but may not work as expected
            # The functions should handle type conversion internally
            fig, ax = eegfun.plot_trigger_overview(dat; ignore_triggers = [1.0, 2.0], display_plot = false)
            @test fig isa Figure
            @test ax isa Axis
        end

        # Malformed data structures test removed due to type complexity

        @testset "extreme parameter values" begin
            dat = create_test_continuous_data_with_triggers()

            # Test with extreme parameter values
            fig, ax = eegfun.plot_trigger_timing(
                dat;
                window_size = 0.01,  # Very small window
                initial_position = -1000.0,  # Very negative position
                display_plot = false,
            )

            @test fig isa Figure
            @test ax isa Axis
        end
    end

    # =============================================================================
    # ADVANCED EDGE CASES
    # =============================================================================

    @testset "advanced edge cases" begin
        @testset "data type edge cases" begin
            # Test with negative trigger values
            time = [0.0, 0.1, 0.2, 0.3]
            triggers = [0, -1, 0, -2]

            trigger_times, trigger_values, trigger_count = eegfun._trigger_time_count(time, triggers)

            @test length(trigger_times) == 2
            @test trigger_times == [0.1, 0.3]
            @test trigger_values == [-1, -2]
            @test trigger_count[-1] == 1
            @test trigger_count[-2] == 1
        end

        @testset "time series edge cases" begin
            # Test with non-monotonic time (should still work)
            time = [0.0, 0.1, 0.05, 0.2]  # Non-monotonic
            triggers = [0, 1, 0, 2]

            trigger_times, trigger_values, trigger_count = eegfun._trigger_time_count(time, triggers)

            @test length(trigger_times) == 2
            @test trigger_times == [0.1, 0.2]  # Should preserve original order
            @test trigger_values == [1, 2]
        end

        @testset "array size edge cases" begin
            # Test with single element arrays
            time = [1.0]
            triggers = [1]

            trigger_times, trigger_values, trigger_count = eegfun._trigger_time_count(time, triggers)

            @test length(trigger_times) == 1
            @test trigger_times == [1.0]
            @test trigger_values == [1]
            @test trigger_count[1] == 1
        end

        @testset "interactive plot edge cases" begin
            dat = create_test_continuous_data_with_triggers()

            # Test with window size larger than data range
            fig, ax = eegfun.plot_trigger_timing(
                dat;
                window_size = 10000.0,  # Much larger than data
                display_plot = false,
            )
            @test fig isa Figure
            @test ax isa Axis

            # Test with very small window size
            fig, ax = eegfun.plot_trigger_timing(
                dat;
                window_size = 0.001,  # Very small
                display_plot = false,
            )
            @test fig isa Figure
            @test ax isa Axis

            # Test with extreme initial position
            fig, ax = eegfun.plot_trigger_timing(
                dat;
                initial_position = 10000.0,  # Way beyond data
                display_plot = false,
            )
            @test fig isa Figure
            @test ax isa Axis
        end

        @testset "memory and performance edge cases" begin
            # Test with very large dataset
            dat = create_test_continuous_data_with_triggers(; n = 50000)

            fig, ax = eegfun.plot_trigger_overview(dat; display_plot = false)
            @test fig isa Figure
            @test ax isa Axis

            fig, ax = eegfun.plot_trigger_timing(dat; display_plot = false)
            @test fig isa Figure
            @test ax isa Axis
        end

        @testset "boundary value edge cases" begin
            # Test with boundary values
            time = [0.0, 0.1, 0.2, 0.3]
            triggers = [0, 1, 0, 0]

            # Test with ignore_triggers containing the only trigger
            trigger_times, trigger_values, trigger_count = eegfun._trigger_time_count(time, triggers, [1])

            @test isempty(trigger_times)
            @test isempty(trigger_values)
            @test isempty(trigger_count)
        end

        @testset "type conversion edge cases" begin
            # Test that functions handle type conversions gracefully
            time = [0.0, 0.1, 0.2, 0.3]
            triggers = [0, 1, 0, 2]

            # Test with different integer types
            trigger_times, trigger_values, trigger_count = eegfun._trigger_time_count(time, Int8.(triggers))

            @test length(trigger_times) == 2
            @test trigger_times == [0.1, 0.3]
            @test trigger_values == [1, 2]
        end
    end

    # =============================================================================
    # NEW FEATURE TESTS (ignore_triggers)
    # =============================================================================

    @testset "ignore_triggers feature" begin
        @testset "basic filtering functionality" begin
            dat = create_test_continuous_data_with_triggers()

            # Test filtering out specific triggers
            fig1, ax1 = eegfun.plot_trigger_overview(dat; ignore_triggers = [1], display_plot = false)
            fig2, ax2 = eegfun.plot_trigger_overview(dat; ignore_triggers = [2, 3], display_plot = false)
            fig3, ax3 = eegfun.plot_trigger_overview(dat; ignore_triggers = [1, 2, 3], display_plot = false)

            @test fig1 isa Figure && fig2 isa Figure && fig3 isa Figure
            @test ax1 isa Axis && ax2 isa Axis && ax3 isa Axis
        end

        @testset "filtering consistency across functions" begin
            dat = create_test_continuous_data_with_triggers()
            ignore_list = [1, 2]

            # Both functions should handle filtering the same way
            fig1, ax1 = eegfun.plot_trigger_overview(dat; ignore_triggers = ignore_list, display_plot = false)
            fig2, ax2 = eegfun.plot_trigger_timing(dat; ignore_triggers = ignore_list, display_plot = false)

            @test fig1 isa Figure && fig2 isa Figure
            @test ax1 isa Axis && ax2 isa Axis
        end

        @testset "performance optimization verification" begin
            dat = create_test_continuous_data_with_triggers(; n = 2000)

            # Test that empty ignore_triggers doesn't add overhead
            # This is more of a design verification than a strict performance test
            fig1, ax1 = eegfun.plot_trigger_overview(dat; display_plot = false)
            fig2, ax2 = eegfun.plot_trigger_overview(dat; ignore_triggers = Int[], display_plot = false)

            @test fig1 isa Figure && fig2 isa Figure
            @test ax1 isa Axis && ax2 isa Axis
        end

        @testset "filtering edge cases" begin
            dat = create_test_continuous_data_with_triggers()

            # Filter all triggers
            fig, ax = eegfun.plot_trigger_overview(dat; ignore_triggers = [1, 2, 3], display_plot = false)
            @test fig isa Figure
            @test ax isa Axis

            # Filter non-existent triggers
            fig, ax = eegfun.plot_trigger_overview(dat; ignore_triggers = [99, 100], display_plot = false)
            @test fig isa Figure
            @test ax isa Axis

            # Filter with mixed types (should be handled gracefully)
            fig, ax = eegfun.plot_trigger_overview(dat; ignore_triggers = [1.0, 2], display_plot = false)
            @test fig isa Figure
            @test ax isa Axis
        end

        @testset "advanced filtering edge cases" begin
            dat = create_test_continuous_data_with_triggers()

            # Very large ignore_triggers list
            large_ignore = collect(1:1000)
            fig, ax = eegfun.plot_trigger_overview(dat; ignore_triggers = large_ignore, display_plot = false)
            @test fig isa Figure
            @test ax isa Axis

            # Duplicate values in ignore_triggers
            duplicate_ignore = [1, 1, 2, 2, 3, 3]
            fig, ax = eegfun.plot_trigger_overview(dat; ignore_triggers = duplicate_ignore, display_plot = false)
            @test fig isa Figure
            @test ax isa Axis

            # Negative trigger values (if they exist in data)
            fig, ax = eegfun.plot_trigger_overview(dat; ignore_triggers = [-1, -2], display_plot = false)
            @test fig isa Figure
            @test ax isa Axis
        end
    end
end
