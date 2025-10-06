using Test
using DataFrames
using eegfun
using OrderedCollections

@testset "triggers" begin

    @testset "_clean_triggers" begin
        @testset "basic onset detection" begin
            # Sustained trigger signal
            input = [0, 1, 1, 1, 0, 0, 2, 2, 0, 0]
            expected = [0, 1, 0, 0, 0, 0, 2, 0, 0, 0]
            result = eegfun._clean_triggers(input)
            @test result == expected
        end

        @testset "single trigger events" begin
            # Already clean trigger data
            input = [0, 1, 0, 0, 2, 0, 0, 3, 0, 0]
            expected = [0, 1, 0, 0, 2, 0, 0, 3, 0, 0]
            result = eegfun._clean_triggers(input)
            @test result == expected
        end

        @testset "no triggers" begin
            # All zeros
            input = [0, 0, 0, 0, 0]
            expected = [0, 0, 0, 0, 0]
            result = eegfun._clean_triggers(input)
            @test result == expected
        end

        @testset "multiple consecutive different triggers" begin
            # Different trigger values in sequence
            input = [0, 1, 2, 3, 0, 0]
            expected = [0, 1, 2, 3, 0, 0]
            result = eegfun._clean_triggers(input)
            @test result == expected
        end

        @testset "edge cases" begin
            # Single element
            @test eegfun._clean_triggers([1]) == [1]
            @test eegfun._clean_triggers([0]) == [0]

            # Start with trigger
            input = [1, 1, 0, 0]
            expected = [1, 0, 0, 0]
            @test eegfun._clean_triggers(input) == expected

            # End with trigger
            input = [0, 0, 1, 1]
            expected = [0, 0, 1, 0]
            @test eegfun._clean_triggers(input) == expected
        end

        @testset "boundary conditions" begin
            # Single sample
            @test eegfun._clean_triggers([1]) == [1]
            @test eegfun._clean_triggers([0]) == [0]

            # Two samples
            @test eegfun._clean_triggers([0, 1]) == [0, 1]
            @test eegfun._clean_triggers([1, 1]) == [1, 0]
            @test eegfun._clean_triggers([1, 0]) == [1, 0]
        end

        @testset "type handling" begin
            # Different integer types
            triggers_int8 = Int8[0, 1, 1, 0, 2]
            triggers_int16 = Int16[0, 1, 1, 0, 2]
            triggers_int32 = Int32[0, 1, 1, 0, 2]

            # All should work with cleaning
            @test eegfun._clean_triggers(triggers_int8) == [0, 1, 0, 0, 2]
            @test eegfun._clean_triggers(triggers_int16) == [0, 1, 0, 0, 2]
            @test eegfun._clean_triggers(triggers_int32) == [0, 1, 0, 0, 2]
        end
    end

    @testset "trigger_count" begin
        @testset "ContinuousData trigger counting" begin
            # Create test data
            triggers = [0, 1, 1, 0, 2, 2, 2, 0, 1, 0]
            time = collect(0:9) ./ 100.0
            df = DataFrame(time = time, triggers = triggers, A = zeros(10), B = zeros(10))
            layout = eegfun.Layout(DataFrame(label = [:A, :B], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)
            dat = eegfun.ContinuousData(df, layout, 100, eegfun.AnalysisInfo())

            # Test with printing disabled
            result = eegfun.trigger_count(dat, print_table = false)
            @test isa(result, DataFrame)
            @test size(result, 1) == 2  # Two unique non-zero triggers
            @test result.trigger == [1, 2]
            @test result.count == [3, 3]  # Raw counts before cleaning
            @test !("triggers_info" in names(result))  # No triggers_info column when not present
        end

        @testset "ContinuousData with triggers_info" begin
            # Create test data with triggers_info
            triggers = [0, 1, 1, 0, 2, 2, 2, 0, 1, 0]
            triggers_info = ["", "S 1", "", "", "S 2", "", "", "", "S 1", ""]
            time = collect(0:9) ./ 100.0
            df =
                DataFrame(time = time, triggers = triggers, triggers_info = triggers_info, A = zeros(10), B = zeros(10))
            layout = eegfun.Layout(DataFrame(label = [:A, :B], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)
            dat = eegfun.ContinuousData(df, layout, 100, eegfun.AnalysisInfo())

            # Test with printing disabled
            result = eegfun.trigger_count(dat, print_table = false)
            @test isa(result, DataFrame)
            @test size(result, 1) == 2  # Two unique non-zero triggers
            @test result.trigger == [1, 2]
            @test result.count == [3, 3]  # Raw counts before cleaning
            @test "triggers_info" in names(result)  # triggers_info column should be present
            @test result.triggers_info == ["S 1", "S 2"]  # Should show the info for each trigger
        end

        @testset "empty triggers" begin
            # Create test data with no triggers
            triggers = [0, 0, 0, 0, 0]
            time = collect(0:4) ./ 100.0
            df = DataFrame(time = time, triggers = triggers, A = zeros(5), B = zeros(5))
            layout = eegfun.Layout(DataFrame(label = [:A, :B], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)
            dat = eegfun.ContinuousData(df, layout, 100, eegfun.AnalysisInfo())

            result = eegfun.trigger_count(dat, print_table = false)
            @test isa(result, DataFrame)
            @test size(result, 1) == 0  # No triggers
        end

        @testset "data structure integration" begin
            # Test with ContinuousData structure
            triggers = [0, 1, 0, 2, 0, 1, 0]
            time = collect(0:6) ./ 100.0
            df = DataFrame(time = time, triggers = triggers, channel1 = randn(7), channel2 = randn(7))
            layout = eegfun.Layout(
                DataFrame(label = [:channel1, :channel2], inc = [0.0, 90.0], azi = [0.0, 0.0]),
                nothing,
                nothing,
            )
            dat = eegfun.ContinuousData(df, layout, 100, eegfun.AnalysisInfo())

            # Test trigger counting
            count_df = eegfun.trigger_count(dat, print_table = false)
            @test size(count_df, 1) == 2
            @test count_df.trigger == [1, 2]
            @test count_df.count == [2, 1]
        end
    end

    @testset "_trigger_count_impl unified helper function" begin
        @testset "single dataset counting" begin
            triggers = [0, 1, 1, 0, 2, 2, 2, 0, 1, 0]
            result = eegfun._trigger_count_impl([triggers], ["count"], print_table = false)
            @test isa(result, DataFrame)
            @test size(result, 1) == 2
            @test result.trigger == [1, 2]
            @test result.count == [3, 3]
        end

        @testset "multiple dataset counting" begin
            raw_triggers = [0, 1, 1, 1, 0, 2, 2, 0, 0]
            cleaned_triggers = eegfun._clean_triggers(raw_triggers)
            result = eegfun._trigger_count_impl(
                [raw_triggers, cleaned_triggers],
                ["raw_count", "cleaned_count"],
                print_table = false,
            )
            @test isa(result, DataFrame)
            @test size(result, 1) == 2
            @test result.trigger == [1, 2]
            @test result.raw_count == [3, 2]      # Raw counts
            @test result.cleaned_count == [1, 1]  # Cleaned counts (onset only)
        end

        @testset "no triggers" begin
            triggers = [0, 0, 0, 0]
            result = eegfun._trigger_count_impl([triggers], ["count"], print_table = false)
            @test isa(result, DataFrame)
            @test size(result, 1) == 0
            @test names(result) == ["trigger", "count"]
        end

        @testset "with triggers_info" begin
            triggers = [0, 1, 1, 0, 2, 2, 0]
            triggers_info = ["", "S 1", "", "", "S 2", "", ""]
            result =
                eegfun._trigger_count_impl([triggers], ["count"], print_table = false, triggers_info = triggers_info)
            @test isa(result, DataFrame)
            @test size(result, 1) == 2
            @test result.trigger == [1, 2]
            @test result.count == [2, 2]
            @test "triggers_info" in names(result)
            @test result.triggers_info == ["S 1", "S 2"]
        end

        @testset "custom headers and notes" begin
            triggers = [0, 1, 0, 2, 0]
            result = eegfun._trigger_count_impl(
                [triggers],
                ["custom_count"],
                print_table = false,
                title = "Custom Title",
                headers = ["ID", "Custom Count"],
                note = "Test note",
            )
            @test isa(result, DataFrame)
            @test names(result) == ["trigger", "custom_count"]
            @test size(result, 1) == 2
        end
    end




    @testset "search_sequence" begin
        @testset "basic functionality" begin
            # Test basic trigger finding
            triggers = [0, 1, 0, 1, 0, 2, 0]
            indices = eegfun.search_sequence(triggers, 1)
            @test indices == [2, 4]  # Both occurrences of trigger 1

            # Test with onset detection built-in
            sustained_triggers = [0, 1, 1, 1, 0, 1, 0]
            indices = eegfun.search_sequence(sustained_triggers, 1)
            @test indices == [2, 6]  # Only onset positions

            # Test non-existent trigger
            indices = eegfun.search_sequence(triggers, 99)
            @test indices == Int[]
        end

        @testset "edge cases" begin
            # Empty array
            @test eegfun.search_sequence(Int[], 1) == Int[]

            # Single element
            @test eegfun.search_sequence([1], 1) == [1]
            @test eegfun.search_sequence([0], 1) == Int[]

            # All same trigger
            @test eegfun.search_sequence([1, 1, 1], 1) == [1]  # Only first onset
        end

        @testset "onset detection correctness" begin
            # Test the specific case that was failing
            @test eegfun.search_sequence([1,2,3,1,0,1], 1) == [1, 4, 6]  
            
            # Test various onset patterns
            @test eegfun.search_sequence([0,1,1,0,2,2,0,1], 1) == [2, 8]  # Two onsets of trigger 1
            @test eegfun.search_sequence([0,1,1,0,2,2,0,1], 2) == [5]     # One onset of trigger 2
            
            # Test consecutive different triggers
            @test eegfun.search_sequence([0,1,2,3,0,1,2,3], 1) == [2, 6]  # Two onsets
            @test eegfun.search_sequence([0,1,2,3,0,1,2,3], 2) == [3, 7]  # Two onsets
            @test eegfun.search_sequence([0,1,2,3,0,1,2,3], 3) == [4, 8]  # Two onsets
            
            # Test sustained triggers
            @test eegfun.search_sequence([0,1,1,1,0,2,2,0], 1) == [2]  # Only onset, not sustained
            @test eegfun.search_sequence([0,1,1,1,0,2,2,0], 2) == [6]  # Only onset, not sustained
            
            # Test mixed patterns
            @test eegfun.search_sequence([1,0,1,0,1,0], 1) == [1, 3, 5]  # All onsets
            @test eegfun.search_sequence([1,1,0,1,1,0], 1) == [1, 4]    # Only first of each sustained block
        end

        @testset "boundary conditions" begin
            # Start with trigger
            @test eegfun.search_sequence([1,0,1,0], 1) == [1, 3]
            
            # End with trigger
            @test eegfun.search_sequence([0,1,0,1], 1) == [2, 4]
            
            # Start and end with same trigger
            @test eegfun.search_sequence([1,0,0,1], 1) == [1, 4]
            
            # All zeros
            @test eegfun.search_sequence([0,0,0,0], 1) == Int[]
            
            # Single non-zero element
            @test eegfun.search_sequence([5], 5) == [1]
            @test eegfun.search_sequence([5], 1) == Int[]
        end

        @testset "complex patterns" begin
            # Alternating pattern
            @test eegfun.search_sequence([1,2,1,2,1,2], 1) == [1, 3, 5]
            @test eegfun.search_sequence([1,2,1,2,1,2], 2) == [2, 4, 6]
            
            # Multiple different triggers
            @test eegfun.search_sequence([1,2,3,1,2,3,1,2,3], 1) == [1, 4, 7]
            @test eegfun.search_sequence([1,2,3,1,2,3,1,2,3], 2) == [2, 5, 8]
            @test eegfun.search_sequence([1,2,3,1,2,3,1,2,3], 3) == [3, 6, 9]
            
            # Sustained blocks of different lengths
            @test eegfun.search_sequence([0,1,1,0,2,2,2,0,3], 1) == [2]
            @test eegfun.search_sequence([0,1,1,0,2,2,2,0,3], 2) == [5]
            @test eegfun.search_sequence([0,1,1,0,2,2,2,0,3], 3) == [9]
        end
    end

    @testset "search_sequence (ranges)" begin
        @testset "basic range searching" begin
            triggers = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

            # Single range
            indices = eegfun.search_sequence(triggers, [3:5])
            @test indices == [3, 4, 5]  # Positions of triggers 3, 4, 5

            # Multiple ranges
            indices = eegfun.search_sequence(triggers, [1:3, 8:10])
            @test indices == [1, 2, 3, 8, 9, 10]

            # Overlapping ranges (no duplicates)
            indices = eegfun.search_sequence(triggers, [2:4, 3:5])
            @test indices == [2, 3, 4, 5]
        end

        @testset "edge cases" begin
            triggers = [1, 2, 3, 4, 5]

            # Empty ranges
            indices = eegfun.search_sequence(triggers, UnitRange{Int}[])
            @test indices == Int[]

            # Range outside data
            indices = eegfun.search_sequence(triggers, [10:15])
            @test indices == Int[]

            # Single value range
            indices = eegfun.search_sequence(triggers, [3:3])
            @test indices == [3]
        end
    end

    @testset "search_sequence (single sequence)" begin
        @testset "exact sequences" begin
            # More realistic trigger data with proper onsets
            triggers = [0, 1, 2, 3, 0, 1, 2, 4, 0, 1, 2, 3, 0]

            # Find exact sequence [1, 2, 3]
            sequence = Vector{Union{Int,Symbol,UnitRange{Int}}}([1, 2, 3])
            indices = eegfun.search_sequence(triggers, sequence)
            @test indices == [2, 10]  # Two occurrences

            # Find sequence [1, 2, 4]
            sequence = Vector{Union{Int,Symbol,UnitRange{Int}}}([1, 2, 4])
            indices = eegfun.search_sequence(triggers, sequence)
            @test indices == [6]  # One occurrence
        end

        @testset "wildcard sequences" begin
            # More realistic trigger data with proper onsets
            triggers = [0, 1, 2, 3, 0, 1, 5, 3, 0, 1, 7, 3, 0]

            # Wildcard sequence [1, :any, 3]
            sequence = Vector{Union{Int,Symbol,UnitRange{Int}}}([1, :any, 3])
            indices = eegfun.search_sequence(triggers, sequence)
            @test indices == [2, 6, 10]  # All three occurrences
        end

        @testset "range sequences" begin
            # More realistic trigger data with proper onsets
            triggers = [0, 1, 2, 3, 0, 1, 4, 3, 0, 1, 5, 3, 0]

            # Range sequence [1, 2:5, 3]
            sequence = Vector{Union{Int,Symbol,UnitRange{Int}}}([1, 2:5, 3])
            indices = eegfun.search_sequence(triggers, sequence)
            @test indices == [2, 6, 10]  # All match since 2, 4, 5 are in range 2:5
        end

        @testset "single element sequences" begin
            triggers = [1, 2, 3, 1, 2, 3]

            # Single integer
            sequence = Vector{Union{Int,Symbol,UnitRange{Int}}}([2])
            indices = eegfun.search_sequence(triggers, sequence)
            @test indices == [2, 5]  # Both occurrences of trigger 2

            # Single range
            sequence = Vector{Union{Int,Symbol,UnitRange{Int}}}([2:3])
            indices = eegfun.search_sequence(triggers, sequence)
            @test indices == [2, 3, 5, 6]  # All occurrences of triggers 2 and 3
        end

        @testset "error handling" begin
            triggers = [1, 2, 3, 4, 5]

            # First element cannot be wildcard
            sequence1 = Vector{Union{Int,Symbol,UnitRange{Int}}}([:any, 2])
            @test_throws Exception eegfun.search_sequence(triggers, sequence1)

            # Single wildcard not supported
            sequence2 = Vector{Union{Int,Symbol,UnitRange{Int}}}([:any])
            @test_throws Exception eegfun.search_sequence(triggers, sequence2)

            # Unsupported trigger type
            sequence3 = Vector{Union{Int,Symbol,UnitRange{Int},String}}([1, "invalid"])
            @test_throws Exception eegfun.search_sequence(triggers, sequence3)
        end

        @testset "boundary conditions" begin
            triggers = [1, 2]

            # Sequence longer than array
            sequence = Vector{Union{Int,Symbol,UnitRange{Int}}}([1, 2, 3, 4])
            indices = eegfun.search_sequence(triggers, sequence)
            @test indices == Int[]

            # Sequence matches entire array
            sequence = Vector{Union{Int,Symbol,UnitRange{Int}}}([1, 2])
            indices = eegfun.search_sequence(triggers, sequence)
            @test indices == [1]
        end
    end

    @testset "search_sequence" begin
        @testset "multiple sequences (OR logic)" begin
            # More realistic trigger data with proper onsets (0 between sequences)
            triggers = [0, 1, 2, 3, 0, 1, 4, 3, 0, 1, 5, 1, 0]

            # Multiple exact sequences
            sequences = [
                Vector{Union{Int,Symbol,UnitRange{Int}}}([1, 2, 3]),
                Vector{Union{Int,Symbol,UnitRange{Int}}}([1, 4, 3]),
            ]
            indices = eegfun.search_sequence(triggers, sequences)
            @test sort(indices) == [2, 6]  # Both sequences found at positions 2 and 6

            # Mix of wildcards and exact
            sequences = [
                Vector{Union{Int,Symbol,UnitRange{Int}}}([1, :any, 3]),
                Vector{Union{Int,Symbol,UnitRange{Int}}}([1, 5, 1]),
            ]
            indices = eegfun.search_sequence(triggers, sequences)
            @test sort(indices) == [2, 6, 10]  # Three matches at positions 2, 6, 10
        end

        @testset "overlapping sequences" begin
            triggers = [1, 2, 3, 4, 5]

            # Overlapping sequences should return unique indices
            sequences = [
                Vector{Union{Int,Symbol,UnitRange{Int}}}([1, 2, 3]),
                Vector{Union{Int,Symbol,UnitRange{Int}}}([2, 3, 4]),
                Vector{Union{Int,Symbol,UnitRange{Int}}}([3, 4, 5]),
            ]
            indices = eegfun.search_sequence(triggers, sequences)
            @test sort(indices) == [1, 2, 3]  # Unique starting positions
        end

        @testset "no matches" begin
            triggers = [1, 2, 3, 4, 5]
            sequences = [
                Vector{Union{Int,Symbol,UnitRange{Int}}}([6, 7, 8]),
                Vector{Union{Int,Symbol,UnitRange{Int}}}([9, 10, 11]),
            ]
            indices = eegfun.search_sequence(triggers, sequences)
            @test indices == Int[]
        end

        @testset "complex mixed sequences" begin
            # More realistic trigger data with proper onsets
            triggers = [0, 1, 2, 3, 0, 1, 5, 3, 0, 2, 7, 8, 0, 1, 6, 3, 0]

            # Mix of ranges, wildcards, and exact values
            sequences = [
                Vector{Union{Int,Symbol,UnitRange{Int}}}([1, 2:5, 3]),      # Should match positions 2, 6
                Vector{Union{Int,Symbol,UnitRange{Int}}}([2, :any, 8]),     # Should match position 10
                Vector{Union{Int,Symbol,UnitRange{Int}}}([1, 6, 3]),         # Should match position 14
            ]
            indices = eegfun.search_sequence(triggers, sequences)
            @test sort(unique(indices)) == [2, 6, 10, 14]
        end
    end

    @testset "integration tests" begin
        @testset "realistic EEG workflow with sequences" begin
            # Create realistic test data with trigger sequences
            triggers = [0, 1, 2, 3, 0, 0, 1, 5, 3, 0, 1, 2, 4, 0, 0]
            time_data = collect(0:14) ./ 100.0

            # Test basic search functions
            @test eegfun.search_sequence(triggers, 1) == [2, 7, 11]  # All trigger 1 onsets

            # Test sequence searching
            seq_123 = eegfun.search_sequence(triggers, Vector{Union{Int,Symbol,UnitRange{Int}}}([1, 2, 3]))
            @test seq_123 == [2]  # Only one [1,2,3] sequence

            # Test wildcard sequences
            seq_1_any_3 =
                eegfun.search_sequence(triggers, Vector{Union{Int,Symbol,UnitRange{Int}}}([1, :any, 3]))
            @test sort(seq_1_any_3) == [2, 7]  # [1,2,3] and [1,5,3]

            # Test range sequences
            seq_1_range_3 =
                eegfun.search_sequence(triggers, Vector{Union{Int,Symbol,UnitRange{Int}}}([1, 2:5, 3]))
            @test sort(seq_1_range_3) == [2, 7]  # Both sequences match

            # Test multiple sequences
            sequences = [
                Vector{Union{Int,Symbol,UnitRange{Int}}}([1, 2, 3]),
                Vector{Union{Int,Symbol,UnitRange{Int}}}([1, 2, 4]),
            ]
            multiple_seqs = eegfun.search_sequence(triggers, sequences)
            @test sort(multiple_seqs) == [2, 11]  # Both sequences found
        end

        @testset "performance with larger data" begin
            # Test with larger dataset
            n_samples = 1000
            large_triggers = zeros(Int, n_samples)

            # Create some sequences
            large_triggers[100:102] = [1, 2, 3]
            large_triggers[200:202] = [1, 5, 3]
            large_triggers[300:302] = [2, 3, 4]
            large_triggers[400:401] = [1, 7]

            # Test search functions don't error and return reasonable results
            indices_1 = eegfun.search_sequence(large_triggers, 1)
            @test length(indices_1) >= 2  # At least 2 occurrences of trigger 1

            seq_indices =
                eegfun.search_sequence(large_triggers, Vector{Union{Int,Symbol,UnitRange{Int}}}([1, :any, 3]))
            @test length(seq_indices) >= 2  # At least [1,2,3] and [1,5,3]

            large_sequences = [
                Vector{Union{Int,Symbol,UnitRange{Int}}}([1, 2, 3]),
                Vector{Union{Int,Symbol,UnitRange{Int}}}([2, 3, 4]),
            ]
            multi_seq_indices = eegfun.search_sequence(large_triggers, large_sequences)
            @test length(multi_seq_indices) >= 2  # At least two different sequences
        end

        @testset "integration with ContinuousData" begin
            # Test that search functions work with real data structures
            triggers = [0, 1, 2, 3, 0, 1, 5, 3, 0, 0]
            time = collect(0:9) ./ 100.0
            df = DataFrame(time = time, triggers = triggers, channel = randn(10))
            layout = eegfun.Layout(DataFrame(label = [:channel], inc = [0.0], azi = [0.0]), nothing, nothing)
            dat = eegfun.ContinuousData(df, layout, 100, eegfun.AnalysisInfo())

            # Test that we can extract triggers and use search functions
            trigger_data = dat.data.triggers

            # Basic searching
            indices_1 = eegfun.search_sequence(trigger_data, 1)
            @test indices_1 == [2, 6]

            # Sequence searching
            seq_indices =
                eegfun.search_sequence(trigger_data, Vector{Union{Int,Symbol,UnitRange{Int}}}([1, :any, 3]))
            @test sort(seq_indices) == [2, 6]  # Both [1,2,3] and [1,5,3]

            # Verify trigger counting still works
            count_result = eegfun.trigger_count(dat, print_table = false)
            @test size(count_result, 1) == 4  # Four unique triggers: 1, 2, 3, 5
        end
    end

    @testset "error handling and edge cases" begin
        @testset "missing triggers column" begin
            # Test ContinuousData without triggers column
            df = DataFrame(time = [0.0, 0.1, 0.2], channel = [1.0, 2.0, 3.0])
            layout = eegfun.Layout(DataFrame(label = [:channel], inc = [0.0], azi = [0.0]), nothing, nothing)
            dat = eegfun.ContinuousData(df, layout, 100, eegfun.AnalysisInfo())

            @test_throws AssertionError eegfun.trigger_count(dat)
        end

        @testset "consistent output types" begin
            # Ensure consistent DataFrame output regardless of input
            triggers1 = [0, 1, 0, 2, 0]
            triggers2 = Int[]  # Empty

            result1 = eegfun._trigger_count_impl([triggers1], ["count"], print_table = false)
            result2 = eegfun._trigger_count_impl([triggers2], ["count"], print_table = false)

            @test isa(result1, DataFrame)
            @test isa(result2, DataFrame)
            @test names(result1) == names(result2)  # Same column structure
        end

        @testset "search function robustness" begin
            # Test with various edge cases
            empty_triggers = Int[]
            single_trigger = [1]
            zero_triggers = [0, 0, 0]

            # search_sequence should handle empty arrays
            @test eegfun.search_sequence(empty_triggers, 1) == Int[]
            @test eegfun.search_sequence(zero_triggers, 1) == Int[]
            @test eegfun.search_sequence(single_trigger, 1) == [1]

            # search_trigger_ranges should handle empty ranges
            @test eegfun.search_sequence(single_trigger, UnitRange{Int}[]) == Int[]

            # search_single_sequence should handle empty arrays
            sequence = Vector{Union{Int,Symbol,UnitRange{Int}}}([1, 2])
            @test eegfun.search_sequence(empty_triggers, sequence) == Int[]

            # search_sequences should handle empty sequence list
            empty_sequences = Vector{Vector{Union{Int,Symbol,UnitRange{Int}}}}()
            @test eegfun.search_sequence(single_trigger, empty_sequences) == Int[]
        end
    end
end
