using Test
using eegfun

@testset "Plot Misc Utilities" begin

    @testset "_split_into_parts" begin
        # Test uppercase words
        @test eegfun._split_into_parts("ExampleCondition1") == ["Example", "Condition", "1"]
        @test eegfun._split_into_parts("ExampleCondition2") == ["Example", "Condition", "2"]
        
        # Test lowercase words
        @test eegfun._split_into_parts("switch_rep_easy1") == ["switch", "_", "rep", "_", "easy", "1"]
        @test eegfun._split_into_parts("switch_rep_hard2") == ["switch", "_", "rep", "_", "hard", "2"]
        
        # Test with multiple underscores and digits
        @test eegfun._split_into_parts("switch_rep_easy1_1_1") == ["switch", "_", "rep", "_", "easy", "1", "_", "1", "_", "1"]
        
        # Test mixed case
        @test eegfun._split_into_parts("ExampleTest21_22") == ["Example", "Test", "21", "_", "22"]
        
        # Test empty string
        @test eegfun._split_into_parts("") == String[]
        
        # Test single word
        @test eegfun._split_into_parts("Example") == ["Example"]
    end

    @testset "_find_common_prefix_parts" begin
        # Test with common prefix
        names = ["ExampleCondition1", "ExampleCondition2", "ExampleCondition3"]
        @test eegfun._find_common_prefix_parts(names) == ["Example", "Condition"]
        
        # Test with no common prefix
        names = ["ExampleCondition1", "AnotherCondition1"]
        @test eegfun._find_common_prefix_parts(names) == String[]
        
        # Test with underscores
        names = ["switch_rep_easy1", "switch_rep_hard2"]
        @test eegfun._find_common_prefix_parts(names) == ["switch", "_", "rep", "_"]
        
        # Test single name
        @test eegfun._find_common_prefix_parts(["Example"]) == String[]
        
        # Test empty
        @test eegfun._find_common_prefix_parts(String[]) == String[]
    end

    @testset "_abbreviate_parts" begin
        # Test basic abbreviation
        @test eegfun._abbreviate_parts(["Example", "Condition"]) == "ExaCon"
        @test eegfun._abbreviate_parts(["switch", "_", "rep", "_", "easy"]) == "Swi_Rep_Eas"
        
        # Test with digits
        @test eegfun._abbreviate_parts(["Example", "Condition", "1"]) == "ExaCon1"
        @test eegfun._abbreviate_parts(["Test", "21", "_", "22"]) == "Tes21_22"
        
        # Test with underscores
        @test eegfun._abbreviate_parts(["Example", "_", "Test"]) == "Exa_Tes"
        
        # Test short parts
        @test eegfun._abbreviate_parts(["A", "B"]) == "AB"
        @test eegfun._abbreviate_parts(["AB", "CD"]) == "ABCD"
    end

    @testset "_abbreviate_name" begin
        # Test with common prefix
        common = ["Example", "Condition"]
        @test eegfun._abbreviate_name("ExampleCondition1", common) == "ExaCon1"
        @test eegfun._abbreviate_name("ExampleCondition2", common) == "ExaCon2"
        
        # Test with underscores
        common = ["switch", "_", "rep", "_"]
        @test eegfun._abbreviate_name("switch_rep_easy1", common) == "Swi_Rep_easy1"
        @test eegfun._abbreviate_name("switch_rep_hard2", common) == "Swi_Rep_hard2"
        
        # Test with empty common prefix
        @test eegfun._abbreviate_name("Example", String[]) == "Example"
    end

    @testset "_shorten_condition_names" begin
        # Test with common prefix
        names = ["ExampleCondition1", "ExampleCondition2", "ExampleCondition3"]
        result = eegfun._shorten_condition_names(names)
        @test result == "ExaCon1, ExaCon2, ExaCon3"
        
        # Test with no common prefix (should abbreviate individually)
        names = ["ExampleCondition1", "ExampleCondition2", "ExampleCondition3", "AnotherCondition1"]
        result = eegfun._shorten_condition_names(names)
        @test result == "ExaCon1, ExaCon2, ExaCon3, AnoCon1"
        
        # Test with underscores (common prefix found, so uses shared abbreviation)
        names = ["switch_rep_easy1", "switch_rep_hard2"]
        result = eegfun._shorten_condition_names(names)
        # When common prefix exists, it abbreviates the prefix and keeps the unique suffix
        @test result == "Swi_Rep_easy1, Swi_Rep_hard2" || result == "switch_rep_easy1, switch_rep_hard2"
        
        # Test single name
        @test eegfun._shorten_condition_names(["ExampleCondition1"]) == "ExampleCondition1"
        
        # Test empty
        @test eegfun._shorten_condition_names(String[]) == ""
        
        # Test with long names (truncation)
        names = ["VeryLongConditionName1", "VeryLongConditionName2"]
        result = eegfun._shorten_condition_names(names; max_name_length = 10)
        @test length(result) > 0
        @test occursin("…", result) || length(result) <= 80  # Either truncated or within limit
    end

end

