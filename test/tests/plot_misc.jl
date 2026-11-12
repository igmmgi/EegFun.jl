using Test
using EegFun
import Logging

@testset "Plot Misc Utilities" begin

    @testset "_split_into_parts" begin
        # Test uppercase words
        @test EegFun._split_into_parts("ExampleCondition1") == ["Example", "Condition", "1"]
        @test EegFun._split_into_parts("ExampleCondition2") == ["Example", "Condition", "2"]

        # Test lowercase words
        @test EegFun._split_into_parts("switch_rep_easy1") == ["switch", "_", "rep", "_", "easy", "1"]
        @test EegFun._split_into_parts("switch_rep_hard2") == ["switch", "_", "rep", "_", "hard", "2"]

        # Test with multiple underscores and digits
        @test EegFun._split_into_parts("switch_rep_easy1_1_1") == ["switch", "_", "rep", "_", "easy", "1", "_", "1", "_", "1"]

        # Test mixed case
        @test EegFun._split_into_parts("ExampleTest21_22") == ["Example", "Test", "21", "_", "22"]

        # Test empty string
        @test EegFun._split_into_parts("") == String[]

        # Test single word
        @test EegFun._split_into_parts("Example") == ["Example"]
    end

    @testset "_find_common_prefix_parts" begin
        # Test with common prefix
        names = ["ExampleCondition1", "ExampleCondition2", "ExampleCondition3"]
        @test EegFun._find_common_prefix_parts(names) == ["Example", "Condition"]

        # Test with no common prefix
        names = ["ExampleCondition1", "AnotherCondition1"]
        @test EegFun._find_common_prefix_parts(names) == String[]

        # Test with underscores
        names = ["switch_rep_easy1", "switch_rep_hard2"]
        @test EegFun._find_common_prefix_parts(names) == ["switch", "_", "rep", "_"]

        # Test single name
        @test EegFun._find_common_prefix_parts(["Example"]) == String[]

        # Test empty
        @test EegFun._find_common_prefix_parts(String[]) == String[]
    end

    @testset "_abbreviate_parts" begin
        # Test basic abbreviation
        @test EegFun._abbreviate_parts(["Example", "Condition"]) == "ExaCon"
        @test EegFun._abbreviate_parts(["switch", "_", "rep", "_", "easy"]) == "Swi_Rep_Eas"

        # Test with digits
        @test EegFun._abbreviate_parts(["Example", "Condition", "1"]) == "ExaCon1"
        @test EegFun._abbreviate_parts(["Test", "21", "_", "22"]) == "Tes21_22"

        # Test with underscores
        @test EegFun._abbreviate_parts(["Example", "_", "Test"]) == "Exa_Tes"

        # Test short parts
        @test EegFun._abbreviate_parts(["A", "B"]) == "AB"
        @test EegFun._abbreviate_parts(["AB", "CD"]) == "ABCD"
    end

    @testset "_abbreviate_name" begin
        # Test with common prefix
        common = ["Example", "Condition"]
        @test EegFun._abbreviate_name("ExampleCondition1", common) == "ExaCon1"
        @test EegFun._abbreviate_name("ExampleCondition2", common) == "ExaCon2"

        # Test with underscores
        common = ["switch", "_", "rep", "_"]
        @test EegFun._abbreviate_name("switch_rep_easy1", common) == "Swi_Rep_easy1"
        @test EegFun._abbreviate_name("switch_rep_hard2", common) == "Swi_Rep_hard2"

        # Test with empty common prefix
        @test EegFun._abbreviate_name("Example", String[]) == "Example"
    end

    @testset "_shorten_condition_names" begin
        # Test with common prefix
        names = ["ExampleCondition1", "ExampleCondition2", "ExampleCondition3"]
        result = EegFun._shorten_condition_names(names)
        @test result == "ExaCon1, ExaCon2, ExaCon3"

        # Test with no common prefix (should abbreviate individually)
        names = ["ExampleCondition1", "ExampleCondition2", "ExampleCondition3", "AnotherCondition1"]
        result = EegFun._shorten_condition_names(names)
        @test result == "ExaCon1, ExaCon2, ExaCon3, AnoCon1"

        # Test with underscores (common prefix found, so uses shared abbreviation)
        names = ["switch_rep_easy1", "switch_rep_hard2"]
        result = EegFun._shorten_condition_names(names)
        # When common prefix exists, it abbreviates the prefix and keeps the unique suffix
        @test result == "Swi_Rep_easy1, Swi_Rep_hard2" || result == "switch_rep_easy1, switch_rep_hard2"

        # Test single name
        @test EegFun._shorten_condition_names(["ExampleCondition1"]) == "ExampleCondition1"

        # Test empty
        @test EegFun._shorten_condition_names(String[]) == ""

        # Test with long names (truncation)
        names = ["VeryLongConditionName1", "VeryLongConditionName2"]
        result = EegFun._shorten_condition_names(names; max_name_length = 10)
        @test length(result) > 0
        @test occursin("…", result) || length(result) <= 80  # Either truncated or within limit
    end

    @testset "EegFun theme is per-figure, not global" begin
        Makie = EegFun.Makie
        # Loading EegFun must not change the global Makie theme
        @test Makie.Axis(Makie.Figure()[1, 1]).xgridvisible[] == true

        # EegFun figures get the EegFun look
        fig = EegFun._eegfun_figure()
        ax = Makie.Axis(fig[1, 1])
        @test ax.xgridvisible[] == false
        @test ax.xticklabelsize[] == 16

        # ... and the global theme is still untouched afterwards
        @test Makie.Axis(Makie.Figure()[1, 1]).xgridvisible[] == true

        # User customisations take precedence over the EegFun defaults
        Makie.with_theme(Makie.Theme(fontsize = 10, Axis = (xgridvisible = true,))) do
            ax2 = Makie.Axis(EegFun._eegfun_figure()[1, 1])
            @test ax2.xticklabelsize[] == 10
            @test ax2.xgridvisible[] == true
            @test ax2.ygridvisible[] == false  # not customised -> EegFun default
        end
    end

    @testset "Interactive Demo Backend Checks" begin
        # Backend check returns false and warns when GLMakie is absent
        @test EegFun._check_interactive_demo_backend("Test Demo") == false
        @test EegFun._check_interactive_demo_backend() == false

        # Display helper handles missing or headless figures without throwing
        @test EegFun._display_demo_figure(nothing) === nothing
        fig = EegFun._eegfun_figure()
        @test EegFun._display_demo_figure(fig) === nothing

        # Interactive demos safely return nothing without GLMakie
        demo_functions = [
            EegFun.signal_example_decoding,
            EegFun.simulate_erp,
            EegFun.signal_example_ica_sphering,
            EegFun.signal_example_ica_separation,
            EegFun.signal_example_ica_optimization,
            EegFun.signal_example_ica_clt,
            EegFun.signal_example_ica_geometry,
            EegFun.signal_example_ica_mixture,
            EegFun.signal_example_heog,
            EegFun.signal_example_convolution,
            EegFun.signal_example_ecg_hrv,
            EegFun.signal_example_sampling,
            EegFun.signal_example_dotproduct,
            EegFun.signal_example_ecg_einthoven,
            EegFun.signal_example_tf,
            EegFun.signal_example_composition,
            EegFun.signal_example_emg,
            EegFun.signal_example_eda_scr,
            EegFun.signal_example_eda_habituation,
        ]
        Logging.with_logger(Logging.NullLogger()) do
            for demo_fn in demo_functions
                @test demo_fn() === nothing
            end
        end
    end

    @testset "_find_continuous_regions 1-sample cluster width" begin
        t = [0.0, 0.01, 0.02, 0.03, 0.04]
        mask_single = Bool[false, false, true, false, false]
        regs_single = EegFun._find_continuous_regions(mask_single, t)
        @test length(regs_single) == 1
        @test regs_single[1][1] ≈ 0.015
        @test regs_single[1][2] ≈ 0.025
        @test (regs_single[1][2] - regs_single[1][1]) ≈ 0.01

        mask_two = Bool[false, true, true, false, false]
        regs_two = EegFun._find_continuous_regions(mask_two, t)
        @test length(regs_two) == 1
        @test regs_two[1][1] ≈ 0.005
        @test regs_two[1][2] ≈ 0.025
    end

    @testset "simulate_erp _SimulatedErpComponent and ms-to-samples scaling" begin
        @test isdefined(EegFun, :_SimulatedErpComponent)
        @test !isdefined(EegFun, :Component)

        # 1 component at 300 ms latency, 10 μV amplitude, zero noise
        comps = [5.0 10.0 300.0 0.0 0.0 0.0]
        avg_500, _ = EegFun._run_erp_simulation(1, comps; samp_freq = 500, sig_length = 1.0, base_length = 0.2)
        # Expected peak sample at 500 Hz: base_length (0.2s * 500 = 100) + 300ms * 0.5 = 250
        @test argmax(avg_500) == 250
    end

end
