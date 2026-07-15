using Test
using EegFun
using DataFrames

@testset "Pipeline Utilities" begin

    # ────────────────────────────────────────────────────────────
    # 1. Section formatting helpers
    # ────────────────────────────────────────────────────────────
    @testset "section formatting" begin
        s = EegFun.section("Test Section")
        @test contains(s, "Test Section")
        @test contains(s, "-")

        sub = EegFun.subsection("Sub Test")
        @test contains(sub, "Sub Test")

        subsub = EegFun.subsubsection("Sub Sub")
        @test contains(subsub, "Sub Sub")
        @test startswith(subsub, "\n# ")
    end

    @testset "_center_title" begin
        title = EegFun._center_title("Hello", 40)
        @test length(title) == 40
        @test contains(title, "Hello")
    end

    @testset "_flag_symbol" begin
        s = EegFun._flag_symbol("is_extreme_value", 200)
        @test s == :is_extreme_value_200

        s2 = EegFun._flag_symbol("is_artifact_value", 100.0)
        @test s2 == Symbol("is_artifact_value_100.0")
    end

    # ────────────────────────────────────────────────────────────
    # 2. check_files_exist
    # ────────────────────────────────────────────────────────────
    @testset "check_files_exist" begin
        # Create temp files
        dir = mktempdir()
        f1 = joinpath(dir, "file1.txt")
        f2 = joinpath(dir, "file2.txt")
        write(f1, "test")
        write(f2, "test")

        @test EegFun.check_files_exist([f1, f2]) == true
        @test EegFun.check_files_exist([f1, "/nonexistent/file.txt"]) == false
        @test EegFun.check_files_exist(String[]) == true
    end

    # ────────────────────────────────────────────────────────────
    # 3. ContinuousRepairInfo
    # ────────────────────────────────────────────────────────────
    @testset "ContinuousRepairInfo" begin
        info = EegFun.create_continuous_repair_info(:neighbor_interpolation; name = "test_repair")

        @test info isa EegFun.ContinuousRepairInfo
        @test info.name == "test_repair"
        @test info.method == :neighbor_interpolation
        @test isempty(info.repaired)
        @test isempty(info.skipped)

        # Show method
        io = IOBuffer()
        show(io, info)
        output = String(take!(io))
        @test contains(output, "ContinuousRepairInfo")
        @test contains(output, "test_repair")
    end

    # ────────────────────────────────────────────────────────────
    # 4. ArtifactInfo
    # ────────────────────────────────────────────────────────────
    @testset "ArtifactInfo" begin
        info = EegFun.ArtifactInfo()

        @test info isa EegFun.ArtifactInfo
        @test isempty(info.continuous_repairs)
        @test isempty(info.epoch_rejections)
        @test isnothing(info.ica_components)

        # Show method
        io = IOBuffer()
        show(io, info)
        output = String(take!(io))
        @test contains(output, "ArtifactInfo")
    end

    # ────────────────────────────────────────────────────────────
    # 5. ChannelRepairInfo
    # ────────────────────────────────────────────────────────────
    @testset "ChannelRepairInfo" begin
        info = EegFun.ChannelRepairInfo()

        @test isnothing(info.continuous)
        @test isempty(info.epochs)

        io = IOBuffer()
        show(io, info)
        output = String(take!(io))
        @test contains(output, "ChannelRepairInfo")
    end

    # ────────────────────────────────────────────────────────────
    # 6. Pipeline template generation
    # ────────────────────────────────────────────────────────────
    @testset "generate_pipeline_template" begin
        dir = mktempdir()
        output_file = joinpath(dir, "test_pipeline.jl")

        result = EegFun.generate_pipeline_template(
            output_file,
            "test_preprocess";
            options = EegFun.PipelineTemplateOptions(num_steps = 2, subsections_per_step = 1),
        )

        @test isfile(result)
        content = read(result, String)
        @test contains(content, "test_preprocess")
        @test contains(content, "using EegFun")
        @test contains(content, "Step 1")
        @test contains(content, "Step 2")
    end

    # ────────────────────────────────────────────────────────────
    # 7. preprocess error on missing config
    # ────────────────────────────────────────────────────────────
    @testset "preprocess - missing config" begin
        @test_throws Exception EegFun.preprocess("/nonexistent/config.toml")
    end

    # ────────────────────────────────────────────────────────────
    # 8. preprocess end-to-end integration test
    # ────────────────────────────────────────────────────────────
    @testset "preprocess - end-to-end" begin
        data_dir = joinpath(dirname(dirname(@__DIR__)), "resources", "data")
        test_file = joinpath(data_dir, "bdf", "example1.bdf")

        if isfile(test_file)
            # Create a temporary directory for the pipeline config and output
            mktempdir() do temp_dir
                # Create epoch condition file
                epoch_file = joinpath(temp_dir, "epochs.toml")
                write(
                    epoch_file,
                    """
  [epochs]
  [[epochs.conditions]]
  name = "Test Condition"
  trigger_sequences = [1]
  """,
                )

                # Create pipeline config file
                pipeline_file = joinpath(temp_dir, "pipeline.toml")
                write(
                    pipeline_file,
                    """
[files.input]
directory = "$(joinpath(data_dir, "bdf"))"
raw_data_files = "example1\\\\.bdf"
layout_file = "biosemi72.csv"
epoch_condition_file = "epochs.toml"

[files.output]
directory = "output"
save_continuous_data_original = false
save_continuous_data_cleaned = false

[preprocess]
epoch_start = -0.1
epoch_end = 0.5
""",
                )

                # Run the pipeline (this shouldn't throw any errors)
                # Redirect stdout/stderr to avoid spamming the test output
                redirect_stdio(stdout = devnull, stderr = devnull) do
                    EegFun.preprocess(pipeline_file; base_dir = temp_dir)
                end

                # Check that the output directory was created and contains the expected files
                out_dir = joinpath(temp_dir, "output")
                @test isdir(out_dir)
                @test isfile(joinpath(out_dir, "example1_epochs_good.jld2"))
                @test isfile(joinpath(out_dir, "example1_erps_good.jld2"))
                @test isfile(joinpath(out_dir, "example1_artifact_info.jld2"))
            end
        else
            EegFun.@minimal_warning "Test file not found: $test_file (skipping end-to-end pipeline test)"
        end
    end

    # ────────────────────────────────────────────────────────────
    # Cleanup: remove any preprocess_log files generated during tests
    # ────────────────────────────────────────────────────────────
    for f in filter(f -> startswith(f, "preprocess_log_") && endswith(f, ".txt"), readdir())
        rm(f, force = true)
    end

end # @testset "Pipeline Utilities"
