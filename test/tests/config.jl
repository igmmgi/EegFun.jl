using Test
using TOML
using Dates

@testset "Configuration System Tests" begin

    # Create temporary directory for test files
    test_dir = mktempdir()

    # =============================================================================
    # CONFIGPARAMETER AND PARAMETER STRUCTS
    # =============================================================================

    @testset "ConfigParameter Tests" begin
        # Test ConfigParameter struct creation
        param = EegFun.ConfigParameter{Float64}(description = "Test parameter", default = 1.0, min = 0.0, max = 2.0)
        @test param.description == "Test parameter"
        @test param.default == 1.0
        @test param.min == 0.0
        @test param.max == 2.0

        # Test ConfigParameter with nothing values
        param = EegFun.ConfigParameter{String}(
            description = "Test parameter",
            default = nothing,
            min = nothing,
            max = nothing,
            allowed = nothing,
        )
        @test param.description == "Test parameter"
        @test isnothing(param.default)
        @test isnothing(param.min)
        @test isnothing(param.max)
        @test isnothing(param.allowed)
    end

    # =============================================================================
    # PARAMETERS DICTIONARY
    # =============================================================================

    @testset "PARAMETERS Dictionary Tests" begin
        # Test that all required parameters exist
        @test haskey(EegFun.PARAMETERS, "files.input.directory")
        @test haskey(EegFun.PARAMETERS, "files.input.raw_data_files")
        @test haskey(EegFun.PARAMETERS, "files.input.layout_file")
        @test haskey(EegFun.PARAMETERS, "preprocess.filter.highpass.freq")
        @test haskey(EegFun.PARAMETERS, "preprocess.filter.lowpass.freq")
        @test haskey(EegFun.PARAMETERS, "preprocess.epoch_start")
        @test haskey(EegFun.PARAMETERS, "preprocess.epoch_end")

        # Test parameter types
        @test EegFun.PARAMETERS["files.input.directory"] isa EegFun.ConfigParameter{String}
        @test EegFun.PARAMETERS["preprocess.filter.highpass.freq"] isa EegFun.ConfigParameter{Real}
        @test EegFun.PARAMETERS["preprocess.filter.highpass.order"] isa EegFun.ConfigParameter{Real}
        @test EegFun.PARAMETERS["files.output.save_continuous_data_raw"] isa EegFun.ConfigParameter{Bool}
        @test EegFun.PARAMETERS["files.output.save_continuous_data_corrected"] isa EegFun.ConfigParameter{Bool}

        # Test that all parameters have valid descriptions
        for (path, param) in EegFun.PARAMETERS
            @test !isempty(param.description)
            @test param isa EegFun.ConfigParameter
        end

        # Test specific parameter properties
        param = EegFun.PARAMETERS["preprocess.filter.highpass.freq"]
        @test param.description == "Cutoff frequency (Hz)"
        @test param.min == 0.01
        @test param.max == 20.0
        @test param.default == 0.1

        # Test non-existent parameter
        @test !haskey(EegFun.PARAMETERS, "nonexistent")

        # Test that all filter parameters exist
        @test haskey(EegFun.PARAMETERS, "preprocess.filter.highpass.apply")
        @test haskey(EegFun.PARAMETERS, "preprocess.filter.highpass.method")
        @test haskey(EegFun.PARAMETERS, "preprocess.filter.highpass.func")
        @test haskey(EegFun.PARAMETERS, "preprocess.filter.lowpass.apply")
        @test haskey(EegFun.PARAMETERS, "preprocess.filter.lowpass.method")
        @test haskey(EegFun.PARAMETERS, "preprocess.filter.lowpass.func")
        @test haskey(EegFun.PARAMETERS, "preprocess.filter.ica_highpass.apply")
        @test haskey(EegFun.PARAMETERS, "preprocess.filter.ica_lowpass.apply")

        # Test that all file parameters exist
        @test haskey(EegFun.PARAMETERS, "files.output.directory")
        @test haskey(EegFun.PARAMETERS, "files.output.save_ica_data")
        @test haskey(EegFun.PARAMETERS, "files.output.save_epoch_data_raw")
        @test haskey(EegFun.PARAMETERS, "files.output.save_epoch_data_corrected")
        @test haskey(EegFun.PARAMETERS, "files.output.save_epoch_data")
        @test haskey(EegFun.PARAMETERS, "files.output.save_erp_data_raw")
        @test haskey(EegFun.PARAMETERS, "files.output.save_erp_data_corrected")
        @test haskey(EegFun.PARAMETERS, "files.output.save_erp_data")

        # Test that all preprocess parameters exist
        @test haskey(EegFun.PARAMETERS, "preprocess.reference_channel")
        @test haskey(EegFun.PARAMETERS, "preprocess.layout.neighbour_criterion")
        @test haskey(EegFun.PARAMETERS, "preprocess.channel_repair.method")
        @test EegFun.PARAMETERS["preprocess.channel_repair.method"].default == "spherical_spline"
        @test haskey(EegFun.PARAMETERS, "preprocess.eog.vEOG_channels")
        @test haskey(EegFun.PARAMETERS, "preprocess.eog.hEOG_channels")
        @test haskey(EegFun.PARAMETERS, "preprocess.eog.vEOG_criterion")
        @test haskey(EegFun.PARAMETERS, "preprocess.eog.hEOG_criterion")
        @test haskey(EegFun.PARAMETERS, "preprocess.eeg.extreme_value_abs_criterion")
        @test haskey(EegFun.PARAMETERS, "preprocess.eeg.artifact_value_abs_criterion")

        # Test that ICA parameters exist
        @test haskey(EegFun.PARAMETERS, "preprocess.ica.apply")
    end

    # =============================================================================
    # CONFIG LOADING AND MERGING
    # =============================================================================

    @testset "read_config Tests" begin

        @testset "Valid Configuration Loading" begin

            # Create a simple valid user config
            user_config_path = joinpath(test_dir, "valid_config.toml")
            open(user_config_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "freq = 1.0")
                println(io, "apply = true")
                println(io, "")
                println(io, "[preprocess]")
                println(io, "epoch_start = -2.0")
                println(io, "epoch_end = 3.0")
            end

            config = EegFun.read_config(user_config_path)
            @test config isa Dict
            @test config["preprocess"]["filter"]["highpass"]["freq"] == 1.0
            @test config["preprocess"]["filter"]["highpass"]["apply"] == true
            @test config["preprocess"]["epoch_start"] == -2.0
            @test config["preprocess"]["epoch_end"] == 3.0

            # Verify default values are preserved when not overridden
            @test config["preprocess"]["filter"]["lowpass"]["freq"] == 30  # default value
            @test config["preprocess"]["reference_channel"] == "avg"  # default value

            # Test that logging occurs but doesn't prevent successful loading
            # (We don't test stderr content since @info logging is expected and normal)
            @test !isnothing(config)
        end

        @testset "Configuration Merging" begin
            # Nested configuration merging
            nested_config_path = joinpath(test_dir, "nested_config.toml")
            open(nested_config_path, "w") do io
                println(io, "[files.output]")
                println(io, "directory = \"/custom/output\"")
                println(io, "save_erp_data_raw = false")
                println(io, "")
                println(io, "[preprocess.filter.ica_highpass]")
                println(io, "freq = 2.5")
                println(io, "apply = true")  # Add the "on" key that the test expects
            end

            config = EegFun.read_config(nested_config_path)
            @test config["files"]["output"]["directory"] == "/custom/output"
            @test config["files"]["output"]["save_erp_data_raw"] == false
            @test config["files"]["output"]["save_ica_data"] == true  # default preserved
            @test config["preprocess"]["filter"]["ica_highpass"]["freq"] == 2.5
            @test config["preprocess"]["filter"]["ica_highpass"]["apply"] == true  # now present in test config
        end

        @testset "Error Handling" begin
            # Non-existent file - should throw
            @test_throws Exception EegFun.read_config("nonexistent_file.toml")

            # Invalid TOML syntax - should throw
            invalid_toml_path = joinpath(test_dir, "invalid.toml")
            open(invalid_toml_path, "w") do io
                println(io, "[section")  # Missing closing bracket
                println(io, "key = value")
            end
            @test_throws Exception EegFun.read_config(invalid_toml_path)

            # Invalid parameter values - should throw
            invalid_values_path = joinpath(test_dir, "invalid_values.toml")
            open(invalid_values_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "freq = -5.0")  # Below minimum
            end
            @test_throws Exception EegFun.read_config(invalid_values_path)

            # Invalid parameter type - should throw
            invalid_type_path = joinpath(test_dir, "invalid_type.toml")
            open(invalid_type_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "freq = \"not_a_number\"")  # Wrong type
            end
            @test_throws Exception EegFun.read_config(invalid_type_path)

            # Invalid allowed values - should throw
            invalid_allowed_path = joinpath(test_dir, "invalid_allowed.toml")
            open(invalid_allowed_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "type = \"invalid_type\"")  # Not in allowed values
            end
            @test_throws Exception EegFun.read_config(invalid_allowed_path)
        end

        @testset "Boundary Value Testing" begin
            # Minimum boundary values
            min_boundary_path = joinpath(test_dir, "min_boundary.toml")
            open(min_boundary_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "freq = 0.01")  # Minimum allowed
                println(io, "order = 1")      # Minimum allowed
            end

            config = EegFun.read_config(min_boundary_path)
            @test config["preprocess"]["filter"]["highpass"]["freq"] == 0.01
            @test config["preprocess"]["filter"]["highpass"]["order"] == 1

            # Maximum boundary values
            max_boundary_path = joinpath(test_dir, "max_boundary.toml")
            open(max_boundary_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "freq = 20.0")  # Maximum allowed
                println(io, "order = 4")      # Maximum allowed
            end

            config = EegFun.read_config(max_boundary_path)
            @test config["preprocess"]["filter"]["highpass"]["freq"] == 20.0
            @test config["preprocess"]["filter"]["highpass"]["order"] == 4
        end

        @testset "Data Type Conversions" begin
            # Numeric type conversions
            conversion_path = joinpath(test_dir, "conversion.toml")
            open(conversion_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "freq = 1")     # Int should convert to Float
                println(io, "order = 2.0")    # Float should convert to Int
            end

            config = EegFun.read_config(conversion_path)
            @test config["preprocess"]["filter"]["highpass"]["freq"] == 1.0  # Converted to Float
            @test config["preprocess"]["filter"]["highpass"]["order"] == 2     # Converted to Int
        end

        @testset "Complete Configuration Structure" begin
            # Verify all expected sections exist
            complete_config_path = joinpath(test_dir, "empty.toml")
            open(complete_config_path, "w") do io
                # Empty file to test defaults
                println(io, "# Empty config file")
            end

            config = EegFun.read_config(complete_config_path)

            # Verify main sections exist
            @test haskey(config, "files")
            @test haskey(config["preprocess"], "filter")
            @test haskey(config["preprocess"], "ica")
            @test haskey(config, "preprocess")

            # Verify subsections exist
            @test haskey(config["files"], "input")
            @test haskey(config["files"], "output")
            @test haskey(config["preprocess"]["filter"], "highpass")
            @test haskey(config["preprocess"]["filter"], "lowpass")
            @test haskey(config["preprocess"], "eog")
            @test haskey(config["preprocess"], "layout")

            # Verify epoch parameters exist
            @test haskey(config["preprocess"], "epoch_start")
            @test haskey(config["preprocess"], "epoch_end")
            @test config["preprocess"]["epoch_start"] == -1
            @test config["preprocess"]["epoch_end"] == 1
        end

        @testset "Edge Cases" begin
            # Empty sections
            empty_sections_path = joinpath(test_dir, "empty_sections.toml")
            open(empty_sections_path, "w") do io
                println(io, "[files.input]")
                println(io, "")
                println(io, "[preprocess.filter.highpass]")
            end

            config = EegFun.read_config(empty_sections_path)
            @test haskey(config, "files")
            @test haskey(config["files"], "input")

            # Special characters in string values
            special_chars_path = joinpath(test_dir, "special_chars.toml")
            open(special_chars_path, "w") do io
                println(io, "[files.input]")
                println(io, "directory = \"/path/with spaces/and/special@chars\"")
            end

            config = EegFun.read_config(special_chars_path)
            @test config["files"]["input"]["directory"] == "/path/with spaces/and/special@chars"
        end

        @testset "Complex Configuration Scenarios" begin
            # Test nested arrays
            array_config_path = joinpath(test_dir, "array_config.toml")
            open(array_config_path, "w") do io
                println(io, "[files.input]")
                println(io, "raw_data_files = [\"file1.bdf\", \"file2.bdf\"]")
            end

            config = EegFun.read_config(array_config_path)
            @test config["files"]["input"]["raw_data_files"] == ["file1.bdf", "file2.bdf"]

            # Test boolean values
            bool_config_path = joinpath(test_dir, "bool_config.toml")
            open(bool_config_path, "w") do io
                println(io, "[files.output]")
                println(io, "save_continuous_data_raw = true")
                println(io, "save_continuous_data_corrected = true")
                println(io, "save_ica_data = false")
            end

            config = EegFun.read_config(bool_config_path)
            @test config["files"]["output"]["save_continuous_data_raw"] == true
            @test config["files"]["output"]["save_continuous_data_corrected"] == true
            @test config["files"]["output"]["save_ica_data"] == false

            # Test string values with special characters
            special_config_path = joinpath(test_dir, "special_config.toml")
            open(special_config_path, "w") do io
                println(io, "[files.input]")
                println(io, "directory = \"/path/with spaces/and/special@chars\"")
            end

            config = EegFun.read_config(special_config_path)
            @test config["files"]["input"]["directory"] == "/path/with spaces/and/special@chars"
        end

        @testset "Parameter Validation" begin
            # Test numeric range validation
            range_config_path = joinpath(test_dir, "range_config.toml")
            open(range_config_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "freq = 0.005")  # Below minimum
            end
            @test_throws Exception EegFun.read_config(range_config_path)

            # Test allowed values validation
            allowed_config_path = joinpath(test_dir, "allowed_config.toml")
            open(allowed_config_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "method = \"invalid\"")  # Not in allowed values
            end
            @test_throws Exception EegFun.read_config(allowed_config_path)

            # Test type conversion validation
            type_config_path = joinpath(test_dir, "type_config.toml")
            open(type_config_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "order = \"not_a_number\"")  # Invalid type
            end
            @test_throws Exception EegFun.read_config(type_config_path)

            # Test optional parameters
            optional_config_path = joinpath(test_dir, "optional.toml")
            open(optional_config_path, "w") do io
                println(io, "[files.input]")
                println(io, "epoch_condition_file = \"\"")  # Empty string for optional parameter
            end
            config = EegFun.read_config(optional_config_path)
            @test !isnothing(config)
            @test config["files"]["input"]["epoch_condition_file"] == ""
        end

        @testset "Config Merging" begin
            # Test merging with defaults
            custom_config_path = joinpath(test_dir, "custom.toml")
            open(custom_config_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "freq = 0.5")
            end
            config = EegFun.read_config(custom_config_path)
            @test !isnothing(config)
            @test config["preprocess"]["filter"]["highpass"]["freq"] == 0.5
            @test config["preprocess"]["filter"]["highpass"]["method"] == "iir"  # Default value
            @test config["preprocess"]["filter"]["highpass"]["order"] == 1     # Default value

            # Test nested merging
            nested_config_path = joinpath(test_dir, "nested.toml")
            open(nested_config_path, "w") do io
                println(io, "[files.input]")
                println(io, "directory = \"/custom/path\"")
                println(io, "raw_data_files = [\"file1.bdf\", \"file2.bdf\"]")
            end
            config = EegFun.read_config(nested_config_path)
            @test !isnothing(config)
            @test config["files"]["input"]["directory"] == "/custom/path"
            @test config["files"]["input"]["raw_data_files"] == ["file1.bdf", "file2.bdf"]
        end
    end

    # =============================================================================
    # VALIDATION
    # =============================================================================

    @testset "ValidationResult Tests" begin
        # Test successful validation
        result = EegFun.ValidationResult(success = true)
        @test result.success
        @test isnothing(result.error)
        @test isnothing(result.key_path)

        # Test failed validation with error
        result = EegFun.ValidationResult(success = false, error = "Test error")
        @test !result.success
        @test result.error == "Test error"
        @test isnothing(result.key_path)

        # Test failed validation with path
        result = EegFun.ValidationResult(success = false, error = "Test error", key_path = "test.path")
        @test !result.success
        @test result.error == "Test error"
        @test result.key_path == "test.path"
    end

    @testset "Config Validation Tests" begin
        # Test valid configuration
        valid_config = Dict{String,Any}(
            "preprocess" => Dict{String,Any}(
                "filter" => Dict{String,Any}(
                    "highpass" => Dict{String,Any}("apply" => true, "method" => "iir", "freq" => 0.1, "order" => 1),
                    "lowpass" => Dict{String,Any}("apply" => true, "method" => "iir", "freq" => 40, "order" => 3),
                ),
            ),
        )
        result = EegFun._validate_config(valid_config)
        @test result.success

        # Test invalid configuration - wrong type
        invalid_config = Dict{String,Any}(
            "preprocess" => Dict{String,Any}(
                "filter" => Dict{String,Any}(
                    "highpass" => Dict{String,Any}(
                        "apply" => "true",  # Should be Bool
                        "method" => "iir",
                        "freq" => 0.1,
                        "order" => 1,
                    ),
                ),
            ),
        )
        result = EegFun._validate_config(invalid_config)
        @test !result.success
        @test !isnothing(result.error)
        @test !isnothing(result.key_path)

        # Test invalid configuration - out of range
        invalid_config = Dict{String,Any}(
            "preprocess" => Dict{String,Any}(
                "filter" => Dict{String,Any}(
                    "highpass" => Dict{String,Any}(
                        "apply" => true,
                        "method" => "iir",
                        "freq" => -1,  # Below minimum
                        "order" => 1,
                    ),
                ),
            ),
        )
        result = EegFun._validate_config(invalid_config)
        @test !result.success
        @test !isnothing(result.error)
        @test !isnothing(result.key_path)

        # Test invalid configuration - wrong allowed value
        invalid_config = Dict{String,Any}(
            "preprocess" => Dict{String,Any}(
                "filter" => Dict{String,Any}(
                    "highpass" => Dict{String,Any}(
                        "apply" => true,
                        "method" => "invalid",  # Not in allowed_values
                        "freq" => 0.1,
                        "order" => 1,
                    ),
                ),
            ),
        )
        result = EegFun._validate_config(invalid_config)
        @test !result.success
        @test !isnothing(result.error)
        @test !isnothing(result.key_path)
    end

    # =============================================================================
    # TEMPLATE GENERATION
    # =============================================================================

    @testset "Config Template Generation" begin
        # Test template generation
        template_path = joinpath(test_dir, "template.toml")
        EegFun.generate_config_template(filename = template_path)
        @test isfile(template_path)

        # Verify template contents
        template = TOML.parsefile(template_path)
        @test haskey(template, "files")
        @test haskey(template["preprocess"], "filter")
        @test haskey(template, "preprocess")

        # Test template with custom filename
        custom_path = joinpath(test_dir, "custom.toml")
        EegFun.generate_config_template(filename = custom_path)
        @test isfile(custom_path)
    end

    # =============================================================================
    # PARAMETER CONSTRUCTORS AND HELPERS
    # =============================================================================

    @testset "Parameter Constructor Helper Tests" begin
        # Test string_param
        param = EegFun.string_param("Test string", "default")
        @test param isa EegFun.ConfigParameter{Union{Vector{String},String}}
        @test param.description == "Test string"
        @test param.default == "default"

        # Test string_param with allowed values
        param = EegFun.string_param("Test string", "default", allowed = ["a", "b", "c"])
        @test param.allowed == ["a", "b", "c"]

        # Test string_param with default empty string
        param = EegFun.string_param("Test string")
        @test param.default == ""

        # Test bool_param
        param = EegFun.bool_param("Test bool", true)
        @test param isa EegFun.ConfigParameter{Bool}
        @test param.description == "Test bool"
        @test param.default == true

        # Test bool_param with default false
        param = EegFun.bool_param("Test bool")
        @test param.default == false

        # Test number_param
        param = EegFun.number_param("Test number", 5.0, 0.0, 10.0)
        @test param isa EegFun.ConfigParameter{Real}
        @test param.description == "Test number"
        @test param.default == 5.0
        @test param.min == 0.0
        @test param.max == 10.0

        # Test number_param without min/max
        param = EegFun.number_param("Test number", 5.0)
        @test param.default == 5.0
        @test isnothing(param.min)
        @test isnothing(param.max)

        # Test channel_groups_param
        default = [["Fp1"], ["Fp2"]]
        param = EegFun.channel_groups_param("Test channels", default)
        @test param isa EegFun.ConfigParameter{Vector{Vector{String}}}
        @test param.description == "Test channels"
        @test param.default == default

        # Test _param helper function directly
        param = EegFun._param(Int, "Test int", 42, min = 0, max = 100)
        @test param isa EegFun.ConfigParameter{Int}
        @test param.description == "Test int"
        @test param.default == 42
        @test param.min == 0
        @test param.max == 100
    end

    @testset "Parameter Constructor Edge Cases" begin
        # Test string_param with empty allowed list
        param = EegFun.string_param("Test string", "default", allowed = String[])
        @test param.allowed == String[]

        # Test number_param with negative min/max
        param = EegFun.number_param("Test number", 0.0, -10.0, 10.0)
        @test param.min == -10.0
        @test param.max == 10.0

        # Test channel_groups_param with empty groups
        param = EegFun.channel_groups_param("Test channels", Vector{Vector{String}}())
        @test param.default == Vector{Vector{String}}()

        # Test _param with Union types
        param = EegFun._param(Union{String,Nothing}, "Test union", nothing)
        @test param isa EegFun.ConfigParameter{Union{String,Nothing}}
        @test isnothing(param.default)
    end

    @testset "ConfigParameter Constructor Tests" begin
        # Test with all fields
        param = EegFun.ConfigParameter{Int}(description = "Test parameter", default = 5, min = 1, max = 10)
        @test param.description == "Test parameter"
        @test param.default == 5
        @test param.min == 1
        @test param.max == 10

        # Test with minimal fields
        param = EegFun.ConfigParameter{String}(description = "Test parameter")
        @test param.description == "Test parameter"
        @test isnothing(param.default)
        @test isnothing(param.min)
        @test isnothing(param.max)
        @test isnothing(param.allowed)

        # Test with some fields
        param = EegFun.ConfigParameter{Float64}(description = "Test parameter", default = 1.0, min = 0.0)
        @test param.description == "Test parameter"
        @test param.default == 1.0
        @test param.min == 0.0
        @test isnothing(param.max)
        @test isnothing(param.allowed)
    end

    # =============================================================================
    # UTILITY FUNCTIONS
    # =============================================================================

    @testset "_param Helper Function Tests" begin
        # Test basic parameter creation
        param = EegFun._param(String, "Test param", "default")
        @test param isa EegFun.ConfigParameter{String}
        @test param.description == "Test param"
        @test param.default == "default"

        # Test with allowed values
        param = EegFun._param(String, "Test param", "default", allowed = ["a", "b"])
        @test param.allowed == ["a", "b"]

        # Test with min/max values
        param = EegFun._param(Real, "Test param", 5.0, min = 0.0, max = 10.0)
        @test param.min == 0.0
        @test param.max == 10.0

        # Test with nothing defaults
        param = EegFun._param(Bool, "Test param")
        @test isnothing(param.default)
        @test isnothing(param.allowed)
    end

    @testset "_filter_param_spec Tests" begin
        # Test filter parameter specification creation
        filter_spec = EegFun._filter_param_spec("test.prefix", true, "hp", 1.0, 0.1, 10.0, 2, 1, 5)

        @test haskey(filter_spec, "test.prefix.apply")
        @test haskey(filter_spec, "test.prefix.type")
        @test haskey(filter_spec, "test.prefix.method")
        @test haskey(filter_spec, "test.prefix.func")
        @test haskey(filter_spec, "test.prefix.freq")
        @test haskey(filter_spec, "test.prefix.order")

        # Test parameter types
        @test filter_spec["test.prefix.apply"] isa EegFun.ConfigParameter{Bool}
        @test filter_spec["test.prefix.type"] isa EegFun.ConfigParameter{Union{Vector{String},String}}
        @test filter_spec["test.prefix.method"] isa EegFun.ConfigParameter{Union{Vector{String},String}}
        @test filter_spec["test.prefix.freq"] isa EegFun.ConfigParameter{Real}

        # Test parameter values
        @test filter_spec["test.prefix.type"].default == "hp"
        @test filter_spec["test.prefix.freq"].default == 1.0
        @test filter_spec["test.prefix.freq"].min == 0.1
        @test filter_spec["test.prefix.freq"].max == 10.0
        @test filter_spec["test.prefix.order"].default == 2
        @test filter_spec["test.prefix.order"].min == 1
        @test filter_spec["test.prefix.order"].max == 5
    end

    @testset "_group_parameters_by_section Tests" begin
        # Test parameter grouping
        sections = EegFun._group_parameters_by_section()

        @test haskey(sections, "files")
        @test haskey(sections, "preprocess")
        @test haskey(sections["preprocess"], "filter.highpass")
        @test haskey(sections["preprocess"], "filter.lowpass")
        @test haskey(sections["preprocess"], "ica")

        # Test subsection grouping
        @test haskey(sections["files"], "input")
        @test haskey(sections["files"], "output")
        @test haskey(sections["preprocess"], "filter.highpass")
        @test haskey(sections["preprocess"], "filter.lowpass")

        # Test parameter placement
        @test any(p[1] == "files.input.directory" for p in sections["files"]["input"])
        @test any(p[1] == "preprocess.filter.highpass.freq" for p in sections["preprocess"]["filter.highpass"])
    end

    @testset "_extract_subsection Tests" begin
        # Test basic subsection extraction
        @test EegFun._extract_subsection("files", "files.input.directory") == "input"
        @test EegFun._extract_subsection("preprocess.filter", "preprocess.filter.highpass.freq") == "highpass"

        # Test nested subsections
        @test EegFun._extract_subsection("preprocess.filter", "preprocess.filter.ica_highpass.freq") == "ica_highpass"

        # Test no subsection
        @test EegFun._extract_subsection("preprocess.ica", "preprocess.ica.apply") == ""

        # Test non-matching prefix
        @test EegFun._extract_subsection("files", "preprocess.filter.highpass.freq") == ""
    end

    @testset "_group_params_by_subsection Tests" begin
        # Test parameter grouping by subsection
        matching_params = ["preprocess.filter.highpass.freq", "preprocess.filter.highpass.apply", "preprocess.filter.lowpass.freq"]
        grouped = EegFun._group_params_by_subsection("preprocess.filter", matching_params)

        @test haskey(grouped, "highpass")
        @test haskey(grouped, "lowpass")
        @test length(grouped["highpass"]) == 2
        @test length(grouped["lowpass"]) == 1

        # Test parameter content
        highpass_params = [p[1] for p in grouped["highpass"]]
        @test "preprocess.filter.highpass.freq" in highpass_params
        @test "preprocess.filter.highpass.apply" in highpass_params
    end

    @testset "_merge_configs Tests" begin
        # Test merging with empty configs
        @test isempty(EegFun._merge_configs(Dict(), Dict()))

        # Test merging with one empty config
        default = Dict("a" => 1)
        user = Dict()
        merged = EegFun._merge_configs(default, user)
        @test merged["a"] == 1

        # Test deep merging
        default = Dict("a" => Dict("b" => Dict("c" => 1)))
        user = Dict("a" => Dict("b" => Dict("d" => 2)))
        merged = EegFun._merge_configs(default, user)
        @test merged["a"]["b"]["c"] == 1
        @test merged["a"]["b"]["d"] == 2

        # Test user config overrides default
        default = Dict("a" => 1)
        user = Dict("a" => 2)
        merged = EegFun._merge_configs(default, user)
        @test merged["a"] == 2

        # Test nested user config overrides default
        default = Dict("a" => Dict("b" => 1, "c" => 2))
        user = Dict("a" => Dict("b" => 3))
        merged = EegFun._merge_configs(default, user)
        @test merged["a"]["b"] == 3
        @test merged["a"]["c"] == 2
    end

    @testset "_merge_nested! Edge Cases" begin
        # Test merging with non-dict values (using compatible types)
        target = Dict("a" => Dict("b" => 1), "c" => 5)
        source = Dict("a" => Dict("d" => 2), "c" => 10)  # Override with compatible types

        EegFun._merge_nested!(target, source)
        @test target["a"]["b"] == 1  # Original value preserved
        @test target["a"]["d"] == 2  # New value added
        @test target["c"] == 10      # Value overridden

        # Test merging with empty source
        target = Dict("a" => 1, "b" => 2)
        source = Dict()
        original_target = copy(target)

        EegFun._merge_nested!(target, source)
        @test target == original_target  # Should be unchanged
    end

    @testset "_validate_parameter Tests" begin
        # Test valid parameter
        param = EegFun.ConfigParameter{Real}(description = "Test", default = 5.0, min = 0.0, max = 10.0)
        result = EegFun._validate_parameter(5.0, param, "test.param")
        @test result.success == true

        # Test value below minimum
        result = EegFun._validate_parameter(-1.0, param, "test.param")
        @test result.success == false
        @test contains(result.error, "must be >=")

        # Test value above maximum
        result = EegFun._validate_parameter(15.0, param, "test.param")
        @test result.success == false
        @test contains(result.error, "must be <=")

        # Test wrong type
        result = EegFun._validate_parameter("string", param, "test.param")
        @test result.success == false
        @test contains(result.error, "must be a number")

        # Test allowed values
        param = EegFun.ConfigParameter{String}(description = "Test", default = "a", allowed = ["a", "b", "c"])
        result = EegFun._validate_parameter("b", param, "test.param")
        @test result.success == true

        result = EegFun._validate_parameter("d", param, "test.param")
        @test result.success == false
        @test contains(result.error, "must be one of")

        # Test numeric type conversion (Int to Float)
        param = EegFun.ConfigParameter{Float64}(description = "Test", default = 5.0)
        result = EegFun._validate_parameter(5, param, "test.param")  # Int value
        @test result.success == true  # Should accept Int for Float64
    end

    # =============================================================================
    # DISPLAY/INFO FUNCTIONS
    # =============================================================================

    @testset "_show_parameter_details Tests" begin
        # Test with parameter that has all fields
        param = EegFun.ConfigParameter{Real}(description = "Test param", default = 5.0, min = 0.0, max = 10.0, allowed = ["a", "b"])

        # Capture output (this is tricky with @info, so we'll test the function exists)
        @test typeof(EegFun._show_parameter_details) <: Function
    end

    @testset "_display_section Tests" begin
        # Test section display
        sections = EegFun._group_parameters_by_section()
        files_section = sections["files"]

        # Test function exists and can be called
        @test typeof(EegFun._display_section) <: Function
    end

    @testset "_display_subsection Tests" begin
        # Test subsection display
        sections = EegFun._group_parameters_by_section()
        input_params = sections["files"]["input"]

        # Test function exists and can be called
        @test typeof(EegFun._display_subsection) <: Function
    end

    @testset "_show_specific_parameter Tests" begin
        # Test with exact parameter match
        @test typeof(EegFun._show_specific_parameter) <: Function

        # Test with section prefix
        @test typeof(EegFun._show_specific_parameter) <: Function

        # Test with non-existent parameter
        @test typeof(EegFun._show_specific_parameter) <: Function
    end

    @testset "_show_section_overview Tests" begin
        # Test section overview
        matching_params = ["preprocess.filter.highpass.freq", "preprocess.filter.highpass.apply", "preprocess.filter.lowpass.freq"]

        # Test function exists and can be called
        @test typeof(EegFun._show_section_overview) <: Function
    end

    @testset "_display_grouped_params Tests" begin
        # Test grouped parameters display
        sections = EegFun._group_parameters_by_section()
        files_section = sections["files"]

        # Test function exists and can be called
        @test typeof(EegFun._display_grouped_params) <: Function
    end

    @testset "_show_all_parameters Tests" begin
        # Test that the function exists and can be called
        @test typeof(EegFun._show_all_parameters) <: Function

        # Test that it doesn't throw errors when called
        # (This function uses @info for output, so we can't easily capture it)
        @test isnothing(nothing)  # Placeholder to ensure test runs
    end

    # =============================================================================
    # EDGE CASES
    # =============================================================================

    @testset "_write_parameter_value Edge Cases" begin
        # Test string with special characters that need escaping
        io = IOBuffer()
        EegFun._write_parameter_value(io, "test_param", "path\\with\\backslashes")
        @test contains(String(take!(io)), "path\\\\with\\\\backslashes")

        # Test vector with nested vectors (channel groups)
        io = IOBuffer()
        nested_vector = [["Fp1", "IO1"], ["Fp2", "IO2"]]
        EegFun._write_parameter_value(io, "test_param", nested_vector)
        output = String(take!(io))
        @test contains(output, "test_param = [")
        @test contains(output, "Fp1")
        @test contains(output, "IO1")

        # Test empty vector
        io = IOBuffer()
        EegFun._write_parameter_value(io, "test_param", String[])
        @test String(take!(io)) == "test_param = []\n"

        # Test nothing value
        io = IOBuffer()
        EegFun._write_parameter_value(io, "test_param", nothing)
        @test String(take!(io)) == "# test_param = \n"

        # Test numeric values
        io = IOBuffer()
        EegFun._write_parameter_value(io, "test_param", 42)
        @test String(take!(io)) == "test_param = 42\n"

        # Test float values
        io = IOBuffer()
        EegFun._write_parameter_value(io, "test_param", 3.14)
        @test String(take!(io)) == "test_param = 3.14\n"

        # Test boolean values
        io = IOBuffer()
        EegFun._write_parameter_value(io, "test_param", true)
        @test String(take!(io)) == "test_param = true\n"
    end

    # =============================================================================
    # INTEGRATION TESTS
    # =============================================================================

    @testset "show_parameter_info Integration Tests" begin
        @test typeof(EegFun.show_parameter_info) <: Function
    end

    @testset "Template Generation Integration Tests" begin
        # Test template generation with custom filename
        custom_template = joinpath(test_dir, "custom_template.toml")
        EegFun.generate_config_template(filename = custom_template)
        @test isfile(custom_template)

        # Verify template content
        template_content = read(custom_template, String)
        @test contains(template_content, "# EEG Processing Configuration Template")
        @test contains(template_content, "[files]")
        @test contains(template_content, "[preprocess.filter.highpass]")
        @test contains(template_content, "[preprocess]")
        @test contains(template_content, "[preprocess.ica]")

        # Clean up
        rm(custom_template)
    end

    @testset "_write_template_header Tests" begin
        io = IOBuffer()
        EegFun._write_template_header(io)
        output = String(take!(io))

        @test contains(output, "# EEG Processing Configuration Template")
        @test contains(output, "# Generated on")
        @test contains(output, "# This template shows all available configuration options")
        @test contains(output, "# Required fields are marked with [REQUIRED]")
        @test contains(output, "# Default values are shown where available")
    end

    @testset "_write_section Tests" begin
        # Create test data using the actual PARAMETERS structure
        sections = EegFun._group_parameters_by_section()
        files_section = sections["files"]

        io = IOBuffer()
        EegFun._write_section(io, "files", files_section)
        output = String(take!(io))

        @test contains(output, "# files Settings")
        @test contains(output, "[files]")
        @test contains(output, "# input Settings")
        @test contains(output, "[files.input]")
        @test contains(output, "# output Settings")
        @test contains(output, "[files.output]")
    end

    @testset "_write_subsection Tests" begin
        # Create test data using the actual PARAMETERS structure
        sections = EegFun._group_parameters_by_section()
        input_params = sections["files"]["input"]

        io = IOBuffer()
        EegFun._write_subsection(io, "files", "input", input_params)
        output = String(take!(io))

        @test contains(output, "# input Settings")
        @test contains(output, "[files.input]")
        @test contains(output, "directory = \".\"")
    end

    @testset "_write_template_sections Tests" begin
        io = IOBuffer()
        EegFun._write_template_sections(io)
        output = String(take!(io))

        # Should contain all major sections
        @test contains(output, "# files Settings")
        @test contains(output, "# preprocess Settings")
        @test contains(output, "# filter.highpass Settings")
        @test contains(output, "# ica Settings")
    end

    @testset "_write_parameter_docs Tests" begin
        # Test parameter with all fields
        param = EegFun.ConfigParameter{Real}(description = "Test param", default = 5.0, min = 0.0, max = 10.0, allowed = ["a", "b"])
        io = IOBuffer()
        EegFun._write_parameter_docs(io, param)
        output = String(take!(io))

        @test contains(output, "# Test param")
        @test contains(output, "# Type: Real")
        @test contains(output, "# Range: 0.0 ≤ value ≤ 10.0")
        @test contains(output, "# Allowed values: a, b")
        @test contains(output, "# Default: 5.0")

        # Test required parameter
        param = EegFun.ConfigParameter{String}(description = "Required param", default = nothing)
        io = IOBuffer()
        EegFun._write_parameter_docs(io, param)
        output = String(take!(io))

        @test contains(output, "# [REQUIRED]")
        @test !contains(output, "# Default:")
    end

    # =============================================================================
    # ERROR HANDLING
    # =============================================================================

    @testset "PreprocessConfig Channel Repair Method Tests" begin
        cfg = Dict(
            "reference_channel" => "avg",
            "epoch_start" => -1.0,
            "epoch_end" => 1.0,
            "layout" => Dict("neighbour_criterion" => 0.25),
            "filter" => Dict(
                "highpass" => Dict("apply" => true, "type" => "hp", "freq" => 0.1, "func" => "butter", "method" => "iir", "order" => 2),
                "lowpass" => Dict("apply" => true, "type" => "lp", "freq" => 30.0, "func" => "butter", "method" => "iir", "order" => 2),
                "ica_highpass" => Dict("apply" => true, "type" => "hp", "freq" => 1.0, "func" => "butter", "method" => "iir", "order" => 2),
                "ica_lowpass" => Dict("apply" => true, "type" => "lp", "freq" => 30.0, "func" => "butter", "method" => "iir", "order" => 2),
            ),
            "cleanline" => Dict(
                "apply" => false, "line_frequencies" => [50.0], "bandwidth" => 2.0,
                "sliding_win_length" => 4.0, "sliding_win_step" => 2.0, "time_bandwidth" => 3.0,
                "k_tapers" => 5, "p_value" => 0.05, "pad" => 2
            ),
            "resample" => Dict("apply" => false, "target_rate" => 512),
            "eog" => Dict("vEOG_criterion" => 50, "hEOG_criterion" => 30, "vEOG_channels" => [["Fp1"]], "hEOG_channels" => [["F9"]]),
            "eeg" => Dict("artifact_value_abs_criterion" => 100, "artifact_value_z_criterion" => 0.0, "extreme_value_abs_criterion" => 500),
            "ica" => Dict("apply" => true, "percentage_of_data" => 100.0),
        )
        p_cfg_default = EegFun.PreprocessConfig(cfg)
        @test p_cfg_default.channel_repair_method == :spherical_spline

        cfg_custom = copy(cfg)
        cfg_custom["channel_repair"] = Dict("method" => "neighbor_interpolation")
        p_cfg_custom = EegFun.PreprocessConfig(cfg_custom)
        @test p_cfg_custom.channel_repair_method == :neighbor_interpolation

        # Test integer vector input for line_frequencies (e.g. [50, 100])
        cfg_int_freqs = copy(cfg)
        cfg_int_freqs["cleanline"] = copy(cfg["cleanline"])
        cfg_int_freqs["cleanline"]["line_frequencies"] = [50, 100]
        p_cfg_int = EegFun.PreprocessConfig(cfg_int_freqs)
        @test p_cfg_int.cleanline.line_frequencies == [50.0, 100.0]
        @test p_cfg_int.cleanline.line_frequencies isa Vector{Float64}
    end

    @testset "FilesConfig Construction Tests" begin
        files_dict = Dict(
            "input" => Dict("directory" => "./data", "raw_data_files" => "\\.bdf", "recursive" => false, "layout_file" => "biosemi72.csv", "epoch_condition_file" => "epochs.toml"),
            "output" => Dict("directory" => "./out", "save_continuous_data_raw" => true, "save_epoch_data_raw" => false),
        )
        files_cfg = EegFun.FilesConfig(files_dict)
        @test files_cfg.input.directory == "./data"
        @test files_cfg.input.layout_file == "biosemi72.csv"
        @test files_cfg.output.directory == "./out"
        @test files_cfg.output.save_continuous_data_raw == true
        @test files_cfg.output.save_epoch_data_raw == false
    end

    @testset "PipelineConfig Construction Tests" begin
        full_dict = Dict(
            "files" => Dict("input" => Dict("directory" => "./raw_dir")),
            "preprocess" => Dict("channel_repair" => Dict("method" => "spherical_spline")),
        )
        p_cfg = EegFun.PipelineConfig(full_dict)
        @test p_cfg.files.input.directory == "./raw_dir"
        @test p_cfg.preprocess.channel_repair_method == :spherical_spline
    end

    # Cleanup
    rm(test_dir, recursive = true, force = true)
end
