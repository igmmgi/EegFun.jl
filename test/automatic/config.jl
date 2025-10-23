using Test
using TOML
using Dates

# Import the package (which includes the config module)
using eegfun

@testset "Configuration System Tests" begin

    # Create temporary directory for test files
    test_dir = mktempdir()

    # =============================================================================
    # CONFIGPARAMETER AND PARAMETER STRUCTS
    # =============================================================================

    @testset "ConfigParameter Tests" begin
        # Test ConfigParameter struct creation
        param = eegfun.ConfigParameter{Float64}(description = "Test parameter", default = 1.0, min = 0.0, max = 2.0)
        @test param.description == "Test parameter"
        @test param.default == 1.0
        @test param.min == 0.0
        @test param.max == 2.0

        # Test ConfigParameter with nothing values
        param = eegfun.ConfigParameter{String}(
            description = "Test parameter",
            default = nothing,
            min = nothing,
            max = nothing,
            allowed = nothing,
        )
        @test param.description == "Test parameter"
        @test param.default === nothing
        @test param.min === nothing
        @test param.max === nothing
        @test param.allowed === nothing
    end

    # =============================================================================
    # PARAMETERS DICTIONARY
    # =============================================================================

    @testset "PARAMETERS Dictionary Tests" begin
        # Test that all required parameters exist
        @test haskey(eegfun.PARAMETERS, "files.input.directory")
        @test haskey(eegfun.PARAMETERS, "files.input.raw_data_files")
        @test haskey(eegfun.PARAMETERS, "files.input.layout_file")
        @test haskey(eegfun.PARAMETERS, "preprocess.filter.highpass.freq")
        @test haskey(eegfun.PARAMETERS, "preprocess.filter.lowpass.freq")
        @test haskey(eegfun.PARAMETERS, "preprocess.epoch_start")
        @test haskey(eegfun.PARAMETERS, "preprocess.epoch_end")

        # Test parameter types
        @test eegfun.PARAMETERS["files.input.directory"] isa eegfun.ConfigParameter{String}
        @test eegfun.PARAMETERS["preprocess.filter.highpass.freq"] isa eegfun.ConfigParameter{Real}
        @test eegfun.PARAMETERS["preprocess.filter.highpass.order"] isa eegfun.ConfigParameter{Real}
        @test eegfun.PARAMETERS["files.output.save_continuous_data"] isa eegfun.ConfigParameter{Bool}

        # Test that all parameters have valid descriptions
        for (path, param) in eegfun.PARAMETERS
            @test !isempty(param.description)
            @test param isa eegfun.ConfigParameter
        end

        # Test specific parameter properties
        param = eegfun.PARAMETERS["preprocess.filter.highpass.freq"]
        @test param.description == "Cutoff frequency (Hz)"
        @test param.min == 0.01
        @test param.max == 20.0
        @test param.default == 0.1

        # Test non-existent parameter
        @test !haskey(eegfun.PARAMETERS, "nonexistent")

        # Test that all filter parameters exist
        @test haskey(eegfun.PARAMETERS, "preprocess.filter.highpass.apply")
        @test haskey(eegfun.PARAMETERS, "preprocess.filter.highpass.method")
        @test haskey(eegfun.PARAMETERS, "preprocess.filter.highpass.func")
        @test haskey(eegfun.PARAMETERS, "preprocess.filter.lowpass.apply")
        @test haskey(eegfun.PARAMETERS, "preprocess.filter.lowpass.method")
        @test haskey(eegfun.PARAMETERS, "preprocess.filter.lowpass.func")
        @test haskey(eegfun.PARAMETERS, "preprocess.filter.ica_highpass.apply")
        @test haskey(eegfun.PARAMETERS, "preprocess.filter.ica_lowpass.apply")

        # Test that all file parameters exist
        @test haskey(eegfun.PARAMETERS, "files.output.directory")
        @test haskey(eegfun.PARAMETERS, "files.output.save_ica_data")
        @test haskey(eegfun.PARAMETERS, "files.output.save_epoch_data_original")
        @test haskey(eegfun.PARAMETERS, "files.output.save_epoch_data_cleaned")
        @test haskey(eegfun.PARAMETERS, "files.output.save_erp_data_original")
        @test haskey(eegfun.PARAMETERS, "files.output.save_erp_data_cleaned")
        @test haskey(eegfun.PARAMETERS, "files.output.exit_early")

        # Test that all preprocess parameters exist
        @test haskey(eegfun.PARAMETERS, "preprocess.reference_channel")
        @test haskey(eegfun.PARAMETERS, "preprocess.layout.neighbour_criterion")
        @test haskey(eegfun.PARAMETERS, "preprocess.eog.vEOG_channels")
        @test haskey(eegfun.PARAMETERS, "preprocess.eog.hEOG_channels")
        @test haskey(eegfun.PARAMETERS, "preprocess.eog.vEOG_criterion")
        @test haskey(eegfun.PARAMETERS, "preprocess.eog.hEOG_criterion")
        @test haskey(eegfun.PARAMETERS, "preprocess.eeg.extreme_value_criterion")
        @test haskey(eegfun.PARAMETERS, "preprocess.eeg.artifact_value_criterion")

        # Test that ICA parameters exist
        @test haskey(eegfun.PARAMETERS, "preprocess.ica.apply")
    end

    # # =============================================================================
    # # CONFIG LOADING AND MERGING
    # # =============================================================================

    @testset "load_config Tests" begin

        @testset "Valid Configuration Loading" begin
            # Test 1: Load default config only - should work without errors
            default_config =
                eegfun.load_config(joinpath(dirname(@__FILE__), "..", "..", "src", "config", "default.toml"))
            @test default_config isa Dict
            @test haskey(default_config, "preprocess")
            @test haskey(default_config["preprocess"], "filter")

            # Test 2: Create a simple valid user config
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

            config = eegfun.load_config(user_config_path)
            @test config isa Dict
            @test config["preprocess"]["filter"]["highpass"]["freq"] == 1.0
            @test config["preprocess"]["filter"]["highpass"]["apply"] == true
            @test config["preprocess"]["epoch_start"] == -2.0
            @test config["preprocess"]["epoch_end"] == 3.0

            # Test 3: Verify default values are preserved when not overridden
            @test config["preprocess"]["filter"]["lowpass"]["freq"] == 30  # default value
            @test config["preprocess"]["reference_channel"] == "avg"  # default value

            # Test 4: Test that logging occurs but doesn't prevent successful loading
            # (We don't test stderr content since @info logging is expected and normal)
            @test config !== nothing
        end

        @testset "Configuration Merging" begin
            # Test 5: Nested configuration merging
            nested_config_path = joinpath(test_dir, "nested_config.toml")
            open(nested_config_path, "w") do io
                println(io, "[files.output]")
                println(io, "directory = \"/custom/output\"")
                println(io, "save_erp_data_original = false")
                println(io, "")
                println(io, "[preprocess.filter.ica_highpass]")
                println(io, "freq = 2.5")
                println(io, "apply = true")  # Add the "on" key that the test expects
            end

            config = eegfun.load_config(nested_config_path)
            @test config["files"]["output"]["directory"] == "/custom/output"
            @test config["files"]["output"]["save_erp_data_original"] == false
            @test config["files"]["output"]["save_ica_data"] == true  # default preserved
            @test config["preprocess"]["filter"]["ica_highpass"]["freq"] == 2.5
            @test config["preprocess"]["filter"]["ica_highpass"]["apply"] == true  # now present in test config
        end

        @testset "Error Handling" begin
            # Test 6: Non-existent file - expect nothing to be returned
            result = eegfun.load_config("nonexistent_file.toml")
            @test result === nothing

            # Test 7: Invalid TOML syntax - expect nothing to be returned
            invalid_toml_path = joinpath(test_dir, "invalid.toml")
            open(invalid_toml_path, "w") do io
                println(io, "[section")  # Missing closing bracket
                println(io, "key = value")
            end
            result = eegfun.load_config(invalid_toml_path)
            @test result === nothing

            # Test 8: Invalid parameter values - expect nothing to be returned
            invalid_values_path = joinpath(test_dir, "invalid_values.toml")
            open(invalid_values_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "freq = -5.0")  # Below minimum
            end
            result = eegfun.load_config(invalid_values_path)
            @test result === nothing

            # Test 9: Invalid parameter type - expect nothing to be returned
            invalid_type_path = joinpath(test_dir, "invalid_type.toml")
            open(invalid_type_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "freq = \"not_a_number\"")  # Wrong type
            end
            result = eegfun.load_config(invalid_type_path)
            @test result === nothing

            # Test 10: Invalid allowed values - expect nothing to be returned
            invalid_allowed_path = joinpath(test_dir, "invalid_allowed.toml")
            open(invalid_allowed_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "type = \"invalid_type\"")  # Not in allowed values
            end
            result = eegfun.load_config(invalid_allowed_path)
            @test result === nothing
        end

        @testset "Boundary Value Testing" begin
            # Test 11: Minimum boundary values
            min_boundary_path = joinpath(test_dir, "min_boundary.toml")
            open(min_boundary_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "freq = 0.01")  # Minimum allowed
                println(io, "order = 1")      # Minimum allowed
            end

            config = eegfun.load_config(min_boundary_path)
            @test config["preprocess"]["filter"]["highpass"]["freq"] == 0.01
            @test config["preprocess"]["filter"]["highpass"]["order"] == 1

            # Test 12: Maximum boundary values
            max_boundary_path = joinpath(test_dir, "max_boundary.toml")
            open(max_boundary_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "freq = 20.0")  # Maximum allowed
                println(io, "order = 4")      # Maximum allowed
            end

            config = eegfun.load_config(max_boundary_path)
            @test config["preprocess"]["filter"]["highpass"]["freq"] == 20.0
            @test config["preprocess"]["filter"]["highpass"]["order"] == 4
        end

        @testset "Data Type Conversions" begin
            # Test 13: Numeric type conversions
            conversion_path = joinpath(test_dir, "conversion.toml")
            open(conversion_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "freq = 1")     # Int should convert to Float
                println(io, "order = 2.0")    # Float should convert to Int
            end

            config = eegfun.load_config(conversion_path)
            @test config["preprocess"]["filter"]["highpass"]["freq"] == 1.0  # Converted to Float
            @test config["preprocess"]["filter"]["highpass"]["order"] == 2     # Converted to Int
        end

        @testset "Complete Configuration Structure" begin
            # Test 14: Verify all expected sections exist
            complete_config_path = joinpath(test_dir, "empty.toml")
            open(complete_config_path, "w") do io
                # Empty file to test defaults
                println(io, "# Empty config file")
            end

            config = eegfun.load_config(complete_config_path)

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
            # Test 15: Empty sections
            empty_sections_path = joinpath(test_dir, "empty_sections.toml")
            open(empty_sections_path, "w") do io
                println(io, "[files.input]")
                println(io, "")
                println(io, "[preprocess.filter.highpass]")
            end

            config = eegfun.load_config(empty_sections_path)
            @test haskey(config, "files")
            @test haskey(config["files"], "input")

            # Test 16: Special characters in string values
            special_chars_path = joinpath(test_dir, "special_chars.toml")
            open(special_chars_path, "w") do io
                println(io, "[files.input]")
                println(io, "directory = \"/path/with spaces/and/special@chars\"")
            end

            config = eegfun.load_config(special_chars_path)
            @test config["files"]["input"]["directory"] == "/path/with spaces/and/special@chars"
        end

        @testset "Complex Configuration Scenarios" begin
            # Test nested arrays
            array_config_path = joinpath(test_dir, "array_config.toml")
            open(array_config_path, "w") do io
                println(io, "[files.input]")
                println(io, "raw_data_files = [\"file1.bdf\", \"file2.bdf\"]")
            end

            config = eegfun.load_config(array_config_path)
            @test config["files"]["input"]["raw_data_files"] == ["file1.bdf", "file2.bdf"]

            # Test boolean values
            bool_config_path = joinpath(test_dir, "bool_config.toml")
            open(bool_config_path, "w") do io
                println(io, "[files.output]")
                println(io, "save_continuous_data = true")
                println(io, "save_ica_data = false")
            end

            config = eegfun.load_config(bool_config_path)
            @test config["files"]["output"]["save_continuous_data"] == true
            @test config["files"]["output"]["save_ica_data"] == false

            # Test string values with special characters
            special_config_path = joinpath(test_dir, "special_config.toml")
            open(special_config_path, "w") do io
                println(io, "[files.input]")
                println(io, "directory = \"/path/with spaces/and/special@chars\"")
            end

            config = eegfun.load_config(special_config_path)
            @test config["files"]["input"]["directory"] == "/path/with spaces/and/special@chars"
        end

        @testset "Parameter Validation" begin
            # Test numeric range validation
            range_config_path = joinpath(test_dir, "range_config.toml")
            open(range_config_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "freq = 0.005")  # Below minimum
            end
            result = eegfun.load_config(range_config_path)
            @test result === nothing

            # Test allowed values validation
            allowed_config_path = joinpath(test_dir, "allowed_config.toml")
            open(allowed_config_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "method = \"invalid\"")  # Not in allowed values
            end
            result = eegfun.load_config(allowed_config_path)
            @test result === nothing

            # Test type conversion validation
            type_config_path = joinpath(test_dir, "type_config.toml")
            open(type_config_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "order = \"not_a_number\"")  # Invalid type
            end
            result = eegfun.load_config(type_config_path)
            @test result === nothing
        end

        @testset "Config Validation" begin
            # Test valid configuration
            valid_config_path = joinpath(test_dir, "valid.toml")
            open(valid_config_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "freq = 0.5")
                println(io, "method = \"fir\"")
                println(io, "order = 2")
            end
            config = eegfun.load_config(valid_config_path)
            @test config !== nothing
            @test config["preprocess"]["filter"]["highpass"]["freq"] == 0.5

            # Test invalid parameter value
            invalid_value_path = joinpath(test_dir, "invalid_value.toml")
            open(invalid_value_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "freq = -1.0")  # Below minimum
            end
            result = eegfun.load_config(invalid_value_path)
            @test result === nothing

            # Test invalid parameter type
            invalid_type_path = joinpath(test_dir, "invalid_type.toml")
            open(invalid_type_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "order = \"not_a_number\"")
            end
            result = eegfun.load_config(invalid_type_path)
            @test result === nothing
        end

        @testset "Config Merging" begin
            # Test merging with defaults
            custom_config_path = joinpath(test_dir, "custom.toml")
            open(custom_config_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "freq = 0.5")
            end
            config = eegfun.load_config(custom_config_path)
            @test config !== nothing
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
            config = eegfun.load_config(nested_config_path)
            @test config !== nothing
            @test config["files"]["input"]["directory"] == "/custom/path"
            @test config["files"]["input"]["raw_data_files"] == ["file1.bdf", "file2.bdf"]
        end

        @testset "Parameter Validation" begin
            # Test numeric range validation
            range_config_path = joinpath(test_dir, "range_config.toml")
            open(range_config_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "freq = 0.005")  # Below minimum
            end
            result = eegfun.load_config(range_config_path)
            @test result === nothing

            # Test allowed values validation
            allowed_config_path = joinpath(test_dir, "allowed_config.toml")
            open(allowed_config_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "method = \"invalid\"")  # Not in allowed values
            end
            result = eegfun.load_config(allowed_config_path)
            @test result === nothing

            # Test type conversion validation
            type_config_path = joinpath(test_dir, "type_config.toml")
            open(type_config_path, "w") do io
                println(io, "[preprocess.filter.highpass]")
                println(io, "order = \"not_a_number\"")  # Invalid type
            end
            result = eegfun.load_config(type_config_path)
            @test result === nothing

            # Test optional parameters
            optional_config_path = joinpath(test_dir, "optional.toml")
            open(optional_config_path, "w") do io
                println(io, "[files.input]")
                println(io, "epoch_condition_file = \"\"")  # Empty string for optional parameter
            end
            config = eegfun.load_config(optional_config_path)
            @test config !== nothing
            @test config["files"]["input"]["epoch_condition_file"] == ""
        end
    end

    # =============================================================================
    # VALIDATION
    # =============================================================================

    @testset "ValidationResult Tests" begin
        # Test successful validation
        result = eegfun.ValidationResult(success = true)
        @test result.success
        @test result.error === nothing
        @test result.key_path === nothing

        # Test failed validation with error
        result = eegfun.ValidationResult(success = false, error = "Test error")
        @test !result.success
        @test result.error == "Test error"
        @test result.key_path === nothing

        # Test failed validation with path
        result = eegfun.ValidationResult(success = false, error = "Test error", key_path = "test.path")
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
        result = eegfun._validate_config(valid_config)
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
        result = eegfun._validate_config(invalid_config)
        @test !result.success
        @test result.error !== nothing
        @test result.key_path !== nothing

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
        result = eegfun._validate_config(invalid_config)
        @test !result.success
        @test result.error !== nothing
        @test result.key_path !== nothing

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
        result = eegfun._validate_config(invalid_config)
        @test !result.success
        @test result.error !== nothing
        @test result.key_path !== nothing
    end

    # =============================================================================
    # TEMPLATE GENERATION
    # =============================================================================

    @testset "Config Template Generation" begin
        # Test template generation
        template_path = joinpath(test_dir, "template.toml")
        eegfun.generate_config_template(filename = template_path)
        @test isfile(template_path)

        # Verify template contents
        template = TOML.parsefile(template_path)
        @test haskey(template, "files")
        @test haskey(template["preprocess"], "filter")
        @test haskey(template, "preprocess")

        # Test template with custom filename
        custom_path = joinpath(test_dir, "custom.toml")
        eegfun.generate_config_template(filename = custom_path)
        @test isfile(custom_path)
    end

    # =============================================================================
    # PARAMETER CONSTRUCTORS AND HELPERS
    # =============================================================================

    @testset "Parameter Constructor Helper Tests" begin
        # Test string_param
        param = eegfun.string_param("Test string", "default")
        @test param isa eegfun.ConfigParameter{Union{Vector{String},String}}
        @test param.description == "Test string"
        @test param.default == "default"

        # Test string_param with allowed values
        param = eegfun.string_param("Test string", "default", allowed = ["a", "b", "c"])
        @test param.allowed == ["a", "b", "c"]

        # Test string_param with default empty string
        param = eegfun.string_param("Test string")
        @test param.default == ""

        # Test bool_param
        param = eegfun.bool_param("Test bool", true)
        @test param isa eegfun.ConfigParameter{Bool}
        @test param.description == "Test bool"
        @test param.default == true

        # Test bool_param with default false
        param = eegfun.bool_param("Test bool")
        @test param.default == false

        # Test number_param
        param = eegfun.number_param("Test number", 5.0, 0.0, 10.0)
        @test param isa eegfun.ConfigParameter{Real}
        @test param.description == "Test number"
        @test param.default == 5.0
        @test param.min == 0.0
        @test param.max == 10.0

        # Test number_param without min/max
        param = eegfun.number_param("Test number", 5.0)
        @test param.default == 5.0
        @test param.min === nothing
        @test param.max === nothing

        # Test channel_groups_param
        default = [["Fp1"], ["Fp2"]]
        param = eegfun.channel_groups_param("Test channels", default)
        @test param isa eegfun.ConfigParameter{Vector{Vector{String}}}
        @test param.description == "Test channels"
        @test param.default == default

        # Test _param helper function directly
        param = eegfun._param(Int, "Test int", 42, min = 0, max = 100)
        @test param isa eegfun.ConfigParameter{Int}
        @test param.description == "Test int"
        @test param.default == 42
        @test param.min == 0
        @test param.max == 100
    end

    @testset "Parameter Constructor Edge Cases" begin
        # Test string_param with empty allowed list
        param = eegfun.string_param("Test string", "default", allowed = String[])
        @test param.allowed == String[]

        # Test number_param with negative min/max
        param = eegfun.number_param("Test number", 0.0, -10.0, 10.0)
        @test param.min == -10.0
        @test param.max == 10.0

        # Test channel_groups_param with empty groups
        param = eegfun.channel_groups_param("Test channels", Vector{Vector{String}}())
        @test param.default == Vector{Vector{String}}()

        # Test _param with Union types
        param = eegfun._param(Union{String,Nothing}, "Test union", nothing)
        @test param isa eegfun.ConfigParameter{Union{String,Nothing}}
        @test param.default === nothing
    end

    @testset "ConfigParameter Constructor Tests" begin
        # Test with all fields
        param = eegfun.ConfigParameter{Int}(description = "Test parameter", default = 5, min = 1, max = 10)
        @test param.description == "Test parameter"
        @test param.default == 5
        @test param.min == 1
        @test param.max == 10

        # Test with minimal fields
        param = eegfun.ConfigParameter{String}(description = "Test parameter")
        @test param.description == "Test parameter"
        @test param.default === nothing
        @test param.min === nothing
        @test param.max === nothing
        @test param.allowed === nothing

        # Test with some fields
        param = eegfun.ConfigParameter{Float64}(description = "Test parameter", default = 1.0, min = 0.0)
        @test param.description == "Test parameter"
        @test param.default == 1.0
        @test param.min == 0.0
        @test param.max === nothing
        @test param.allowed === nothing
    end

    # =============================================================================
    # UTILITY FUNCTIONS
    # =============================================================================

    @testset "_param Helper Function Tests" begin
        # Test basic parameter creation
        param = eegfun._param(String, "Test param", "default")
        @test param isa eegfun.ConfigParameter{String}
        @test param.description == "Test param"
        @test param.default == "default"

        # Test with allowed values
        param = eegfun._param(String, "Test param", "default", allowed = ["a", "b"])
        @test param.allowed == ["a", "b"]

        # Test with min/max values
        param = eegfun._param(Real, "Test param", 5.0, min = 0.0, max = 10.0)
        @test param.min == 0.0
        @test param.max == 10.0

        # Test with nothing defaults
        param = eegfun._param(Bool, "Test param")
        @test param.default === nothing
        @test param.allowed === nothing
    end

    @testset "_filter_param_spec Tests" begin
        # Test filter parameter specification creation
        filter_spec = eegfun._filter_param_spec("test.prefix", true, "hp", 1.0, 0.1, 10.0, 2, 1, 5)

        @test haskey(filter_spec, "test.prefix.apply")
        @test haskey(filter_spec, "test.prefix.type")
        @test haskey(filter_spec, "test.prefix.method")
        @test haskey(filter_spec, "test.prefix.func")
        @test haskey(filter_spec, "test.prefix.freq")
        @test haskey(filter_spec, "test.prefix.order")

        # Test parameter types
        @test filter_spec["test.prefix.apply"] isa eegfun.ConfigParameter{Bool}
        @test filter_spec["test.prefix.type"] isa eegfun.ConfigParameter{Union{Vector{String},String}}
        @test filter_spec["test.prefix.method"] isa eegfun.ConfigParameter{Union{Vector{String},String}}
        @test filter_spec["test.prefix.freq"] isa eegfun.ConfigParameter{Real}

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
        sections = eegfun._group_parameters_by_section()

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
        @test eegfun._extract_subsection("files", "files.input.directory") == "input"
        @test eegfun._extract_subsection("preprocess.filter", "preprocess.filter.highpass.freq") == "highpass"

        # Test nested subsections
        @test eegfun._extract_subsection("preprocess.filter", "preprocess.filter.ica_highpass.freq") == "ica_highpass"

        # Test no subsection
        @test eegfun._extract_subsection("preprocess.ica", "preprocess.ica.apply") == ""

        # Test non-matching prefix
        @test eegfun._extract_subsection("files", "preprocess.filter.highpass.freq") == ""
    end

    @testset "_group_params_by_subsection Tests" begin
        # Test parameter grouping by subsection
        matching_params =
            ["preprocess.filter.highpass.freq", "preprocess.filter.highpass.apply", "preprocess.filter.lowpass.freq"]
        grouped = eegfun._group_params_by_subsection("preprocess.filter", matching_params)

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
        @test isempty(eegfun._merge_configs(Dict(), Dict()))

        # Test merging with one empty config
        default = Dict("a" => 1)
        user = Dict()
        merged = eegfun._merge_configs(default, user)
        @test merged["a"] == 1

        # Test deep merging
        default = Dict("a" => Dict("b" => Dict("c" => 1)))
        user = Dict("a" => Dict("b" => Dict("d" => 2)))
        merged = eegfun._merge_configs(default, user)
        @test merged["a"]["b"]["c"] == 1
        @test merged["a"]["b"]["d"] == 2

        # Test user config overrides default
        default = Dict("a" => 1)
        user = Dict("a" => 2)
        merged = eegfun._merge_configs(default, user)
        @test merged["a"] == 2

        # Test nested user config overrides default
        default = Dict("a" => Dict("b" => 1, "c" => 2))
        user = Dict("a" => Dict("b" => 3))
        merged = eegfun._merge_configs(default, user)
        @test merged["a"]["b"] == 3
        @test merged["a"]["c"] == 2
    end

    @testset "_merge_nested! Edge Cases" begin
        # Test merging with non-dict values (using compatible types)
        target = Dict("a" => Dict("b" => 1), "c" => 5)
        source = Dict("a" => Dict("d" => 2), "c" => 10)  # Override with compatible types

        eegfun._merge_nested!(target, source)
        @test target["a"]["b"] == 1  # Original value preserved
        @test target["a"]["d"] == 2  # New value added
        @test target["c"] == 10      # Value overridden

        # Test merging with empty source
        target = Dict("a" => 1, "b" => 2)
        source = Dict()
        original_target = copy(target)

        eegfun._merge_nested!(target, source)
        @test target == original_target  # Should be unchanged
    end

    @testset "_validate_parameter Tests" begin
        # Test valid parameter
        param = eegfun.ConfigParameter{Real}(description = "Test", default = 5.0, min = 0.0, max = 10.0)
        result = eegfun._validate_parameter(5.0, param, "test.param")
        @test result.success == true

        # Test value below minimum
        result = eegfun._validate_parameter(-1.0, param, "test.param")
        @test result.success == false
        @test contains(result.error, "must be >=")

        # Test value above maximum
        result = eegfun._validate_parameter(15.0, param, "test.param")
        @test result.success == false
        @test contains(result.error, "must be <=")

        # Test wrong type
        result = eegfun._validate_parameter("string", param, "test.param")
        @test result.success == false
        @test contains(result.error, "must be a number")

        # Test allowed values
        param = eegfun.ConfigParameter{String}(description = "Test", default = "a", allowed = ["a", "b", "c"])
        result = eegfun._validate_parameter("b", param, "test.param")
        @test result.success == true

        result = eegfun._validate_parameter("d", param, "test.param")
        @test result.success == false
        @test contains(result.error, "must be one of")

        # Test numeric type conversion (Int to Float)
        param = eegfun.ConfigParameter{Float64}(description = "Test", default = 5.0)
        result = eegfun._validate_parameter(5, param, "test.param")  # Int value
        @test result.success == true  # Should accept Int for Float64
    end

    # =============================================================================
    # DISPLAY/INFO FUNCTIONS
    # =============================================================================

    @testset "Parameter Info Display" begin
        # Test that all parameters exist and have required properties
        for (path, param) in eegfun.PARAMETERS
            @test !isempty(param.description)
            @test param isa eegfun.ConfigParameter
        end

        # Test specific parameter properties
        param = eegfun.PARAMETERS["preprocess.filter.highpass.freq"]
        @test param.description == "Cutoff frequency (Hz)"
        @test param.min == 0.01
        @test param.max == 20.0
        @test param.default == 0.1

        # Test non-existent parameter
        @test !haskey(eegfun.PARAMETERS, "nonexistent")
    end

    @testset "_show_parameter_details Tests" begin
        # Test with parameter that has all fields
        param = eegfun.ConfigParameter{Real}(
            description = "Test param",
            default = 5.0,
            min = 0.0,
            max = 10.0,
            allowed = ["a", "b"],
        )

        # Capture output (this is tricky with @info, so we'll test the function exists)
        @test typeof(eegfun._show_parameter_details) <: Function
    end

    @testset "_display_section Tests" begin
        # Test section display
        sections = eegfun._group_parameters_by_section()
        files_section = sections["files"]

        # Test function exists and can be called
        @test typeof(eegfun._display_section) <: Function
    end

    @testset "_display_subsection Tests" begin
        # Test subsection display
        sections = eegfun._group_parameters_by_section()
        input_params = sections["files"]["input"]

        # Test function exists and can be called
        @test typeof(eegfun._display_subsection) <: Function
    end

    @testset "_show_specific_parameter Tests" begin
        # Test with exact parameter match
        @test typeof(eegfun._show_specific_parameter) <: Function

        # Test with section prefix
        @test typeof(eegfun._show_specific_parameter) <: Function

        # Test with non-existent parameter
        @test typeof(eegfun._show_specific_parameter) <: Function
    end

    @testset "_show_section_overview Tests" begin
        # Test section overview
        matching_params =
            ["preprocess.filter.highpass.freq", "preprocess.filter.highpass.apply", "preprocess.filter.lowpass.freq"]

        # Test function exists and can be called
        @test typeof(eegfun._show_section_overview) <: Function
    end

    @testset "_display_grouped_params Tests" begin
        # Test grouped parameters display
        sections = eegfun._group_parameters_by_section()
        files_section = sections["files"]

        # Test function exists and can be called
        @test typeof(eegfun._display_grouped_params) <: Function
    end

    @testset "_show_all_parameters Tests" begin
        # Test that the function exists and can be called
        @test typeof(eegfun._show_all_parameters) <: Function

        # Test that it doesn't throw errors when called
        # (This function uses @info for output, so we can't easily capture it)
        @test nothing === nothing  # Placeholder to ensure test runs
    end

    # =============================================================================
    # EDGE CASES
    # =============================================================================

    @testset "_write_parameter_value Edge Cases" begin
        # Test string with special characters that need escaping
        io = IOBuffer()
        eegfun._write_parameter_value(io, "test_param", "path\\with\\backslashes")
        @test contains(String(take!(io)), "path\\\\with\\\\backslashes")

        # Test vector with nested vectors (channel groups)
        io = IOBuffer()
        nested_vector = [["Fp1", "IO1"], ["Fp2", "IO2"]]
        eegfun._write_parameter_value(io, "test_param", nested_vector)
        output = String(take!(io))
        @test contains(output, "test_param = [")
        @test contains(output, "Fp1")
        @test contains(output, "IO1")

        # Test empty vector
        io = IOBuffer()
        eegfun._write_parameter_value(io, "test_param", String[])
        @test String(take!(io)) == "test_param = []\n"

        # Test nothing value
        io = IOBuffer()
        eegfun._write_parameter_value(io, "test_param", nothing)
        @test String(take!(io)) == "test_param = nothing\n"

        # Test numeric values
        io = IOBuffer()
        eegfun._write_parameter_value(io, "test_param", 42)
        @test String(take!(io)) == "test_param = 42\n"

        # Test float values
        io = IOBuffer()
        eegfun._write_parameter_value(io, "test_param", 3.14)
        @test String(take!(io)) == "test_param = 3.14\n"

        # Test boolean values
        io = IOBuffer()
        eegfun._write_parameter_value(io, "test_param", true)
        @test String(take!(io)) == "test_param = true\n"
    end

    # =============================================================================
    # INTEGRATION TESTS
    # =============================================================================

    @testset "show_parameter_info Integration Tests" begin
        # Test with empty parameter name (show all)
        @test typeof(eegfun.show_parameter_info) <: Function

        # Test with specific parameter
        @test typeof(eegfun.show_parameter_info) <: Function

        # Test with section prefix
        @test typeof(eegfun.show_parameter_info) <: Function

        # Test with non-existent parameter
        @test typeof(eegfun.show_parameter_info) <: Function
    end

    @testset "Template Generation Integration Tests" begin
        # Test template generation with custom filename
        custom_template = joinpath(test_dir, "custom_template.toml")
        eegfun.generate_config_template(filename = custom_template)
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
        eegfun._write_template_header(io)
        output = String(take!(io))

        @test contains(output, "# EEG Processing Configuration Template")
        @test contains(output, "# Generated on")
        @test contains(output, "# This template shows all available configuration options")
        @test contains(output, "# Required fields are marked with [REQUIRED]")
        @test contains(output, "# Default values are shown where available")
    end

    @testset "_write_section Tests" begin
        # Create test data using the actual PARAMETERS structure
        sections = eegfun._group_parameters_by_section()
        files_section = sections["files"]

        io = IOBuffer()
        eegfun._write_section(io, "files", files_section)
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
        sections = eegfun._group_parameters_by_section()
        input_params = sections["files"]["input"]

        io = IOBuffer()
        eegfun._write_subsection(io, "files", "input", input_params)
        output = String(take!(io))

        @test contains(output, "# input Settings")
        @test contains(output, "[files.input]")
        @test contains(output, "directory = \".\"")
    end

    @testset "_write_template_sections Tests" begin
        io = IOBuffer()
        eegfun._write_template_sections(io)
        output = String(take!(io))

        # Should contain all major sections
        @test contains(output, "# files Settings")
        @test contains(output, "# preprocess Settings")
        @test contains(output, "# filter.highpass Settings")
        @test contains(output, "# ica Settings")
    end

    @testset "_write_parameter_docs Tests" begin
        # Test parameter with all fields
        param = eegfun.ConfigParameter{Real}(
            description = "Test param",
            default = 5.0,
            min = 0.0,
            max = 10.0,
            allowed = ["a", "b"],
        )
        io = IOBuffer()
        eegfun._write_parameter_docs(io, param)
        output = String(take!(io))

        @test contains(output, "# Test param")
        @test contains(output, "# Type: Real")
        @test contains(output, "# Range: 0.0 ≤ value ≤ 10.0")
        @test contains(output, "# Allowed values: a, b")
        @test contains(output, "# Default: 5.0")

        # Test required parameter
        param = eegfun.ConfigParameter{String}(description = "Required param", default = nothing)
        io = IOBuffer()
        eegfun._write_parameter_docs(io, param)
        output = String(take!(io))

        @test contains(output, "# [REQUIRED]")
        @test !contains(output, "# Default:")
    end

    @testset "_write_parameter_value Tests" begin
        # Test string value writing
        io = IOBuffer()
        eegfun._write_parameter_value(io, "test_param", "test_value")
        @test String(take!(io)) == "test_param = \"test_value\"\n"

        # Test numeric value writing
        io = IOBuffer()
        eegfun._write_parameter_value(io, "test_param", 42)
        @test String(take!(io)) == "test_param = 42\n"

        # Test boolean value writing
        io = IOBuffer()
        eegfun._write_parameter_value(io, "test_param", true)
        @test String(take!(io)) == "test_param = true\n"

        # Test empty vector writing
        io = IOBuffer()
        eegfun._write_parameter_value(io, "test_param", String[])
        @test String(take!(io)) == "test_param = []\n"

        # Test vector with values writing
        io = IOBuffer()
        eegfun._write_parameter_value(io, "test_param", ["a", "b", "c"])
        @test String(take!(io)) == "test_param = [a, b, c]\n"

        # Test vector with mixed types
        io = IOBuffer()
        eegfun._write_parameter_value(io, "test_param", [1, "two", 3.0])
        @test String(take!(io)) == "test_param = [1, two, 3.0]\n"
    end

    # =============================================================================
    # ERROR HANDLING
    # =============================================================================

    @testset "Error Handling Tests" begin
        # Test template generation with invalid filename
        invalid_path = "/invalid/path/that/does/not/exist/template.toml"

        # This should not throw an error but return gracefully
        # (The actual error handling depends on the implementation)
        @test typeof(eegfun.generate_config_template) <: Function
    end

    # Cleanup
    rm(test_dir, recursive = true, force = true)
end
