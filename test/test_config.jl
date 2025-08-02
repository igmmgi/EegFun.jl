using Test
using TOML
using Dates

# Import the package (which includes the config module)
using eegfun

@testset "Configuration System Tests" begin
    
    # Create temporary directory for test files
    test_dir = mktempdir()
    
    @testset "ConfigParameter Tests" begin
        # Test ConfigParameter struct creation
        param = eegfun.ConfigParameter{Float64}(
            description = "Test parameter",
            default = 1.0,
            min = 0.0,
            max = 2.0,
            allowed_values = [0.5, 1.0, 1.5]
        )
        @test param.description == "Test parameter"
        @test param.default == 1.0
        @test param.min == 0.0
        @test param.max == 2.0
        @test param.allowed_values == [0.5, 1.0, 1.5]

        # Test ConfigParameter with nothing values
        param = eegfun.ConfigParameter{String}(
            description = "Test parameter",
            default = nothing,
            min = nothing,
            max = nothing,
            allowed_values = nothing
        )
        @test param.description == "Test parameter"
        @test param.default === nothing
        @test param.min === nothing
        @test param.max === nothing
        @test param.allowed_values === nothing
    end

    @testset "PARAMETERS Dictionary Tests" begin
        # Test that all required parameters exist
        @test haskey(eegfun.PARAMETERS, "files.input.directory")
        @test haskey(eegfun.PARAMETERS, "files.input.raw_data_files")
        @test haskey(eegfun.PARAMETERS, "files.input.layout_file")
        @test haskey(eegfun.PARAMETERS, "filter.highpass.cutoff")
        @test haskey(eegfun.PARAMETERS, "filter.lowpass.cutoff")
        @test haskey(eegfun.PARAMETERS, "preprocess.epoch_start")
        @test haskey(eegfun.PARAMETERS, "preprocess.epoch_end")

        # Test parameter types
        @test eegfun.PARAMETERS["files.input.directory"] isa eegfun.ConfigParameter{String}
        @test eegfun.PARAMETERS["filter.highpass.cutoff"] isa eegfun.ConfigParameter{Real}
        @test eegfun.PARAMETERS["filter.highpass.order"] isa eegfun.ConfigParameter{Int}
        @test eegfun.PARAMETERS["files.output.save_continuous_data"] isa eegfun.ConfigParameter{Bool}
    end

    @testset "load_config Tests" begin
        
        @testset "Valid Configuration Loading" begin
            # Test 1: Load default config only - should work without errors
            default_config = eegfun.load_config(joinpath(dirname(@__FILE__), "..", "src", "config", "default.toml"))
            @test default_config isa Dict
            @test haskey(default_config, "preprocess")
            @test haskey(default_config, "filter") 
            
            # Test 2: Create a simple valid user config
            user_config_path = joinpath(test_dir, "valid_config.toml")
            open(user_config_path, "w") do io
                println(io, "[filter.highpass]")
                println(io, "cutoff = 1.0")
                println(io, "on = true")
                println(io, "")
                println(io, "[preprocess]")
                println(io, "epoch_start = -2.0")
                println(io, "epoch_end = 3.0")
            end
            
            config = eegfun.load_config(user_config_path)
            @test config isa Dict
            @test config["filter"]["highpass"]["cutoff"] == 1.0
            @test config["filter"]["highpass"]["on"] == true
            @test config["preprocess"]["epoch_start"] == -2.0
            @test config["preprocess"]["epoch_end"] == 3.0
            
            # Test 3: Verify default values are preserved when not overridden
            @test config["filter"]["lowpass"]["cutoff"] == 40  # default value
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
                println(io, "save_erp_data = false")
                println(io, "")
                println(io, "[ica.filter.highpass]")
                println(io, "cutoff = 2.5")
                println(io, "on = true")  # Add the "on" key that the test expects
            end
            
            config = eegfun.load_config(nested_config_path)
            @test config["files"]["output"]["directory"] == "/custom/output"
            @test config["files"]["output"]["save_erp_data"] == false
            @test config["files"]["output"]["save_ica_data"] == true  # default preserved
            @test config["ica"]["filter"]["highpass"]["cutoff"] == 2.5
            @test config["ica"]["filter"]["highpass"]["on"] == true  # now present in test config
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
                println(io, "[filter.highpass]")
                println(io, "cutoff = -5.0")  # Below minimum
            end
            result = eegfun.load_config(invalid_values_path)
            @test result === nothing
            
            # Test 9: Invalid parameter type - expect nothing to be returned
            invalid_type_path = joinpath(test_dir, "invalid_type.toml")
            open(invalid_type_path, "w") do io
                println(io, "[filter.highpass]")
                println(io, "cutoff = \"not_a_number\"")  # Wrong type
            end
            result = eegfun.load_config(invalid_type_path)
            @test result === nothing
            
            # Test 10: Invalid allowed values - expect nothing to be returned
            invalid_allowed_path = joinpath(test_dir, "invalid_allowed.toml")
            open(invalid_allowed_path, "w") do io
                println(io, "[filter.highpass]")
                println(io, "type = \"invalid_type\"")  # Not in allowed values
            end
            result = eegfun.load_config(invalid_allowed_path)
            @test result === nothing
        end
        
        @testset "Boundary Value Testing" begin
            # Test 11: Minimum boundary values
            min_boundary_path = joinpath(test_dir, "min_boundary.toml")
            open(min_boundary_path, "w") do io
                println(io, "[filter.highpass]")
                println(io, "cutoff = 0.01")  # Minimum allowed
                println(io, "order = 1")      # Minimum allowed
            end
            
            config = eegfun.load_config(min_boundary_path)
            @test config["filter"]["highpass"]["cutoff"] == 0.01
            @test config["filter"]["highpass"]["order"] == 1
            
            # Test 12: Maximum boundary values
            max_boundary_path = joinpath(test_dir, "max_boundary.toml")
            open(max_boundary_path, "w") do io
                println(io, "[filter.highpass]")
                println(io, "cutoff = 20.0")  # Maximum allowed
                println(io, "order = 8")      # Maximum allowed
            end
            
            config = eegfun.load_config(max_boundary_path)
            @test config["filter"]["highpass"]["cutoff"] == 20.0
            @test config["filter"]["highpass"]["order"] == 8
        end
        
        @testset "Data Type Conversions" begin
            # Test 13: Numeric type conversions
            conversion_path = joinpath(test_dir, "conversion.toml")
            open(conversion_path, "w") do io
                println(io, "[filter.highpass]")
                println(io, "cutoff = 1")     # Int should convert to Float
                println(io, "order = 2.0")    # Float should convert to Int
            end
            
            config = eegfun.load_config(conversion_path)
            @test config["filter"]["highpass"]["cutoff"] == 1.0  # Converted to Float
            @test config["filter"]["highpass"]["order"] == 2     # Converted to Int
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
            @test haskey(config, "filter")
            @test haskey(config, "ica")
            @test haskey(config, "preprocess")
            
            # Verify subsections exist
            @test haskey(config["files"], "input")
            @test haskey(config["files"], "output")
            @test haskey(config["filter"], "highpass")
            @test haskey(config["filter"], "lowpass")
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
                println(io, "[filter.highpass]")
            end
            
            config = eegfun.load_config(empty_sections_path)
            @test haskey(config, "files")
            @test haskey(config["files"], "input")
            
            # Test 16: Special characters in string values
            special_chars_path = joinpath(test_dir, "special_chars.toml")
            open(special_chars_path, "w") do io
                println(io, "[files.input]")
                println(io, "path = \"/path/with spaces/and/special@chars\"")
            end
            
            config = eegfun.load_config(special_chars_path)
            @test config["files"]["input"]["path"] == "/path/with spaces/and/special@chars"
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
                 println(io, "[filter.highpass]")
                 println(io, "cutoff = 0.005")  # Below minimum
             end
             result = eegfun.load_config(range_config_path)
             @test result === nothing

             # Test allowed values validation
             allowed_config_path = joinpath(test_dir, "allowed_config.toml")
             open(allowed_config_path, "w") do io
                 println(io, "[filter.highpass]")
                 println(io, "type = \"invalid\"")  # Not in allowed values
             end
             result = eegfun.load_config(allowed_config_path)
             @test result === nothing

             # Test type conversion validation
             type_config_path = joinpath(test_dir, "type_config.toml")
             open(type_config_path, "w") do io
                 println(io, "[filter.highpass]")
                 println(io, "order = \"not_a_number\"")  # Invalid type
             end
             result = eegfun.load_config(type_config_path)
             @test result === nothing
         end

         @testset "Config Validation" begin
             # Test valid configuration
             valid_config_path = joinpath(test_dir, "valid.toml")
             open(valid_config_path, "w") do io
                 println(io, "[filter.highpass]")
                 println(io, "cutoff = 0.5")
                 println(io, "type = \"fir\"")
                 println(io, "order = 2")
             end
             config = eegfun.load_config(valid_config_path)
             @test config !== nothing
             @test config["filter"]["highpass"]["cutoff"] == 0.5

             # Test invalid parameter value
             invalid_value_path = joinpath(test_dir, "invalid_value.toml")
             open(invalid_value_path, "w") do io
                 println(io, "[filter.highpass]")
                 println(io, "cutoff = -1.0")  # Below minimum
             end
             result = eegfun.load_config(invalid_value_path)
             @test result === nothing

             # Test invalid parameter type
             invalid_type_path = joinpath(test_dir, "invalid_type.toml")
             open(invalid_type_path, "w") do io
                 println(io, "[filter.highpass]")
                 println(io, "order = \"not_a_number\"")
             end
             result = eegfun.load_config(invalid_type_path)
             @test result === nothing
         end

         @testset "Config Merging" begin
             # Test merging with defaults
             custom_config_path = joinpath(test_dir, "custom.toml")
             open(custom_config_path, "w") do io
                 println(io, "[filter.highpass]")
                 println(io, "cutoff = 0.5")
             end
             config = eegfun.load_config(custom_config_path)
             @test config !== nothing
             @test config["filter"]["highpass"]["cutoff"] == 0.5
             @test config["filter"]["highpass"]["type"] == "fir"  # Default value
             @test config["filter"]["highpass"]["order"] == 1     # Default value

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
                 println(io, "[filter.highpass]")
                 println(io, "cutoff = 0.005")  # Below minimum
             end
             result = eegfun.load_config(range_config_path)
             @test result === nothing

             # Test allowed values validation
             allowed_config_path = joinpath(test_dir, "allowed_config.toml")
             open(allowed_config_path, "w") do io
                 println(io, "[filter.highpass]")
                 println(io, "type = \"invalid\"")  # Not in allowed values
             end
             result = eegfun.load_config(allowed_config_path)
             @test result === nothing

             # Test type conversion validation
             type_config_path = joinpath(test_dir, "type_config.toml")
             open(type_config_path, "w") do io
                 println(io, "[filter.highpass]")
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
    
     @testset "ValidationResult Tests" begin
         # Test successful validation
         result = eegfun.ValidationResult(success=true)
         @test result.success
         @test result.error === nothing
         @test result.path === nothing

         # Test failed validation with error
         result = eegfun.ValidationResult(success=false, error="Test error")
         @test !result.success
         @test result.error == "Test error"
         @test result.path === nothing

         # Test failed validation with path
         result = eegfun.ValidationResult(success=false, error="Test error", path="test.path")
         @test !result.success
         @test result.error == "Test error"
         @test result.path == "test.path"
     end
    
     @testset "Config Template Generation" begin
         # Test template generation
         template_path = joinpath(test_dir, "template.toml")
         eegfun.generate_config_template(template_path)
         @test isfile(template_path)

         # Verify template contents
         template = TOML.parsefile(template_path)
         @test haskey(template, "files")
         @test haskey(template, "filter")
         @test haskey(template, "preprocess")

         # Test template with custom filename
         custom_path = joinpath(test_dir, "custom.toml")
         eegfun.generate_config_template(custom_path)
         @test isfile(custom_path)
     end

     @testset "Parameter Info Display" begin
         # Test that all parameters exist and have required properties
         for (path, param) in eegfun.PARAMETERS
             @test !isempty(param.description)
             @test param isa eegfun.ConfigParameter
         end

         # Test specific parameter properties
         param = eegfun.PARAMETERS["filter.highpass.cutoff"]
         @test param.description == "High-pass filter cutoff frequency (Hz)"
         @test param.min == 0.01
         @test param.max == 20.0
         @test param.default == 0.1

         # Test non-existent parameter
         @test !haskey(eegfun.PARAMETERS, "nonexistent")
     end

     @testset "Config Value Extraction" begin
         # Test simple value extraction
         simple_config = Dict("key" => "value")
         result = eegfun._extract_values(simple_config)
         @test result["key"] == "value"

         # Test nested value extraction
         nested_config = Dict(
             "outer" => Dict(
                 "inner" => Dict(
                     "value" => 42
                 )
             )
         )
         result = eegfun._extract_values(nested_config)
         @test result["outer"]["inner"]["value"] == 42

         # Test mixed value extraction
         mixed_config = Dict(
             "simple" => "value",
             "nested" => Dict(
                 "number" => 42,
                 "array" => [1, 2, 3]
             )
         )
         result = eegfun._extract_values(mixed_config)
         @test result["simple"] == "value"
         @test result["nested"]["number"] == 42
         @test result["nested"]["array"] == [1, 2, 3]
     end

     @testset "Config Merging" begin
         # Test merging with empty config
         default = Dict("key" => "value")
         empty = Dict()
         result = eegfun._merge_configs(default, empty)
         @test result == default

         # Test merging with new values
         default = Dict("key" => "old")
         new = Dict("key" => "new")
         result = eegfun._merge_configs(default, new)
         @test result["key"] == "new"

         # Test deep merging
         default = Dict(
             "outer" => Dict(
                 "inner" => Dict(
                     "value" => "old"
                 )
             )
         )
         new = Dict(
             "outer" => Dict(
                 "inner" => Dict(
                     "value" => "new"
                 )
             )
         )
         result = eegfun._merge_configs(default, new)
         @test result["outer"]["inner"]["value"] == "new"

         # Test partial deep merging
         default = Dict(
             "outer" => Dict(
                 "inner1" => Dict("value" => "old1"),
                 "inner2" => Dict("value" => "old2")
             )
         )
         new = Dict(
             "outer" => Dict(
                 "inner1" => Dict("value" => "new1")
             )
         )
         result = eegfun._merge_configs(default, new)
         @test result["outer"]["inner1"]["value"] == "new1"
         @test result["outer"]["inner2"]["value"] == "old2"

         # Test array merging
         default = Dict("array" => [1, 2, 3])
         new = Dict("array" => [4, 5, 6])
         result = eegfun._merge_configs(default, new)
         @test result["array"] == [4, 5, 6]
     end
    
     @testset "_merge_configs" begin
         # Test simple merge
         default = Dict("a" => 1, "b" => 2)
         user = Dict("b" => 3, "c" => 4)
         merged = eegfun._merge_configs(default, user)
         @test merged["a"] == 1
         @test merged["b"] == 3
         @test merged["c"] == 4

         # Test nested merge
         default = Dict("a" => Dict("x" => 1, "y" => 2))
         user = Dict("a" => Dict("y" => 3, "z" => 4))
         merged = eegfun._merge_configs(default, user)
         @test merged["a"]["x"] == 1
         @test merged["a"]["y"] == 3
         @test merged["a"]["z"] == 4

         # Test deep nested merge
         default = Dict("a" => Dict("b" => Dict("c" => 1)))
         user = Dict("a" => Dict("b" => Dict("d" => 2)))
         merged = eegfun._merge_configs(default, user)
         @test merged["a"]["b"]["c"] == 1
         @test merged["a"]["b"]["d"] == 2
     end

     @testset "_extract_values Tests" begin
         # Test simple values
         config = Dict("a" => 1, "b" => "test")
         extracted = eegfun._extract_values(config)
         @test extracted["a"] == 1
         @test extracted["b"] == "test"

         # Test nested values
         config = Dict("a" => Dict("x" => 1, "y" => 2))
         extracted = eegfun._extract_values(config)
         @test extracted["a"]["x"] == 1
         @test extracted["a"]["y"] == 2

         # Test deep nested values
         config = Dict("a" => Dict("b" => Dict("c" => 1)))
         extracted = eegfun._extract_values(config)
         @test extracted["a"]["b"]["c"] == 1

         # Test empty dictionary
         config = Dict()
         extracted = eegfun._extract_values(config)
         @test isempty(extracted)
     end

     @testset "_validate_parameter Tests" begin
         # Test numeric validation
         param = eegfun.ConfigParameter{Int}(
             description = "Test parameter",
             min = 1,
             max = 10,
             default = 5
         )
         @test eegfun._validate_parameter(5, param, "test").success
         @test !eegfun._validate_parameter(0, param, "test").success
         @test !eegfun._validate_parameter(11, param, "test").success
         @test !eegfun._validate_parameter("5", param, "test").success

         # Test string validation with allowed values
         param = eegfun.ConfigParameter{String}(
             description = "Test parameter",
             allowed_values = ["a", "b", "c"],
             default = "a"
         )
         @test eegfun._validate_parameter("a", param, "test").success
         @test !eegfun._validate_parameter("d", param, "test").success
         @test !eegfun._validate_parameter(1, param, "test").success

         # Test vector validation
         param = eegfun.ConfigParameter{Vector{String}}(
             description = "Test parameter",
             default = ["a", "b"]
         )
         @test eegfun._validate_parameter(["a", "b"], param, "test").success
         @test !eegfun._validate_parameter(["a", 1], param, "test").success

         # Test type conversion
         param = eegfun.ConfigParameter{Float64}(
             description = "Test parameter",
             min = 1.0,
             max = 10.0,
             default = 5.0
         )
         @test eegfun._validate_parameter(5, param, "test").success  # Int to Float64
         @test !eegfun._validate_parameter("5.0", param, "test").success  # String not convertible

         # Test missing values
         param = eegfun.ConfigParameter{Union{String,Nothing}}(
             description = "Test parameter",
             default = nothing
         )
         @test eegfun._validate_parameter(nothing, param, "test").success
         @test eegfun._validate_parameter("test", param, "test").success
     end

     @testset "generate_config_template" begin
         # Test template generation
         template_file = "test_template.toml"
         eegfun.generate_config_template(template_file)
         @test isfile(template_file)
         
         # Read and verify template content
         content = read(template_file, String)
         @test contains(content, "# EEG Processing Configuration Template")
         @test contains(content, "[files]")
         @test contains(content, "[filter]")
         
         # Clean up
         rm(template_file)
     end
    
     @testset "ConfigParameter Constructor Tests" begin
         # Test with all fields
         param = eegfun.ConfigParameter{Int}(
             description = "Test parameter",
             default = 5,
             min = 1,
             max = 10,
             allowed_values = [1, 2, 3]
         )
         @test param.description == "Test parameter"
         @test param.default == 5
         @test param.min == 1
         @test param.max == 10
         @test param.allowed_values == [1, 2, 3]

         # Test with minimal fields
         param = eegfun.ConfigParameter{String}(
             description = "Test parameter"
         )
         @test param.description == "Test parameter"
         @test param.default === nothing
         @test param.min === nothing
         @test param.max === nothing
         @test param.allowed_values === nothing

         # Test with some fields
         param = eegfun.ConfigParameter{Float64}(
             description = "Test parameter",
             default = 1.0,
             min = 0.0
         )
         @test param.description == "Test parameter"
         @test param.default == 1.0
         @test param.min == 0.0
         @test param.max === nothing
         @test param.allowed_values === nothing
     end

     @testset "load_config Error Handling" begin
         # Test non-existent file
         result = eegfun.load_config("nonexistent.toml")
         @test result === nothing
     end

     @testset "Config Validation Tests" begin
         # Test valid configuration
         valid_config = Dict{String,Any}(
             "filter" => Dict{String,Any}(
                 "highpass" => Dict{String,Any}(
                     "on" => true,
                     "type" => "fir",
                     "cutoff" => 0.1,
                     "order" => 1
                 ),
                 "lowpass" => Dict{String,Any}(
                     "on" => true,
                     "type" => "fir",
                     "cutoff" => 40,
                     "order" => 3
                 )
             )
         )
         result = eegfun._validate_config(valid_config)
         @test result.success

         # Test invalid configuration - wrong type
         invalid_config = Dict{String,Any}(
             "filter" => Dict{String,Any}(
                 "highpass" => Dict{String,Any}(
                     "on" => "true",  # Should be Bool
                     "type" => "fir",
                     "cutoff" => 0.1,
                     "order" => 1
                 )
             )
         )
         result = eegfun._validate_config(invalid_config)
         @test !result.success
         @test result.error !== nothing
         @test result.path !== nothing

         # Test invalid configuration - out of range
         invalid_config = Dict{String,Any}(
             "filter" => Dict{String,Any}(
                 "highpass" => Dict{String,Any}(
                     "on" => true,
                     "type" => "fir",
                     "cutoff" => -1,  # Below minimum
                     "order" => 1
                 )
             )
         )
         result = eegfun._validate_config(invalid_config)
         @test !result.success
         @test result.error !== nothing
         @test result.path !== nothing

         # Test invalid configuration - wrong allowed value
         invalid_config = Dict{String,Any}(
             "filter" => Dict{String,Any}(
                 "highpass" => Dict{String,Any}(
                     "on" => true,
                     "type" => "invalid",  # Not in allowed_values
                     "cutoff" => 0.1,
                     "order" => 1
                 )
             )
         )
         result = eegfun._validate_config(invalid_config)
         @test !result.success
         @test result.error !== nothing
         @test result.path !== nothing
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
         default = Dict(
             "a" => Dict(
                 "b" => Dict("c" => 1)
             )
         )
         user = Dict(
             "a" => Dict(
                 "b" => Dict("d" => 2)
             )
         )
         merged = eegfun._merge_configs(default, user)
         @test merged["a"]["b"]["c"] == 1
         @test merged["a"]["b"]["d"] == 2

         # Test user config overrides default
         default = Dict("a" => 1)
         user = Dict("a" => 2)
         merged = eegfun._merge_configs(default, user)
         @test merged["a"] == 2

         # Test nested user config overrides default
         default = Dict(
             "a" => Dict(
                 "b" => 1,
                 "c" => 2
             )
         )
         user = Dict(
             "a" => Dict(
                 "b" => 3
             )
         )
         merged = eegfun._merge_configs(default, user)
         @test merged["a"]["b"] == 3
         @test merged["a"]["c"] == 2
     end

    # Cleanup
    rm(test_dir, recursive=true, force=true)
end 