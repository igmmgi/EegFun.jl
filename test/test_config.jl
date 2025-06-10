using Test
using TOML
using Dates

# Import the package (which includes the config module)
using eegfun

@testset "Configuration System Tests" begin
    
    # Create temporary directory for test files
    test_dir = mktempdir()
    
    @testset "load_config Tests" begin
        
        @testset "Valid Configuration Loading" begin
            # Test 1: Load default config only - should work without errors
            default_config = load_config(joinpath(dirname(@__FILE__), "..", "src", "config", "default.toml"))
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
            
            config = load_config(user_config_path)
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
            
            config = load_config(nested_config_path)
            @test config["files"]["output"]["directory"] == "/custom/output"
            @test config["files"]["output"]["save_erp_data"] == false
            @test config["files"]["output"]["save_ica_data"] == true  # default preserved
            @test config["ica"]["filter"]["highpass"]["cutoff"] == 2.5
            @test config["ica"]["filter"]["highpass"]["on"] == true  # now present in test config
        end
        
        @testset "Error Handling" begin
            # Test 6: Non-existent file - expect nothing to be returned
            result = load_config("nonexistent_file.toml")
            @test result === nothing
            
            # Test 7: Invalid TOML syntax - expect nothing to be returned
            invalid_toml_path = joinpath(test_dir, "invalid.toml")
            open(invalid_toml_path, "w") do io
                println(io, "[section")  # Missing closing bracket
                println(io, "key = value")
            end
            result = load_config(invalid_toml_path)
            @test result === nothing
            
            # Test 8: Invalid parameter values - expect nothing to be returned
            invalid_values_path = joinpath(test_dir, "invalid_values.toml")
            open(invalid_values_path, "w") do io
                println(io, "[filter.highpass]")
                println(io, "cutoff = -5.0")  # Below minimum
            end
            result = load_config(invalid_values_path)
            @test result === nothing
            
            # Test 9: Invalid parameter type - expect nothing to be returned
            invalid_type_path = joinpath(test_dir, "invalid_type.toml")
            open(invalid_type_path, "w") do io
                println(io, "[filter.highpass]")
                println(io, "cutoff = \"not_a_number\"")  # Wrong type
            end
            result = load_config(invalid_type_path)
            @test result === nothing
            
            # Test 10: Invalid allowed values - expect nothing to be returned
            invalid_allowed_path = joinpath(test_dir, "invalid_allowed.toml")
            open(invalid_allowed_path, "w") do io
                println(io, "[filter.highpass]")
                println(io, "type = \"invalid_type\"")  # Not in allowed values
            end
            result = load_config(invalid_allowed_path)
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
            
            config = load_config(min_boundary_path)
            @test config["filter"]["highpass"]["cutoff"] == 0.01
            @test config["filter"]["highpass"]["order"] == 1
            
            # Test 12: Maximum boundary values
            max_boundary_path = joinpath(test_dir, "max_boundary.toml")
            open(max_boundary_path, "w") do io
                println(io, "[filter.highpass]")
                println(io, "cutoff = 20.0")  # Maximum allowed
                println(io, "order = 8")      # Maximum allowed
            end
            
            config = load_config(max_boundary_path)
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
            
            config = load_config(conversion_path)
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
            
            config = load_config(complete_config_path)
            
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
            
            config = load_config(empty_sections_path)
            @test haskey(config, "files")
            @test haskey(config["files"], "input")
            
            # Test 16: Special characters in string values
            special_chars_path = joinpath(test_dir, "special_chars.toml")
            open(special_chars_path, "w") do io
                println(io, "[files.input]")
                println(io, "directory = \"/path/with spaces and symbols!@#\"")
            end
            
            config = load_config(special_chars_path)
            @test config["files"]["input"]["directory"] == "/path/with spaces and symbols!@#"
        end
    end
    
    @testset "ValidationResult Tests" begin
        # Test 17: ValidationResult constructor with keywords
        result = ValidationResult(success=true)
        @test result.success == true
        @test result.error === nothing
        @test result.path === nothing
        
        # Test 18: ValidationResult constructor with error
        result = ValidationResult(success=false, error="Test error", path="test.path")
        @test result.success == false
        @test result.error == "Test error"
        @test result.path == "test.path"
    end
    
    # Cleanup
    rm(test_dir, recursive=true, force=true)
end 