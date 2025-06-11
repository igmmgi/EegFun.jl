using Test
using DataFrames
using OrderedCollections

@testset "Layout Tests" begin
    # Read the BioSemi 64-channel layout
    layout_path = joinpath(@__DIR__, "..", "data", "layouts", "biosemi64.csv")
    test_layout = eegfun.read_layout(layout_path)
    
    @testset "Layout Reading" begin
        @test size(test_layout, 1) == 64  # Should have 64 channels
        @test all(col -> col in propertynames(test_layout), [:label, :inc, :azi])  # Should have required columns
        @test all(test_layout.label .!== missing)  # No missing labels
        @test all(test_layout.inc .!== missing)    # No missing inclinations
        @test all(test_layout.azi .!== missing)    # No missing azimuths
        
        # Test error handling for read_layout
        @test_throws SystemError eegfun.read_layout("nonexistent_file.csv")
    end

    @testset "Polar to Cartesian Conversion" begin
        # Test 2D conversion
        layout_2d = copy(test_layout)
        eegfun.polar_to_cartesian_xy!(layout_2d)
        
        @test all(col -> col in propertynames(layout_2d), [:x2, :y2])
        @test size(layout_2d, 1) == 64
        
        # Test specific coordinate calculations for known positions
        # Fp1 (inc=-92, azi=-72)
        fp1_idx = findfirst(test_layout.label .== :Fp1)
        @test isapprox(layout_2d[fp1_idx, :x2], -43.66, atol=0.01)  # Fp1 x
        @test isapprox(layout_2d[fp1_idx, :y2], 134.39, atol=0.01)  # Fp1 y
        
        # Cz (inc=0, azi=0)
        cz_idx = findfirst(test_layout.label .== :Cz)
        @test isapprox(layout_2d[cz_idx, :x2], 0.0, atol=0.01)     # Cz x
        @test isapprox(layout_2d[cz_idx, :y2], 0.0, atol=0.01)     # Cz y

        # Test 3D conversion
        layout_3d = copy(test_layout)
        eegfun.polar_to_cartesian_xyz!(layout_3d)
        
        @test all(col -> col in propertynames(layout_3d), [:x3, :y3, :z3])
        @test size(layout_3d, 1) == 64
        
        # Test specific coordinate calculations for known positions
        # Fp1 (inc=-92, azi=-72)
        @test isapprox(layout_3d[fp1_idx, :x3], -27.18, atol=0.01)  # Fp1 x
        @test isapprox(layout_3d[fp1_idx, :y3], 83.64, atol=0.01)   # Fp1 y
        @test isapprox(layout_3d[fp1_idx, :z3], -3.07, atol=0.01)   # Fp1 z
        
        # Cz (inc=0, azi=0)
        @test isapprox(layout_3d[cz_idx, :x3], 0.0, atol=0.01)     # Cz x
        @test isapprox(layout_3d[cz_idx, :y3], 0.0, atol=0.01)     # Cz y
        @test isapprox(layout_3d[cz_idx, :z3], 88.0, atol=0.01)    # Cz z
    end

    @testset "Distance Calculations" begin
        # Test 2D distance
        @test isapprox(eegfun.distance_xy(0, 0, 3, 4), 5.0)
        @test isapprox(eegfun.squared_distance_xy(0, 0, 3, 4), 25.0)
        
        # Test 3D distance
        @test isapprox(eegfun.distance_xyz(0, 0, 0, 1, 2, 2), 3.0)
        @test isapprox(eegfun.squared_distance_xyz(0, 0, 0, 1, 2, 2), 9.0)
        
        # Test additional cases
        @test isapprox(eegfun.distance_xy(-1, -1, 1, 1), 2.8284271247461903)
        @test isapprox(eegfun.distance_xyz(-1, -1, -1, 1, 1, 1), 3.4641016151377544)
        
        # Test zero distance
        @test isapprox(eegfun.distance_xy(1, 1, 1, 1), 0.0)
        @test isapprox(eegfun.distance_xyz(1, 1, 1, 1, 1, 1), 0.0)
    end

    @testset "Electrode Neighbours" begin
        # Test 2D neighbours
        layout_xy = copy(test_layout)
        eegfun.polar_to_cartesian_xy!(layout_xy)
        neighbours_2d = eegfun.get_electrode_neighbours_xy(layout_xy, 50.0)
        
        @test neighbours_2d isa OrderedDict{Symbol, eegfun.Neighbours}
        @test haskey(neighbours_2d, :Fp1)
        @test haskey(neighbours_2d, :Cz)
        
        # Test that Fp1 has expected neighbours
        @test :AF7 in neighbours_2d[:Fp1].electrodes
        @test :AF3 in neighbours_2d[:Fp1].electrodes
        
        # Test weight calculations for 2D
        @test length(neighbours_2d[:Fp1].weights) == length(neighbours_2d[:Fp1].electrodes)
        @test isapprox(sum(neighbours_2d[:Fp1].weights), 1.0, atol=1e-10)
        
        # Test 3D neighbours
        layout_xyz = copy(test_layout)
        eegfun.polar_to_cartesian_xyz!(layout_xyz)
        neighbours_3d = eegfun.get_electrode_neighbours_xyz(layout_xyz, 50.0)
        
        @test neighbours_3d isa OrderedDict{Symbol, eegfun.Neighbours}
        @test haskey(neighbours_3d, :Fp1)
        @test haskey(neighbours_3d, :Cz)
        
        # Test weight calculations for 3D
        @test length(neighbours_3d[:Fp1].weights) == length(neighbours_3d[:Fp1].electrodes)
        @test isapprox(sum(neighbours_3d[:Fp1].weights), 1.0, atol=1e-10)
        
        # Test error handling
        @test_throws ArgumentError eegfun.get_electrode_neighbours_xy(test_layout, 50.0)  # Missing x2, y2
        @test_throws ArgumentError eegfun.get_electrode_neighbours_xyz(test_layout, 50.0)  # Missing x3, y3, z3
        @test_throws ArgumentError eegfun.get_electrode_neighbours_xy(layout_xy, -1.0)   # Negative distance
        @test_throws ArgumentError eegfun.get_electrode_neighbours_xyz(layout_xyz, -1.0) # Negative distance
    end

    @testset "Error Handling" begin
        # Test missing columns
        invalid_layout = DataFrame(label = [:Fp1, :Fp2])
        @test_throws ArgumentError eegfun.polar_to_cartesian_xy!(invalid_layout)
        @test_throws ArgumentError eegfun.polar_to_cartesian_xyz!(invalid_layout)
        
        # Test invalid data types
        invalid_types = DataFrame(
            label = [:Fp1, :Fp2],
            inc = ["invalid", "invalid"],
            azi = ["invalid", "invalid"]
        )
        @test_throws ArgumentError eegfun.polar_to_cartesian_xy!(invalid_types)
        @test_throws ArgumentError eegfun.polar_to_cartesian_xyz!(invalid_types)
    end

    @testset "Layout Functions" begin
        # Test read_layout
        @testset "read_layout" begin
            # Test reading a valid layout file
            layout_file = joinpath(@__DIR__, "..", "data", "layouts", "biosemi64.csv")
            layout = eegfun.read_layout(layout_file)
            @test size(layout, 1) > 0
            @test :label in propertynames(layout)
            @test :inc in propertynames(layout)
            @test :azi in propertynames(layout)
            
            # Test error handling for non-existent file
            @test_throws SystemError eegfun.read_layout("nonexistent.csv")
        end

        # Test polar to cartesian conversions
        @testset "polar_to_cartesian" begin
            # Create test layout
            layout = DataFrame(
                label = [:Fp1, :Fp2, :F3, :F4],
                inc = [90.0, 90.0, 45.0, 45.0],
                azi = [0.0, 180.0, 45.0, 135.0]
            )
            
            # Test xy conversion
            eegfun.polar_to_cartesian_xy!(layout)
            @test :x2 in propertynames(layout)
            @test :y2 in propertynames(layout)
            @test all(isfinite.(layout.x2))
            @test all(isfinite.(layout.y2))
            
            # Test xyz conversion
            eegfun.polar_to_cartesian_xyz!(layout)
            @test :x3 in propertynames(layout)
            @test :y3 in propertynames(layout)
            @test :z3 in propertynames(layout)
            @test all(isfinite.(layout.x3))
            @test all(isfinite.(layout.y3))
            @test all(isfinite.(layout.z3))
            
            # Test error handling for missing columns
            invalid_layout = DataFrame(label = [:Fp1, :Fp2])
            @test_throws ArgumentError eegfun.polar_to_cartesian_xy!(invalid_layout)
            @test_throws ArgumentError eegfun.polar_to_cartesian_xyz!(invalid_layout)
            
            # Test error handling for non-numeric values
            invalid_layout = DataFrame(
                label = [:Fp1, :Fp2],
                inc = ["90", "90"],
                azi = ["0", "180"]
            )
            @test_throws ArgumentError eegfun.polar_to_cartesian_xy!(invalid_layout)
            @test_throws ArgumentError eegfun.polar_to_cartesian_xyz!(invalid_layout)
        end

        # Test distance calculations
        @testset "distance_calculations" begin
            # Test 2D distance
            @test eegfun.distance_xy(0, 0, 3, 4) ≈ 5.0
            @test eegfun.squared_distance_xy(0, 0, 3, 4) ≈ 25.0
            
            # Test 3D distance
            @test eegfun.distance_xyz(0, 0, 0, 1, 2, 2) ≈ 3.0
            @test eegfun.squared_distance_xyz(0, 0, 0, 1, 2, 2) ≈ 9.0
            
            # Test edge cases
            @test eegfun.distance_xy(0, 0, 0, 0) ≈ 0.0
            @test eegfun.distance_xyz(0, 0, 0, 0, 0, 0) ≈ 0.0
            @test eegfun.squared_distance_xy(0, 0, 0, 0) ≈ 0.0
            @test eegfun.squared_distance_xyz(0, 0, 0, 0, 0, 0) ≈ 0.0
        end

        # Test electrode neighbours
        @testset "electrode_neighbours" begin
            # Create test layout
            layout = DataFrame(
                label = [:Fp1, :Fp2, :F3, :F4],
                inc = [90.0, 90.0, 45.0, 45.0],
                azi = [0.0, 180.0, 45.0, 135.0]
            )
            eegfun.polar_to_cartesian_xy!(layout)
            
            # Test xy neighbours
            neighbours = eegfun.get_electrode_neighbours_xy(layout, 100.0)
            @test neighbours isa OrderedDict{Symbol, eegfun.Neighbours}
            @test length(neighbours) == size(layout, 1)
            
            # Test error handling for missing columns
            invalid_layout = DataFrame(label = [:Fp1, :Fp2])
            @test_throws ArgumentError eegfun.get_electrode_neighbours_xy(invalid_layout, 100.0)
            
            # Test error handling for invalid distance criterion
            @test_throws ArgumentError eegfun.get_electrode_neighbours_xy(layout, -1.0)
            
            # Test xyz neighbours
            eegfun.polar_to_cartesian_xyz!(layout)
            neighbours = eegfun.get_electrode_neighbours_xyz(layout, 100.0)
            @test neighbours isa OrderedDict{Symbol, eegfun.Neighbours}
            @test length(neighbours) == size(layout, 1)
            
            # Test error handling for missing columns
            invalid_layout = DataFrame(label = [:Fp1, :Fp2])
            @test_throws ArgumentError eegfun.get_electrode_neighbours_xyz(invalid_layout, 100.0)
            
            # Test error handling for invalid distance criterion
            @test_throws ArgumentError eegfun.get_electrode_neighbours_xyz(layout, -1.0)
        end
    end
end 