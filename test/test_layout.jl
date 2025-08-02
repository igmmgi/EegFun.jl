using Test
using DataFrames
using OrderedCollections

@testset "Layout Tests" begin
    # Read the BioSemi 64-channel layout
    layout_path = joinpath(@__DIR__, "..", "data", "layouts", "biosemi64.csv")
    test_layout = eegfun.read_layout(layout_path)
    
    @testset "Layout Reading" begin
        @test nrow(test_layout.data) == 64  # Should have 64 channels
        @test all(col -> col in propertynames(test_layout.data), [:label, :inc, :azi])  # Should have required columns
        @test all(test_layout.data.label .!== missing)  # No missing labels
        @test all(test_layout.data.inc .!== missing)    # No missing inclinations
        @test all(test_layout.data.azi .!== missing)    # No missing azimuths
        
        # Test error handling for read_layout
        result = eegfun.read_layout("nonexistent_file.csv")
        @test result === nothing
    end

    @testset "Polar to Cartesian Conversion" begin
        # Test 2D conversion
        layout_2d = copy(test_layout)
        eegfun.polar_to_cartesian_xy!(layout_2d)
        
        @test all(col -> col in propertynames(layout_2d.data), [:x2, :y2])
        @test nrow(layout_2d.data) == 64
        
        # Test specific coordinate calculations for known positions
        # Fp1 (inc=-92, azi=-72)
        fp1_idx = findfirst(test_layout.data.label .== :Fp1)
        @test isapprox(layout_2d.data[fp1_idx, :x2], -43.66, atol=0.01)  # Fp1 x
        @test isapprox(layout_2d.data[fp1_idx, :y2], 134.39, atol=0.01)  # Fp1 y
        
        # Cz (inc=0, azi=0)
        cz_idx = findfirst(test_layout.data.label .== :Cz)
        @test isapprox(layout_2d.data[cz_idx, :x2], 0.0, atol=0.01)     # Cz x
        @test isapprox(layout_2d.data[cz_idx, :y2], 0.0, atol=0.01)     # Cz y

        # Test 3D conversion
        layout_3d = copy(test_layout)
        eegfun.polar_to_cartesian_xyz!(layout_3d)
        
        @test all(col -> col in propertynames(layout_3d.data), [:x3, :y3, :z3])
        @test nrow(layout_3d.data) == 64
        
        # Test specific coordinate calculations for known positions
        # Fp1 (inc=-92, azi=-72)
        @test isapprox(layout_3d.data[fp1_idx, :x3], -27.18, atol=0.01)  # Fp1 x
        @test isapprox(layout_3d.data[fp1_idx, :y3], 83.64, atol=0.01)   # Fp1 y
        @test isapprox(layout_3d.data[fp1_idx, :z3], -3.07, atol=0.01)   # Fp1 z
        
        # Cz (inc=0, azi=0)
        @test isapprox(layout_3d.data[cz_idx, :x3], 0.0, atol=0.01)     # Cz x
        @test isapprox(layout_3d.data[cz_idx, :y3], 0.0, atol=0.01)     # Cz y
        @test isapprox(layout_3d.data[cz_idx, :z3], 88.0, atol=0.01)    # Cz z
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
        eegfun.get_layout_neighbours_xy!(layout_xy, 50.0)
        
        @test !isnothing(layout_xy.neighbours)
        @test haskey(layout_xy.neighbours, :Fp1)
        @test haskey(layout_xy.neighbours, :Cz)
        
        # Test that Fp1 has expected neighbours
        @test :AF7 in layout_xy.neighbours[:Fp1].electrodes
        @test :AF3 in layout_xy.neighbours[:Fp1].electrodes
        
        # Test weight calculations for 2D
        @test length(layout_xy.neighbours[:Fp1].weights) == length(layout_xy.neighbours[:Fp1].electrodes)
        @test isapprox(sum(layout_xy.neighbours[:Fp1].weights), 1.0, atol=1e-10)
        
        # Test 3D neighbours
        layout_xyz = copy(test_layout)
        eegfun.polar_to_cartesian_xyz!(layout_xyz)
        eegfun.get_layout_neighbours_xyz!(layout_xyz, 50.0)
        
        @test !isnothing(layout_xyz.neighbours)
        @test haskey(layout_xyz.neighbours, :Fp1)
        @test haskey(layout_xyz.neighbours, :Cz)
        
        # Test weight calculations for 3D
        @test length(layout_xyz.neighbours[:Fp1].weights) == length(layout_xyz.neighbours[:Fp1].electrodes)
        @test isapprox(sum(layout_xyz.neighbours[:Fp1].weights), 1.0, atol=1e-10)
        
        # Test error handling
        result = eegfun.get_layout_neighbours_xy!(test_layout, -1.0) # Negative distance
        @test result === nothing
    end

    @testset "Error Handling" begin
        # Test missing columns
        invalid_layout = eegfun.Layout(DataFrame(label = [:Fp1, :Fp2]), nothing, nothing)
        @test_throws ArgumentError eegfun.polar_to_cartesian_xy!(invalid_layout)
        @test_throws ArgumentError eegfun.polar_to_cartesian_xyz!(invalid_layout)
        
        # Test invalid data types
        invalid_types = eegfun.Layout(DataFrame(
            label = [:Fp1, :Fp2],
            inc = ["invalid", "invalid"],
            azi = ["invalid", "invalid"]
        ), nothing, nothing)
        @test_throws ArgumentError eegfun.polar_to_cartesian_xy!(invalid_types)
        @test_throws ArgumentError eegfun.polar_to_cartesian_xyz!(invalid_types)
    end

    @testset "Layout Functions" begin
        # Test read_layout
        @testset "read_layout" begin
            # Test reading a valid layout file
            layout_file = joinpath(@__DIR__, "..", "data", "layouts", "biosemi64.csv")
            layout = eegfun.read_layout(layout_file)
            @test size(layout.data, 1) > 0
            @test :label in propertynames(layout.data)
            @test :inc in propertynames(layout.data)
            @test :azi in propertynames(layout.data)
            
            # Test error handling for non-existent file
            result = eegfun.read_layout("nonexistent.csv")
            @test result === nothing
        end

        # Test polar to cartesian conversions
        @testset "polar_to_cartesian" begin
            # Create test layout
            layout = eegfun.Layout(DataFrame(
                label = [:Fp1, :Fp2, :F3, :F4],
                inc = [90.0, 90.0, 45.0, 45.0],
                azi = [0.0, 180.0, 45.0, 135.0]
            ), nothing, nothing)
            
            # Test xy conversion
            eegfun.polar_to_cartesian_xy!(layout)
            @test :x2 in propertynames(layout.data)
            @test :y2 in propertynames(layout.data)
            @test all(isfinite.(layout.data.x2))
            @test all(isfinite.(layout.data.y2))
            
            # Test xyz conversion
            eegfun.polar_to_cartesian_xyz!(layout)
            @test :x3 in propertynames(layout.data)
            @test :y3 in propertynames(layout.data)
            @test :z3 in propertynames(layout.data)
            @test all(isfinite.(layout.data.x3))
            @test all(isfinite.(layout.data.y3))
            @test all(isfinite.(layout.data.z3))
            
            # Test error handling for missing columns
            invalid_layout = eegfun.Layout(DataFrame(label = [:Fp1, :Fp2]), nothing, nothing)
            @test_throws ArgumentError eegfun.polar_to_cartesian_xy!(invalid_layout)
            @test_throws ArgumentError eegfun.polar_to_cartesian_xyz!(invalid_layout)
            
            # Test error handling for non-numeric values
            invalid_layout = eegfun.Layout(DataFrame(
                label = [:Fp1, :Fp2],
                inc = ["90", "90"],
                azi = ["0", "180"]
            ), nothing, nothing)
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
            layout = eegfun.Layout(DataFrame(
                label = [:Fp1, :Fp2, :F3, :F4],
                inc = [90.0, 90.0, 45.0, 45.0],
                azi = [0.0, 180.0, 45.0, 135.0]
            ), nothing, nothing)
            eegfun.polar_to_cartesian_xy!(layout)
            
            # Test xy neighbours
            eegfun.get_layout_neighbours_xy!(layout, 100.0)
            @test !isnothing(layout.neighbours)
            @test length(layout.neighbours) == size(layout.data, 1)
            
            # Test error handling for invalid distance criterion
            result = eegfun.get_layout_neighbours_xy!(layout, -1.0)
            @test result === nothing
            
            # Test xyz neighbours
            eegfun.polar_to_cartesian_xyz!(layout)
            eegfun.get_layout_neighbours_xyz!(layout, 100.0)
            @test !isnothing(layout.neighbours)
            @test length(layout.neighbours) == size(layout.data, 1)
            
            # Test error handling for invalid distance criterion
            result = eegfun.get_layout_neighbours_xyz!(layout, -1.0)
            @test result === nothing
        end
    end
end 