using Test
using DataFrames
using OrderedCollections

@testset "Layout Tests" begin
    # Read the BioSemi 64-channel layout
    layout_path = joinpath(@__DIR__, "..", "..", "resources", "layouts", "biosemi", "biosemi64.csv")
    test_layout = EegFun.read_layout(layout_path)

    @testset "Layout Reading" begin
        @test nrow(test_layout.data) == 64  # Should have 64 channels
        @test all(col -> col in propertynames(test_layout.data), [:label, :inc, :azi])  # Should have required columns
        @test all(test_layout.data.label .!== missing)  # No missing labels
        @test all(test_layout.data.inc .!== missing)    # No missing inclinations
        @test all(test_layout.data.azi .!== missing)    # No missing azimuths

        # Test error handling for read_layout
        @test_throws Exception EegFun.read_layout("nonexistent_file.csv")
    end

    @testset "Polar to Cartesian Conversion" begin
        # Test 2D conversion
        layout_2d = copy(test_layout)
        EegFun.polar_to_cartesian_xy!(layout_2d)

        @test all(col -> col in propertynames(layout_2d.data), [:x2, :y2])
        @test nrow(layout_2d.data) == 64

        # Test specific coordinate calculations for known positions
        # With preserve_radial_distance=true (default)
        # Fp1 (inc=-92, azi=-72)
        fp1_idx = findfirst(test_layout.data.label .== :Fp1)
        @test isapprox(layout_2d.data[fp1_idx, :x2], -0.316, atol = 0.01)  # Fp1 x
        @test isapprox(layout_2d.data[fp1_idx, :y2], 0.972, atol = 0.01)   # Fp1 y

        # Cz (inc=0, azi=0) - at vertex, preserves radial distance means stays at origin
        cz_idx = findfirst(test_layout.data.label .== :Cz)
        @test isapprox(layout_2d.data[cz_idx, :x2], 0.0, atol = 0.01)     # Cz x
        @test isapprox(layout_2d.data[cz_idx, :y2], 0.0, atol = 0.01)   # Cz y

        # Test 3D conversion
        layout_3d = copy(test_layout)
        EegFun.polar_to_cartesian_xyz!(layout_3d)

        @test all(col -> col in propertynames(layout_3d.data), [:x3, :y3, :z3])
        @test nrow(layout_3d.data) == 64

        # Test specific coordinate calculations for known positions
        # Fp1 (inc=-92, azi=-72) - normalized coordinates
        @test isapprox(layout_3d.data[fp1_idx, :x3], -0.309, atol = 0.01)  # Fp1 x
        @test isapprox(layout_3d.data[fp1_idx, :y3], 0.951, atol = 0.01)   # Fp1 y
        @test isapprox(layout_3d.data[fp1_idx, :z3], -0.324, atol = 0.01)  # Fp1 z

        # Cz (inc=0, azi=0) - normalized coordinates
        @test isapprox(layout_3d.data[cz_idx, :x3], 0.0, atol = 0.01)     # Cz x
        @test isapprox(layout_3d.data[cz_idx, :y3], 0.0, atol = 0.01)     # Cz y
        @test isapprox(layout_3d.data[cz_idx, :z3], 0.712, atol = 0.01)   # Cz z
    end

    @testset "Distance Calculations" begin
        # Test 2D distance
        @test isapprox(EegFun.distance_xy(0, 0, 3, 4), 5.0)
        @test isapprox(EegFun.squared_distance_xy(0, 0, 3, 4), 25.0)

        # Test 3D distance
        @test isapprox(EegFun.distance_xyz(0, 0, 0, 1, 2, 2), 3.0)
        @test isapprox(EegFun.squared_distance_xyz(0, 0, 0, 1, 2, 2), 9.0)

        # Test additional cases
        @test isapprox(EegFun.distance_xy(-1, -1, 1, 1), 2.8284271247461903)
        @test isapprox(EegFun.distance_xyz(-1, -1, -1, 1, 1, 1), 3.4641016151377544)

        # Test zero distance
        @test isapprox(EegFun.distance_xy(1, 1, 1, 1), 0.0)
        @test isapprox(EegFun.distance_xyz(1, 1, 1, 1, 1, 1), 0.0)
    end

    @testset "Electrode Neighbours" begin
        # Test 2D neighbours
        layout_xy = copy(test_layout)
        EegFun.polar_to_cartesian_xy!(layout_xy)
        EegFun.get_neighbours_xy!(layout_xy, 50.0)

        @test !isnothing(layout_xy.neighbours)
        @test haskey(layout_xy.neighbours, :Fp1)
        @test haskey(layout_xy.neighbours, :Cz)

        # Test that Fp1 has expected neighbours
        @test :AF7 in layout_xy.neighbours[:Fp1].channels
        @test :AF3 in layout_xy.neighbours[:Fp1].channels

        # Test weight calculations for 2D
        @test length(layout_xy.neighbours[:Fp1].weights) == length(layout_xy.neighbours[:Fp1].channels)
        @test isapprox(sum(layout_xy.neighbours[:Fp1].weights), 1.0, atol = 1e-10)

        # Test 3D neighbours
        layout_xyz = copy(test_layout)
        EegFun.polar_to_cartesian_xyz!(layout_xyz)
        EegFun.get_neighbours_xyz!(layout_xyz, 50.0)

        @test !isnothing(layout_xyz.neighbours)
        @test haskey(layout_xyz.neighbours, :Fp1)
        @test haskey(layout_xyz.neighbours, :Cz)

        # Test weight calculations for 3D
        @test length(layout_xyz.neighbours[:Fp1].weights) == length(layout_xyz.neighbours[:Fp1].channels)
        @test isapprox(sum(layout_xyz.neighbours[:Fp1].weights), 1.0, atol = 1e-10)

        # Test error handling
        result = EegFun.get_neighbours_xy!(test_layout, -1.0) # Negative distance
        @test result === nothing
    end

    @testset "Error Handling" begin
        # Test missing columns
        invalid_layout = EegFun.Layout(DataFrame(label = [:Fp1, :Fp2]), nothing, nothing)
        @test_throws ArgumentError EegFun.polar_to_cartesian_xy!(invalid_layout)
        @test_throws ArgumentError EegFun.polar_to_cartesian_xyz!(invalid_layout)

        # Test invalid data types
        invalid_types =
            EegFun.Layout(DataFrame(label = [:Fp1, :Fp2], inc = ["invalid", "invalid"], azi = ["invalid", "invalid"]), nothing, nothing)
        @test_throws ArgumentError EegFun.polar_to_cartesian_xy!(invalid_types)
        @test_throws ArgumentError EegFun.polar_to_cartesian_xyz!(invalid_types)
    end

    @testset "Layout Functions" begin
        # Test read_layout
        @testset "read_layout" begin
            # Test reading a valid layout file
            layout_file = joinpath(@__DIR__, "..", "..", "resources", "layouts", "biosemi", "biosemi64.csv")
            layout = EegFun.read_layout(layout_file)
            @test size(layout.data, 1) > 0
            @test :label in propertynames(layout.data)
            @test :inc in propertynames(layout.data)
            @test :azi in propertynames(layout.data)

            # Test error handling for non-existent file
            @test_throws Exception EegFun.read_layout("nonexistent.csv")
        end

        # Test polar to cartesian conversions
        @testset "polar_to_cartesian" begin
            # Create test layout
            layout = EegFun.Layout(
                DataFrame(label = [:Fp1, :Fp2, :F3, :F4], inc = [90.0, 90.0, 45.0, 45.0], azi = [0.0, 180.0, 45.0, 135.0]),
                nothing,
                nothing,
            )

            # Test xy conversion
            EegFun.polar_to_cartesian_xy!(layout)
            @test :x2 in propertynames(layout.data)
            @test :y2 in propertynames(layout.data)
            @test all(isfinite.(layout.data.x2))
            @test all(isfinite.(layout.data.y2))

            # Test xyz conversion
            EegFun.polar_to_cartesian_xyz!(layout)
            @test :x3 in propertynames(layout.data)
            @test :y3 in propertynames(layout.data)
            @test :z3 in propertynames(layout.data)
            @test all(isfinite.(layout.data.x3))
            @test all(isfinite.(layout.data.y3))
            @test all(isfinite.(layout.data.z3))

            # Test error handling for missing columns
            invalid_layout = EegFun.Layout(DataFrame(label = [:Fp1, :Fp2]), nothing, nothing)
            @test_throws ArgumentError EegFun.polar_to_cartesian_xy!(invalid_layout)
            @test_throws ArgumentError EegFun.polar_to_cartesian_xyz!(invalid_layout)

            # Test error handling for non-numeric values
            invalid_layout = EegFun.Layout(DataFrame(label = [:Fp1, :Fp2], inc = ["90", "90"], azi = ["0", "180"]), nothing, nothing)
            @test_throws ArgumentError EegFun.polar_to_cartesian_xy!(invalid_layout)
            @test_throws ArgumentError EegFun.polar_to_cartesian_xyz!(invalid_layout)
        end

        # Test distance calculations
        @testset "distance_calculations" begin
            # Test 2D distance
            @test EegFun.distance_xy(0, 0, 3, 4) ≈ 5.0
            @test EegFun.squared_distance_xy(0, 0, 3, 4) ≈ 25.0

            # Test 3D distance
            @test EegFun.distance_xyz(0, 0, 0, 1, 2, 2) ≈ 3.0
            @test EegFun.squared_distance_xyz(0, 0, 0, 1, 2, 2) ≈ 9.0

            # Test edge cases
            @test EegFun.distance_xy(0, 0, 0, 0) ≈ 0.0
            @test EegFun.distance_xyz(0, 0, 0, 0, 0, 0) ≈ 0.0
            @test EegFun.squared_distance_xy(0, 0, 0, 0) ≈ 0.0
            @test EegFun.squared_distance_xyz(0, 0, 0, 0, 0, 0) ≈ 0.0
        end

        # Test electrode neighbours
        @testset "electrode_neighbours" begin
            # Create test layout
            layout = EegFun.Layout(
                DataFrame(label = [:Fp1, :Fp2, :F3, :F4], inc = [90.0, 90.0, 45.0, 45.0], azi = [0.0, 180.0, 45.0, 135.0]),
                nothing,
                nothing,
            )
            EegFun.polar_to_cartesian_xy!(layout)

            # Test xy neighbours
            EegFun.get_neighbours_xy!(layout, 100.0)
            @test !isnothing(layout.neighbours)
            @test length(layout.neighbours) == size(layout.data, 1)

            # Test error handling for invalid distance criterion
            result = EegFun.get_neighbours_xy!(layout, -1.0)
            @test result === nothing

            # Test xyz neighbours
            EegFun.polar_to_cartesian_xyz!(layout)
            EegFun.get_neighbours_xyz!(layout, 100.0)
            @test !isnothing(layout.neighbours)
            @test length(layout.neighbours) == size(layout.data, 1)

            # Test error handling for invalid distance criterion
            result = EegFun.get_neighbours_xyz!(layout, -1.0)
            @test result === nothing
        end
    end

    @testset "Channel Renaming" begin
        # Create a simple test layout
        test_layout = EegFun.Layout(
            DataFrame(label = [:Fp1, :Fp2, :F3, :F4, :Cz], inc = [90.0, 90.0, 45.0, 45.0, 0.0], azi = [0.0, 180.0, 45.0, 135.0, 0.0]),
            nothing,
            nothing,
        )

        @testset "rename_channel! - Single channel" begin
            layout = copy(test_layout)
            original_labels = copy(layout.data.label)

            # Rename a single channel
            EegFun.rename_channel!(layout, Dict(:Fp1 => :Fpz))

            @test :Fpz in layout.data.label
            @test :Fp1 ∉ layout.data.label
            @test length(layout.data.label) == length(original_labels)

            # Check that other channels remain unchanged
            unchanged_channels = setdiff(original_labels, [:Fp1])
            for ch in unchanged_channels
                @test ch in layout.data.label
            end
        end

        @testset "rename_channel! - Multiple channels" begin
            layout = copy(test_layout)
            original_labels = copy(layout.data.label)

            # Rename multiple channels to unique names
            rename_dict = Dict(:Fp1 => :Fpz, :Fp2 => :Fp2_new, :F3 => :F3_new)
            EegFun.rename_channel!(layout, rename_dict)

            @test :Fpz in layout.data.label
            @test :Fp2_new in layout.data.label
            @test :F3_new in layout.data.label
            @test :Fp1 ∉ layout.data.label
            @test :Fp2 ∉ layout.data.label
            @test :F3 ∉ layout.data.label

            # Check that other channels remain unchanged
            unchanged_channels = setdiff(original_labels, [:Fp1, :Fp2, :F3])
            for ch in unchanged_channels
                @test ch in layout.data.label
            end
        end

        @testset "rename_channel! - Multiple channels to same name (prevented)" begin
            layout = copy(test_layout)

            # When multiple channels would be renamed to the same name, an error should be thrown
            rename_dict = Dict(:Fp1 => :Fpz, :Fp2 => :Fpz)
            @test_throws Any EegFun.rename_channel!(layout, rename_dict)

            # Layout should remain unchanged
            @test layout.data.label == test_layout.data.label
        end

        @testset "rename_channel! - Non-existent channels" begin
            layout = copy(test_layout)
            original_labels = copy(layout.data.label)

            # Try to rename channels that don't exist
            EegFun.rename_channel!(layout, Dict(:NonExistent => :NewName))

            # Layout should remain unchanged
            @test layout.data.label == original_labels
        end

        @testset "rename_channel! - Empty rename dict" begin
            layout = copy(test_layout)
            original_labels = copy(layout.data.label)

            # Empty rename dictionary should do nothing
            EegFun.rename_channel!(layout, Dict{Symbol,Symbol}())

            # Layout should remain unchanged
            @test layout.data.label == original_labels
        end

        @testset "rename_channel - Non-mutating version" begin
            layout = copy(test_layout)
            original_layout = copy(test_layout)

            # Use non-mutating version
            new_layout = EegFun.rename_channel(layout, Dict(:Fp1 => :Fpz))

            # Original should be unchanged
            @test original_layout.data.label == test_layout.data.label
            @test :Fp1 in original_layout.data.label
            @test :Fpz ∉ original_layout.data.label

            # New layout should have the changes
            @test :Fpz in new_layout.data.label
            @test :Fp1 ∉ new_layout.data.label
            @test length(new_layout.data.label) == length(original_layout.data.label)
        end

        @testset "rename_channel! - Neighbour cache clearing" begin
            layout = copy(test_layout)

            # Add some neighbours first
            EegFun.polar_to_cartesian_xy!(layout)
            EegFun.get_neighbours_xy!(layout, 100.0)

            # Verify neighbours exist
            @test !isnothing(layout.neighbours)
            @test haskey(layout.neighbours, :Fp1)

            # Rename a channel
            EegFun.rename_channel!(layout, Dict(:Fp1 => :Fpz))

            # Neighbours should be cleared
            @test isnothing(layout.neighbours)
        end

        @testset "rename_channel! - Coordinate preservation" begin
            layout = copy(test_layout)

            # Add coordinates first
            EegFun.polar_to_cartesian_xy!(layout)
            EegFun.polar_to_cartesian_xyz!(layout)

            # Store original coordinates for Fp1
            fp1_idx = findfirst(layout.data.label .== :Fp1)
            original_x2 = layout.data[fp1_idx, :x2]
            original_y2 = layout.data[fp1_idx, :y2]
            original_x3 = layout.data[fp1_idx, :x3]
            original_y3 = layout.data[fp1_idx, :y3]
            original_z3 = layout.data[fp1_idx, :z3]

            # Rename Fp1 to Fpz
            EegFun.rename_channel!(layout, Dict(:Fp1 => :Fpz))

            # Find the renamed channel
            fpz_idx = findfirst(layout.data.label .== :Fpz)

            # Coordinates should be preserved
            @test layout.data[fpz_idx, :x2] == original_x2
            @test layout.data[fpz_idx, :y2] == original_y2
            @test layout.data[fpz_idx, :x3] == original_x3
            @test layout.data[fpz_idx, :y3] == original_y3
            @test layout.data[fpz_idx, :z3] == original_z3
        end

        @testset "rename_channel! - Edge cases" begin
            layout = copy(test_layout)

            # Test renaming to the same name (should be a no-op)
            EegFun.rename_channel!(layout, Dict(:Fp1 => :Fp1))
            @test :Fp1 in layout.data.label

            # Test renaming all channels
            all_channels = copy(layout.data.label)
            rename_dict = Dict(old => Symbol("new_$(old)") for old in all_channels)
            EegFun.rename_channel!(layout, rename_dict)

            # All original names should be gone, all new names should exist
            for old in all_channels
                @test old ∉ layout.data.label
            end
            for old in all_channels
                @test Symbol("new_$(old)") in layout.data.label
            end
        end

        @testset "rename_channel! - Swap behavior" begin
            layout = copy(test_layout)
            original_labels = copy(layout.data.label)

            # Test swapping Fp1 and Fp2
            swap_dict = Dict(:Fp1 => :Fp2, :Fp2 => :Fp1)
            EegFun.rename_channel!(layout, swap_dict)

            # Should have swapped positions
            @test layout.data.label[1] == :Fp2  # First position now has Fp2
            @test layout.data.label[2] == :Fp1  # Second position now has Fp1

            # Other channels should remain unchanged
            @test layout.data.label[3] == :F3
            @test layout.data.label[4] == :F4
            @test layout.data.label[5] == :Cz

            # Total number of channels should remain the same
            @test length(layout.data.label) == 5

            # Verify it's a proper swap by checking the inverse operation
            inverse_swap = Dict(:Fp2 => :Fp1, :Fp1 => :Fp2)
            EegFun.rename_channel!(layout, inverse_swap)

            # Should be back to original
            @test layout.data.label == original_labels
        end
    end
end
