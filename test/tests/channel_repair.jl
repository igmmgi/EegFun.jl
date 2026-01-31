using Test
using DataFrames
using Statistics
using EegFun

@testset "channel_repair" begin

    # Test 1: Neighbor interpolation on ContinuousData
    @testset "neighbor_interpolation_continuous" begin
        dat = EegFun.create_test_continuous_data(n = 100, n_channels = 4)

        # Add some layout coordinates for neighbor calculation (closer together)
        dat.layout.data.x3 = [0.0, 0.1, 0.0, -0.1]
        dat.layout.data.y3 = [0.1, 0.0, -0.1, 0.0]
        dat.layout.data.z3 = [0.0, 0.0, 0.0, 0.0]

        # Calculate neighbors
        EegFun.get_neighbours_xyz!(dat.layout, 0.5)

        # Store original data
        original_data = copy(dat.data)

        # Repair channel Ch2
        EegFun.repair_channels!(dat, [:Ch2], method = :neighbor_interpolation)

        # Check that Ch2 data changed
        @test !isapprox(dat.data.Ch2, original_data.Ch2, rtol = 1e-10)

        # Check that other channels didn't change
        @test isapprox(dat.data.Ch1, original_data.Ch1, rtol = 1e-10)
        @test isapprox(dat.data.Ch3, original_data.Ch3, rtol = 1e-10)
        @test isapprox(dat.data.Ch4, original_data.Ch4, rtol = 1e-10)
    end

    # Test 2: Spherical spline on ContinuousData
    @testset "spherical_spline_continuous" begin
        dat = EegFun.create_test_continuous_data(n = 100, n_channels = 4)

        # Add some layout coordinates for spherical spline
        dat.layout.data.x3 = [0.0, 1.0, 0.0, -1.0]
        dat.layout.data.y3 = [1.0, 0.0, -1.0, 0.0]
        dat.layout.data.z3 = [0.0, 0.0, 0.0, 0.0]

        # Store original data
        original_data = copy(dat.data)

        # Repair channel Ch2
        EegFun.repair_channels!(dat, [:Ch2], method = :spherical_spline)

        # Check that Ch2 data changed
        @test !isapprox(dat.data.Ch2, original_data.Ch2, rtol = 1e-10)

        # Check that other channels didn't change
        @test isapprox(dat.data.Ch1, original_data.Ch1, rtol = 1e-10)
        @test isapprox(dat.data.Ch3, original_data.Ch3, rtol = 1e-10)
        @test isapprox(dat.data.Ch4, original_data.Ch4, rtol = 1e-10)
    end

    # Test 3: Non-mutating version
    @testset "non_mutating_version" begin
        dat = EegFun.create_test_continuous_data(n = 100, n_channels = 4)

        # Add some layout coordinates (closer together)
        dat.layout.data.x3 = [0.0, 0.1, 0.0, -0.1]
        dat.layout.data.y3 = [0.1, 0.0, -0.1, 0.0]
        dat.layout.data.z3 = [0.0, 0.0, 0.0, 0.0]

        # Calculate neighbors
        EegFun.get_neighbours_xyz!(dat.layout, 0.5)

        # Store original data
        original_data = copy(dat.data)

        # Use non-mutating version
        repaired_dat = EegFun.repair_channels(dat, [:Ch2], method = :neighbor_interpolation)

        # Check that original data is unchanged
        @test isapprox(dat.data.Ch2, original_data.Ch2, rtol = 1e-10)

        # Check that repaired data is different
        @test !isapprox(repaired_dat.data.Ch2, original_data.Ch2, rtol = 1e-10)

        # Check that other channels are the same
        @test isapprox(repaired_dat.data.Ch1, original_data.Ch1, rtol = 1e-10)
        @test isapprox(repaired_dat.data.Ch3, original_data.Ch3, rtol = 1e-10)
        @test isapprox(repaired_dat.data.Ch4, original_data.Ch4, rtol = 1e-10)
    end

    # Test 4: EpochData with neighbor interpolation
    @testset "neighbor_interpolation_epoch" begin
        dat = EegFun.create_test_epoch_data(n = 100, n_epochs = 3, n_channels = 4)

        # Add some layout coordinates (closer together)
        dat.layout.data.x3 = [0.0, 0.1, 0.0, -0.1]
        dat.layout.data.y3 = [0.1, 0.0, -0.1, 0.0]
        dat.layout.data.z3 = [0.0, 0.0, 0.0, 0.0]

        # Calculate neighbors
        EegFun.get_neighbours_xyz!(dat.layout, 0.5)

        # Store original data
        original_data = [copy(epoch) for epoch in dat.data]

        # Repair channel Ch2 in all epochs
        EegFun.repair_channels!(dat, [:Ch2], method = :neighbor_interpolation)

        # Check that Ch2 data changed in all epochs
        for i = 1:3
            @test !isapprox(dat.data[i].Ch2, original_data[i].Ch2, rtol = 1e-10)
        end

        # Check that other channels didn't change
        for i = 1:3
            @test isapprox(dat.data[i].Ch1, original_data[i].Ch1, rtol = 1e-10)
            @test isapprox(dat.data[i].Ch3, original_data[i].Ch3, rtol = 1e-10)
            @test isapprox(dat.data[i].Ch4, original_data[i].Ch4, rtol = 1e-10)
        end
    end

    # Test 5: EpochData with spherical spline
    @testset "spherical_spline_epoch" begin
        dat = EegFun.create_test_epoch_data(n = 100, n_epochs = 3, n_channels = 4)

        # Add some layout coordinates (closer together)
        dat.layout.data.x3 = [0.0, 0.1, 0.0, -0.1]
        dat.layout.data.y3 = [0.1, 0.0, -0.1, 0.0]
        dat.layout.data.z3 = [0.0, 0.0, 0.0, 0.0]

        # Store original data
        original_data = [copy(epoch) for epoch in dat.data]

        # Repair channel Ch2 in all epochs
        EegFun.repair_channels!(dat, [:Ch2], method = :spherical_spline)

        # Check that Ch2 data changed in all epochs
        for i = 1:3
            @test !isapprox(dat.data[i].Ch2, original_data[i].Ch2, rtol = 1e-10)
        end

        # Check that other channels didn't change
        for i = 1:3
            @test isapprox(dat.data[i].Ch1, original_data[i].Ch1, rtol = 1e-10)
            @test isapprox(dat.data[i].Ch3, original_data[i].Ch3, rtol = 1e-10)
            @test isapprox(dat.data[i].Ch4, original_data[i].Ch4, rtol = 1e-10)
        end
    end

    # Test 6: Multiple channels repair
    @testset "multiple_channels" begin
        dat = EegFun.create_test_continuous_data(n = 100, n_channels = 4)

        # Add some layout coordinates (closer together)
        dat.layout.data.x3 = [0.0, 0.1, 0.0, -0.1]
        dat.layout.data.y3 = [0.1, 0.0, -0.1, 0.0]
        dat.layout.data.z3 = [0.0, 0.0, 0.0, 0.0]

        # Calculate neighbors
        EegFun.get_neighbours_xyz!(dat.layout, 0.5)

        # Store original data
        original_data = copy(dat.data)

        # Repair multiple channels
        EegFun.repair_channels!(dat, [:Ch2, :Ch3], method = :neighbor_interpolation)

        # Check that repaired channels changed
        @test !isapprox(dat.data.Ch2, original_data.Ch2, rtol = 1e-10)
        @test !isapprox(dat.data.Ch3, original_data.Ch3, rtol = 1e-10)

        # Check that other channels didn't change
        @test isapprox(dat.data.Ch1, original_data.Ch1, rtol = 1e-10)
        @test isapprox(dat.data.Ch4, original_data.Ch4, rtol = 1e-10)
    end

    # Test 7: Error handling for unknown method
    @testset "error_handling" begin
        dat = EegFun.create_test_continuous_data(n = 100, n_channels = 4)

        # Test unknown method
        @test_throws ArgumentError EegFun.repair_channels!(dat, [:Ch2], method = :unknown_method)
    end

    # Test 8: Custom parameters for spherical spline
    @testset "custom_parameters" begin
        dat = EegFun.create_test_continuous_data(n = 100, n_channels = 4)

        # Add some layout coordinates (closer together)
        dat.layout.data.x3 = [0.0, 0.1, 0.0, -0.1]
        dat.layout.data.y3 = [0.1, 0.0, -0.1, 0.0]
        dat.layout.data.z3 = [0.0, 0.0, 0.0, 0.0]

        # Store original data
        original_data = copy(dat.data)

        # Repair with custom parameters
        EegFun.repair_channels!(dat, [:Ch2], method = :spherical_spline, m = 6, lambda = 1e-6)

        # Check that data changed
        @test !isapprox(dat.data.Ch2, original_data.Ch2, rtol = 1e-10)
    end

    # Test 9: Epoch selection for EpochData
    @testset "epoch_selection" begin
        dat = EegFun.create_test_epoch_data(n = 100, n_epochs = 3, n_channels = 4)

        # Add some layout coordinates (closer together)
        dat.layout.data.x3 = [0.0, 0.1, 0.0, -0.1]
        dat.layout.data.y3 = [0.1, 0.0, -0.1, 0.0]
        dat.layout.data.z3 = [0.0, 0.0, 0.0, 0.0]

        # Calculate neighbors
        EegFun.get_neighbours_xyz!(dat.layout, 0.5)

        # Store original data
        original_data = [copy(epoch) for epoch in dat.data]

        # Repair only first epoch
        EegFun.repair_channels!(dat, [:Ch2], method = :neighbor_interpolation, epoch_selection = EegFun.epochs(1))

        # Check that Ch2 data changed only in first epoch
        @test !isapprox(dat.data[1].Ch2, original_data[1].Ch2, rtol = 1e-10)
        @test isapprox(dat.data[2].Ch2, original_data[2].Ch2, rtol = 1e-10)
        @test isapprox(dat.data[3].Ch2, original_data[3].Ch2, rtol = 1e-10)
    end

end
