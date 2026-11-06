using Test
using DataFrames
using Statistics
using EegFun

@testset "channel_repair" begin

    # Neighbor interpolation on ContinuousData
    @testset "neighbor_interpolation_continuous" begin
        dat = EegFun.create_test_continuous_data(n = 100, n_channels = 4)

        # Add some layout coordinates for neighbor calculation (closer together)
        dat.layout.data.x3 = [0.0, 0.1, 0.0, -0.1]
        dat.layout.data.y3 = [0.1, 0.0, -0.1, 0.0]
        dat.layout.data.z3 = [0.0, 0.0, 0.0, 0.0]

        # Calculate neighbors
        EegFun.get_neighbors_xyz!(dat.layout, 0.5)

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

    # Spherical spline on ContinuousData
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

    # Non-mutating version
    @testset "non_mutating_version" begin
        dat = EegFun.create_test_continuous_data(n = 100, n_channels = 4)

        # Add some layout coordinates (closer together)
        dat.layout.data.x3 = [0.0, 0.1, 0.0, -0.1]
        dat.layout.data.y3 = [0.1, 0.0, -0.1, 0.0]
        dat.layout.data.z3 = [0.0, 0.0, 0.0, 0.0]

        # Calculate neighbors
        EegFun.get_neighbors_xyz!(dat.layout, 0.5)

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

    # EpochData with neighbor interpolation
    @testset "neighbor_interpolation_epoch" begin
        dat = EegFun.create_test_epoch_data(n = 100, n_epochs = 3, n_channels = 4)

        # Add some layout coordinates (closer together)
        dat.layout.data.x3 = [0.0, 0.1, 0.0, -0.1]
        dat.layout.data.y3 = [0.1, 0.0, -0.1, 0.0]
        dat.layout.data.z3 = [0.0, 0.0, 0.0, 0.0]

        # Calculate neighbors
        EegFun.get_neighbors_xyz!(dat.layout, 0.5)

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

    # EpochData with spherical spline
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

    # Multiple channels repair
    @testset "multiple_channels" begin
        dat = EegFun.create_test_continuous_data(n = 100, n_channels = 4)

        # Add some layout coordinates (closer together)
        dat.layout.data.x3 = [0.0, 0.1, 0.0, -0.1]
        dat.layout.data.y3 = [0.1, 0.0, -0.1, 0.0]
        dat.layout.data.z3 = [0.0, 0.0, 0.0, 0.0]

        # Calculate neighbors
        EegFun.get_neighbors_xyz!(dat.layout, 0.5)

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

    # Error handling for unknown method
    @testset "error_handling" begin
        dat = EegFun.create_test_continuous_data(n = 100, n_channels = 4)

        # Test unknown method
        @test_throws EegFun.EegFunError EegFun.repair_channels!(dat, [:Ch2], method = :unknown_method)
    end

    # Custom parameters for spherical spline
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

    # Epoch selection for EpochData
    @testset "epoch_selection" begin
        dat = EegFun.create_test_epoch_data(n = 100, n_epochs = 3, n_channels = 4)

        # Add some layout coordinates (closer together)
        dat.layout.data.x3 = [0.0, 0.1, 0.0, -0.1]
        dat.layout.data.y3 = [0.1, 0.0, -0.1, 0.0]
        dat.layout.data.z3 = [0.0, 0.0, 0.0, 0.0]

        # Calculate neighbors
        EegFun.get_neighbors_xyz!(dat.layout, 0.5)

        # Store original data
        original_data = [copy(epoch) for epoch in dat.data]

        # Repair only first epoch
        EegFun.repair_channels!(dat, [:Ch2], method = :neighbor_interpolation, epoch_selection = EegFun.epochs(1))

        # Check that Ch2 data changed only in first epoch
        @test !isapprox(dat.data[1].Ch2, original_data[1].Ch2, rtol = 1e-10)
        @test isapprox(dat.data[2].Ch2, original_data[2].Ch2, rtol = 1e-10)
        @test isapprox(dat.data[3].Ch2, original_data[3].Ch2, rtol = 1e-10)
    end

    @testset "Vector{<:EegData} methods" begin
        dat1 = EegFun.create_test_continuous_data(n = 100, n_channels = 4)
        dat2 = EegFun.create_test_continuous_data(n = 100, n_channels = 4)
        for d in [dat1, dat2]
            d.layout.data.x3 = [0.0, 0.1, 0.0, -0.1]
            d.layout.data.y3 = [0.1, 0.0, -0.1, 0.0]
            d.layout.data.z3 = [0.0, 0.0, 0.0, 0.0]
            EegFun.get_neighbors_xyz!(d.layout, 0.5)
        end
        orig1_Ch2 = copy(dat1.data.Ch2)
        orig2_Ch2 = copy(dat2.data.Ch2)

        # non-mutating
        repaired_vec = EegFun.repair_channels([dat1, dat2], [:Ch2], method = :neighbor_interpolation)
        @test isapprox(dat1.data.Ch2, orig1_Ch2, rtol = 1e-10)
        @test !isapprox(repaired_vec[1].data.Ch2, orig1_Ch2, rtol = 1e-10)
        @test !isapprox(repaired_vec[2].data.Ch2, orig2_Ch2, rtol = 1e-10)

        # mutating
        EegFun.repair_channels!([dat1, dat2], [:Ch2], method = :neighbor_interpolation)
        @test !isapprox(dat1.data.Ch2, orig1_Ch2, rtol = 1e-10)
        @test !isapprox(dat2.data.Ch2, orig2_Ch2, rtol = 1e-10)
    end

    @testset "String channels and AbstractVector views" begin
        dat1 = EegFun.create_test_continuous_data(n = 100, n_channels = 4)
        dat1.layout.data.x3 = [0.0, 0.1, 0.0, -0.1]
        dat1.layout.data.y3 = [0.1, 0.0, -0.1, 0.0]
        dat1.layout.data.z3 = [0.0, 0.0, 0.0, 0.0]
        EegFun.get_neighbors_xyz!(dat1.layout, 0.5)

        orig_Ch2 = copy(dat1.data.Ch2)
        # Test string channel "Ch2"
        repaired = EegFun.repair_channels(dat1, "Ch2", method = :neighbor_interpolation)
        @test !isapprox(repaired.data.Ch2, orig_Ch2, rtol = 1e-10)

        # Test view dispatch
        dats = [copy(dat1), copy(dat1), copy(dat1)]
        v = @view dats[1:2]
        EegFun.repair_channels!(v, ["Ch2"], method = :neighbor_interpolation)
        @test !isapprox(dats[1].data.Ch2, orig_Ch2, rtol = 1e-10)
        @test !isapprox(dats[2].data.Ch2, orig_Ch2, rtol = 1e-10)
        @test isapprox(dats[3].data.Ch2, orig_Ch2, rtol = 1e-10)
    end

    @testset "Order-independent neighbor repair, bad-neighbor exclusion, ErpData, and analysis_info" begin
        # 1. Exclude co-bad channels from neighbor weights independent of channel order
        df_rep = DataFrame(
            time = [0.0, 0.01, 0.02],
            sample = [1, 2, 3],
            Ch1 = [1000.0, 1000.0, 1000.0], # corrupted bad channel
            Ch2 = [-500.0, -500.0, -500.0], # corrupted bad channel (neighbor of Ch1!)
            Ch3 = [10.0, 20.0, 30.0],       # clean neighbor
            Ch4 = [30.0, 40.0, 50.0],       # clean neighbor
        )
        nb_dict = EegFun.OrderedDict(
            :Ch1 => EegFun.Neighbors([:Ch2, :Ch3, :Ch4], [1.0, 1.0, 1.0], [0.5, 0.25, 0.25]),
            :Ch2 => EegFun.Neighbors([:Ch1, :Ch3, :Ch4], [1.0, 1.0, 1.0], [0.5, 0.25, 0.25]),
            :Ch3 => EegFun.Neighbors([:Ch1, :Ch2, :Ch4], [1.0, 1.0, 1.0], [1 / 3, 1 / 3, 1 / 3]),
            :Ch4 => EegFun.Neighbors([:Ch1, :Ch2, :Ch3], [1.0, 1.0, 1.0], [1 / 3, 1 / 3, 1 / 3]),
        )
        lay_rep = EegFun.Layout(DataFrame(label = [:Ch1, :Ch2, :Ch3, :Ch4], inc = zeros(4), azi = zeros(4)), nb_dict, 0.5, :xyz)
        cont_ab = EegFun.ContinuousData("rep_ab", copy(df_rep), lay_rep, 100, EegFun.AnalysisInfo())
        cont_ba = EegFun.ContinuousData("rep_ba", copy(df_rep), lay_rep, 100, EegFun.AnalysisInfo())

        EegFun.repair_channels!(cont_ab, [:Ch1, :Ch2]; method = :neighbor_interpolation)
        EegFun.repair_channels!(cont_ba, [:Ch2, :Ch1]; method = :neighbor_interpolation)

        # Equal re-normalized clean weights (0.25/(0.25+0.25) = 0.5 each on Ch3 and Ch4) -> exact mean [20, 30, 40]
        @test cont_ab.data.Ch1 ≈ [20.0, 30.0, 40.0]
        @test cont_ab.data.Ch2 ≈ [20.0, 30.0, 40.0]
        # Order-independent!
        @test cont_ab.data.Ch1 ≈ cont_ba.data.Ch1
        @test cont_ab.data.Ch2 ≈ cont_ba.data.Ch2
        @test sort(cont_ab.analysis_info.repaired_channels) == [:Ch1, :Ch2]

        # 2. All neighbors bad -> warns, skips channel cleanly without corrupting with NaN, and omits from analysis_info
        cont_allbad = EegFun.ContinuousData("rep_skip", copy(df_rep), lay_rep, 100, EegFun.AnalysisInfo())
        @test_logs (:warn,) match_mode = :any EegFun.repair_channels!(
            cont_allbad,
            [:Ch1, :Ch2, :Ch3, :Ch4];
            method = :neighbor_interpolation,
        )
        @test cont_allbad.data.Ch1 == df_rep.Ch1
        @test isempty(cont_allbad.analysis_info.repaired_channels)

        # 3. ErpData support for both :neighbor_interpolation and :spherical_spline
        erp_dat = EegFun.create_test_erp_data(n_channels = 4, seed = 10)
        erp_dat.layout.data.x3 = [0.0, 0.1, 0.0, -0.1]
        erp_dat.layout.data.y3 = [0.1, 0.0, -0.1, 0.0]
        erp_dat.layout.data.z3 = [0.0, 0.0, 0.0, 0.0]
        EegFun.get_neighbors_xyz!(erp_dat.layout, 0.5)
        erp_dat.data[!, :Ch2] .= 999.0

        erp_nb = EegFun.repair_channels(erp_dat, [:Ch2]; method = :neighbor_interpolation)
        @test erp_nb isa EegFun.ErpData
        @test all(abs.(erp_nb.data.Ch2) .< 100.0)
        @test erp_nb.analysis_info.repaired_channels == [:Ch2]

        erp_sp = EegFun.repair_channels(erp_dat, [:Ch2]; method = :spherical_spline)
        @test erp_sp isa EegFun.ErpData
        @test all(abs.(erp_sp.data.Ch2) .< 100.0)
        @test erp_sp.analysis_info.repaired_channels == [:Ch2]
    end

end
