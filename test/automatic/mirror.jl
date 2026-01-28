using DataFrames

@testset "Data Mirroring" begin

    @testset "Mirror EpochData - :pre" begin

        # Create simple test epoch
        time = -0.2:0.1:0.2
        n_samples = length(time)

        epoch1 = DataFrame(time = collect(time), Cz = collect(1.0:n_samples), Pz = collect(n_samples:-1.0:1))

        epoch2 = DataFrame(time = collect(time), Cz = collect(1.0:n_samples) .* 2, Pz = collect(n_samples:-1.0:1) .* 2)

        epochs = EegFun.EpochData(
            "test_data",
            1,
            "condition_1",
            [epoch1, epoch2],
            EegFun.Layout(DataFrame(), nothing, nothing),
            100,
            EegFun.AnalysisInfo(),
        )

        original_length = nrow(epochs.data[1])

        # Mirror :pre
        EegFun.mirror!(epochs, :pre)

        # Check length increased
        @test nrow(epochs.data[1]) > original_length

        # Check time vector continuity
        time_diffs = diff(epochs.data[1].time)
        @test all(time_diffs .≈ time_diffs[1])  # Uniform spacing

        # Unmirror
        EegFun.unmirror!(epochs, :pre)

        # Check restored to original
        @test nrow(epochs.data[1]) == original_length
        @test epochs.data[1].time ≈ collect(time)
        @test epochs.data[1].Cz ≈ collect(1.0:n_samples)
    end


    @testset "Mirror EpochData - :post" begin
        time = -0.2:0.1:0.2
        n_samples = length(time)

        epoch1 = DataFrame(time = collect(time), Cz = collect(1.0:n_samples), Pz = collect(n_samples:-1.0:1))

        epochs = EegFun.EpochData(
            "test_data",
            1,
            "condition_1",
            [epoch1],
            EegFun.Layout(DataFrame(), nothing, nothing),
            100,
            EegFun.AnalysisInfo(),
        )

        original_length = nrow(epochs.data[1])
        original_data = copy(epochs.data[1])

        # Mirror :post
        EegFun.mirror!(epochs, :post)

        # Check length increased
        @test nrow(epochs.data[1]) > original_length

        # Check time vector continuity
        time_diffs = diff(epochs.data[1].time)
        @test all(time_diffs .≈ time_diffs[1])

        # Unmirror
        EegFun.unmirror!(epochs, :post)

        # Check restored
        @test nrow(epochs.data[1]) == original_length
        @test epochs.data[1].time ≈ original_data.time
        @test epochs.data[1].Cz ≈ original_data.Cz
    end


    @testset "Mirror EpochData - :both" begin
        time = -0.2:0.1:0.2
        n_samples = length(time)

        epoch1 = DataFrame(
            time = collect(time),
            Cz = collect(1.0:n_samples),
            Pz = collect(n_samples:-1.0:1),
            trial = fill(1, n_samples),
            condition = fill(1, n_samples),
        )

        epochs = EegFun.EpochData(
            "test_data",
            1,
            "condition_1",
            [epoch1],
            EegFun.Layout(DataFrame(), nothing, nothing),
            100,
            EegFun.AnalysisInfo(),
        )

        original_length = nrow(epochs.data[1])
        original_time = copy(epochs.data[1].time)
        original_cz = copy(epochs.data[1].Cz)

        # Mirror :both
        EegFun.mirror!(epochs, :both)

        # Check length (should be roughly 3× original)
        mirrored_length = nrow(epochs.data[1])
        @test mirrored_length > 2 * original_length

        # Check time vector continuity
        time_diffs = diff(epochs.data[1].time)
        @test all(time_diffs .≈ time_diffs[1])

        # Unmirror
        EegFun.unmirror!(epochs, :both)

        # Check fully restored
        @test nrow(epochs.data[1]) == original_length
        @test epochs.data[1].time ≈ original_time
        @test epochs.data[1].Cz ≈ original_cz
    end


    @testset "Mirror EpochData - Non-mutating" begin
        time = -0.2:0.1:0.2
        n_samples = length(time)

        epoch1 = DataFrame(time = collect(time), Cz = collect(1.0:n_samples))

        epochs_original = EegFun.EpochData(
            "test_data",
            1,
            "condition_1",
            [epoch1],
            EegFun.Layout(DataFrame(), nothing, nothing),
            100,
            EegFun.AnalysisInfo(),
        )
        original_length = nrow(epochs_original.data[1])

        # Non-mutating mirror
        epochs_mirrored = EegFun.mirror(epochs_original, :both)

        # Check original unchanged
        @test nrow(epochs_original.data[1]) == original_length

        # Check mirrored is different
        @test nrow(epochs_mirrored.data[1]) > original_length

        # Non-mutating unmirror
        epochs_unmirrored = EegFun.unmirror(epochs_mirrored, :both)

        # Check unmirrored matches original
        @test nrow(epochs_unmirrored.data[1]) == original_length
        @test epochs_unmirrored.data[1].time ≈ epochs_original.data[1].time
        @test epochs_unmirrored.data[1].Cz ≈ epochs_original.data[1].Cz
    end


    @testset "Mirror ErpData - :both" begin
        time = -0.2:0.1:0.2
        n_samples = length(time)

        erp_df = DataFrame(time = collect(time), Cz = collect(1.0:n_samples), Pz = collect(n_samples:-1.0:1))

        erp = EegFun.ErpData(
            "test_data",
            1,
            "condition_1",
            erp_df,
            EegFun.Layout(DataFrame(), nothing, nothing),
            100,
            EegFun.AnalysisInfo(),
            10,
        )

        original_length = nrow(erp.data)
        original_time = copy(erp.data.time)
        original_cz = copy(erp.data.Cz)

        # Mirror
        EegFun.mirror!(erp, :both)

        # Check length increased
        @test nrow(erp.data) > 2 * original_length

        # Check time continuity
        time_diffs = diff(erp.data.time)
        @test all(time_diffs .≈ time_diffs[1])

        # Unmirror
        EegFun.unmirror!(erp, :both)

        # Check restored
        @test nrow(erp.data) == original_length
        @test erp.data.time ≈ original_time
        @test erp.data.Cz ≈ original_cz
    end


    @testset "Mirror ErpData - Non-mutating" begin
        time = -0.2:0.1:0.2
        n_samples = length(time)

        erp_df = DataFrame(time = collect(time), Cz = collect(1.0:n_samples))

        erp_original = EegFun.ErpData(
            "test_data",
            1,
            "condition_1",
            erp_df,
            EegFun.Layout(DataFrame(), nothing, nothing),
            100,
            EegFun.AnalysisInfo(),
            10,
        )
        original_length = nrow(erp_original.data)

        # Non-mutating mirror
        erp_mirrored = EegFun.mirror(erp_original, :both)

        # Check original unchanged
        @test nrow(erp_original.data) == original_length

        # Check mirrored is different
        @test nrow(erp_mirrored.data) > original_length

        # Non-mutating unmirror
        erp_unmirrored = EegFun.unmirror(erp_mirrored, :both)

        # Check matches original
        @test nrow(erp_unmirrored.data) == original_length
        @test erp_unmirrored.data.time ≈ erp_original.data.time
        @test erp_unmirrored.data.Cz ≈ erp_original.data.Cz
    end


    @testset "Invalid side parameter" begin
        time = -0.2:0.1:0.2
        n_samples = length(time)

        epoch1 = DataFrame(time = collect(time), Cz = collect(1.0:n_samples))

        epochs = EegFun.EpochData(
            "test_data",
            1,
            "condition_1",
            [epoch1],
            EegFun.Layout(DataFrame(), nothing, nothing),
            100,
            EegFun.AnalysisInfo(),
        )

        # Test invalid side
        @test_throws Exception mirror!(epochs, :invalid)
        @test_throws Exception unmirror!(epochs, :middle)
    end


    @testset "Multiple epochs" begin
        time = -0.2:0.1:0.2
        n_samples = length(time)

        # Create 3 different epochs
        epoch1 = DataFrame(time = collect(time), Cz = collect(1.0:n_samples))

        epoch2 = DataFrame(time = collect(time), Cz = collect(1.0:n_samples) .* 2)

        epoch3 = DataFrame(time = collect(time), Cz = collect(n_samples:-1.0:1))

        epochs = EegFun.EpochData(
            "test_data",
            1,
            "condition_1",
            [epoch1, epoch2, epoch3],
            EegFun.Layout(DataFrame(), nothing, nothing),
            100,
            EegFun.AnalysisInfo(),
        )

        original_cz1 = copy(epochs.data[1].Cz)
        original_cz2 = copy(epochs.data[2].Cz)
        original_cz3 = copy(epochs.data[3].Cz)

        # Mirror and unmirror
        EegFun.mirror!(epochs, :both)
        EegFun.unmirror!(epochs, :both)

        # Check all epochs restored correctly
        @test epochs.data[1].Cz ≈ original_cz1
        @test epochs.data[2].Cz ≈ original_cz2
        @test epochs.data[3].Cz ≈ original_cz3
    end


    @testset "Mirroring preserves metadata columns" begin
        time = -0.2:0.1:0.2
        n_samples = length(time)

        epoch1 = DataFrame(
            time = collect(time),
            Cz = collect(1.0:n_samples),
            Pz = collect(n_samples:-1.0:1),
            trial = fill(1, n_samples),
            condition = fill(2, n_samples),
            response = fill("left", n_samples),
        )

        epochs = EegFun.EpochData(
            "test_data",
            2,
            "condition_2",
            [epoch1],
            EegFun.Layout(DataFrame(), nothing, nothing),
            100,
            EegFun.AnalysisInfo(),
        )

        # Mirror
        EegFun.mirror!(epochs, :both)

        # Check metadata preserved
        @test all(epochs.data[1].trial .== 1)
        @test epochs.condition == 2
        @test all(epochs.data[1].response .== "left")

        # Unmirror
        EegFun.unmirror!(epochs, :both)

        # Check metadata still there
        @test all(epochs.data[1].trial .== 1)
        @test all(epochs.data[1].condition .== 2)
        @test all(epochs.data[1].response .== "left")
    end

    @testset "Mirroring pattern validation" begin
        # Test specific mirroring patterns with known data
        time = [0.0, 0.1, 0.2, 0.3, 0.4]  # 5 samples: [0, 1, 2, 3, 4]
        data = [1.0, 2.0, 3.0, 4.0, 5.0]  # Simple ascending pattern

        epoch1 = DataFrame(time = time, Cz = data)
        epochs = EegFun.EpochData(
            "test_data",
            1,
            "condition_1",
            [epoch1],
            EegFun.Layout(DataFrame(), nothing, nothing),
            10,
            EegFun.AnalysisInfo(),
        )

        # Test :pre mirroring
        EegFun.mirror!(epochs, :pre)
        expected_pre = [5.0, 4.0, 3.0, 2.0, 1.0, 2.0, 3.0, 4.0, 5.0]  # [5,4,3,2] + [1,2,3,4,5]
        @test epochs.data[1].Cz ≈ expected_pre

        # Reset and test :post mirroring  
        EegFun.unmirror!(epochs, :pre)
        EegFun.mirror!(epochs, :post)
        expected_post = [1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0]  # [1,2,3,4,5] + [4,3,2,1]
        @test epochs.data[1].Cz ≈ expected_post

        # Reset and test :both mirroring
        EegFun.unmirror!(epochs, :post)
        EegFun.mirror!(epochs, :both)
        expected_both = [5.0, 4.0, 3.0, 2.0, 1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0]
        @test epochs.data[1].Cz ≈ expected_both

        # Test no duplication at boundaries (there should be no duplications)
        diffs = diff(epochs.data[1].Cz)
        @test sum(diffs .== 0) == 0  # No duplications

        # Test symmetry - the mirrored sections should be symmetric around the original data
        n_orig = 5
        n_mirrored = length(epochs.data[1].Cz)
        mid_start = n_orig  # Start of original data in mirrored array (1-indexed)
        mid_end = mid_start + n_orig - 1  # End of original data

        # Pre-mirror should be reverse of the original data (excluding first point)
        pre_mirror = epochs.data[1].Cz[1:(n_orig-1)]  # [5, 4, 3, 2]
        original_data = epochs.data[1].Cz[mid_start:mid_end]  # [1, 2, 3, 4, 5]
        @test pre_mirror ≈ reverse(original_data[2:end])  # reverse([2, 3, 4, 5]) = [5, 4, 3, 2]

        # Post-mirror should be reverse of the original data (excluding last point)
        post_mirror = epochs.data[1].Cz[(mid_end+1):end]  # [4, 3, 2, 1]
        @test post_mirror ≈ reverse(original_data[1:(end-1)])  # reverse([1, 2, 3, 4]) = [4, 3, 2, 1]
    end


    @testset "Roundtrip consistency" begin
        # Test that mirror → unmirror → mirror → unmirror works
        time = -0.2:0.05:0.2
        n_samples = length(time)

        epoch1 = DataFrame(time = collect(time), Cz = randn(n_samples), Pz = randn(n_samples))

        epochs = EegFun.EpochData(
            "test_data",
            1,
            "condition_1",
            [epoch1],
            EegFun.Layout(DataFrame(), nothing, nothing),
            100,
            EegFun.AnalysisInfo(),
        )

        original_cz = copy(epochs.data[1].Cz)
        original_pz = copy(epochs.data[1].Pz)
        original_time = copy(epochs.data[1].time)

        # First roundtrip
        EegFun.mirror!(epochs, :both)
        EegFun.unmirror!(epochs, :both)

        @test epochs.data[1].Cz ≈ original_cz
        @test epochs.data[1].Pz ≈ original_pz
        @test epochs.data[1].time ≈ original_time

        # Second roundtrip
        EegFun.mirror!(epochs, :both)
        EegFun.unmirror!(epochs, :both)

        @test epochs.data[1].Cz ≈ original_cz
        @test epochs.data[1].Pz ≈ original_pz
        @test epochs.data[1].time ≈ original_time
    end
end
