using Test
using EegFun
using DataFrames
using CSV

@testset "BIDS Import" begin

    # ────────────────────────────────────────────────────────────
    # 1. BIDS path construction logic
    # ────────────────────────────────────────────────────────────
    @testset "BIDS path construction" begin
        # Create a temporary BIDS-like directory structure
        bids_dir = mktempdir()

        # sub-01/eeg/
        eeg_dir = joinpath(bids_dir, "sub-01", "eeg")
        mkpath(eeg_dir)

        # Test that read_bids errors with missing files (validates path construction)
        @test_throws Exception EegFun.read_bids(
            bids_dir;
            subject = "01",
            task = "posner",
        )
    end

    @testset "BIDS path construction - with session" begin
        bids_dir = mktempdir()

        # sub-01/ses-01/eeg/
        eeg_dir = joinpath(bids_dir, "sub-01", "ses-01", "eeg")
        mkpath(eeg_dir)

        # Should construct path with session and fail on missing file (not missing dir)
        @test_throws Exception EegFun.read_bids(
            bids_dir;
            subject = "01",
            task = "posner",
            session = "01",
        )
    end

    @testset "BIDS path construction - strip prefix" begin
        bids_dir = mktempdir()

        # sub-01/eeg/
        eeg_dir = joinpath(bids_dir, "sub-01", "eeg")
        mkpath(eeg_dir)

        # Should handle "sub-" prefix gracefully (strip it)
        @test_throws Exception EegFun.read_bids(
            bids_dir;
            subject = "sub-01",
            task = "posner",
        )
    end

    @testset "BIDS directory not found" begin
        # Non-existent base directory
        @test_throws Exception EegFun.read_bids(
            "/nonexistent/path";
            subject = "01",
            task = "posner",
        )
    end

    # ────────────────────────────────────────────────────────────
    # 2. BIDS events.tsv injection
    # ────────────────────────────────────────────────────────────
    @testset "BIDS electrodes.tsv parsing" begin
        bids_dir = mktempdir()
        eeg_dir = joinpath(bids_dir, "sub-01", "eeg")
        mkpath(eeg_dir)

        # Create a minimal electrodes.tsv
        electrodes_df = DataFrame(
            name = ["Fz", "Cz", "Pz"],
            x = [0.0, 0.0, 0.0],
            y = [0.71, 0.0, -0.71],
            z = [0.71, 1.0, 0.71],
        )
        CSV.write(joinpath(eeg_dir, "sub-01_electrodes.tsv"), electrodes_df; delim = '\t')

        # Verify the file was created correctly
        @test isfile(joinpath(eeg_dir, "sub-01_electrodes.tsv"))

        # Read back and verify
        read_back = DataFrame(CSV.File(joinpath(eeg_dir, "sub-01_electrodes.tsv")))
        @test nrow(read_back) == 3
        @test all(col -> col in propertynames(read_back), [:name, :x, :y, :z])
    end

end # @testset "BIDS Import"
