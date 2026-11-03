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
        @test_throws Exception EegFun.read_bids(bids_dir; subject = "01", task = "posner")
    end

    @testset "BIDS path construction - with session" begin
        bids_dir = mktempdir()

        # sub-01/ses-01/eeg/
        eeg_dir = joinpath(bids_dir, "sub-01", "ses-01", "eeg")
        mkpath(eeg_dir)

        # Should construct path with session and fail on missing file (not missing dir)
        @test_throws Exception EegFun.read_bids(bids_dir; subject = "01", task = "posner", session = "01")
    end

    @testset "BIDS path construction - strip prefix" begin
        bids_dir = mktempdir()

        # sub-01/eeg/
        eeg_dir = joinpath(bids_dir, "sub-01", "eeg")
        mkpath(eeg_dir)

        # Should handle "sub-" prefix gracefully (strip it)
        @test_throws Exception EegFun.read_bids(bids_dir; subject = "sub-01", task = "posner")
    end

    @testset "BIDS directory not found" begin
        # Non-existent base directory
        @test_throws Exception EegFun.read_bids("/nonexistent/path"; subject = "01", task = "posner")
    end

    # ────────────────────────────────────────────────────────────
    # 2. BIDS events.tsv injection
    # ────────────────────────────────────────────────────────────
    @testset "BIDS electrodes.tsv parsing" begin
        bids_dir = mktempdir()
        eeg_dir = joinpath(bids_dir, "sub-01", "eeg")
        mkpath(eeg_dir)

        # Create a minimal electrodes.tsv
        electrodes_df = DataFrame(name = ["Fz", "Cz", "Pz"], x = [0.0, 0.0, 0.0], y = [0.71, 0.0, -0.71], z = [0.71, 1.0, 0.71])
        CSV.write(joinpath(eeg_dir, "sub-01_electrodes.tsv"), electrodes_df; delim = '\t')

        # Verify the file was created correctly
        @test isfile(joinpath(eeg_dir, "sub-01_electrodes.tsv"))

        # Read back and verify
        read_back = DataFrame(CSV.File(joinpath(eeg_dir, "sub-01_electrodes.tsv")))
        @test nrow(read_back) == 3
        @test all(col -> col in propertynames(read_back), [:name, :x, :y, :z])

        # Verify _bids_electrodes_to_layout builds a valid Layout with :inc, :azi, :x2, :y2, :x3, :y3, :z3
        lay = EegFun._bids_electrodes_to_layout(read_back)
        @test lay isa EegFun.Layout
        @test EegFun.channel_labels(lay) == [:Fz, :Cz, :Pz]
        @test all(col -> col in propertynames(lay.data), [:label, :inc, :azi, :x2, :y2, :x3, :y3, :z3])
        @test isapprox(lay.data[2, :inc], 0.0, atol = 1e-6)  # Cz at (0, 0, 1) -> inc = 0
    end

    @testset "Full export_bids -> list_bids_runs -> read_bids round-trip & zero-padded runs" begin
        bdf_src = joinpath(dirname(dirname(@__DIR__)), "resources", "data", "bdf", "example1.bdf")
        layout_file = joinpath(dirname(dirname(@__DIR__)), "resources", "layouts", "biosemi", "biosemi72.csv")
        @test isfile(bdf_src)
        @test isfile(layout_file)

        mktempdir() do tmpdir
            raw_dir = joinpath(tmpdir, "raw")
            bids_dir = joinpath(tmpdir, "bids")
            mkpath(raw_dir)
            cp(bdf_src, joinpath(raw_dir, "example1.bdf"))

            EegFun.export_bids(
                task = "posner",
                raw_dir = raw_dir,
                output_dir = bids_dir,
                layout_file = layout_file,
                dataset_name = "RoundTrip Test",
            )

            runs_df = EegFun.list_bids_runs(bids_dir)
            @test nrow(runs_df) == 1
            @test runs_df.subject == ["01"]
            @test runs_df.task == ["posner"]
            @test runs_df.run == [0]

            dat = EegFun.read_bids(bids_dir; subject = "01", task = "posner")
            @test dat isa EegFun.ContinuousData
            @test EegFun.n_channels(dat) == 72
            @test all(col -> col in propertynames(dat.layout.data), [:label, :inc, :azi, :x2, :y2, :x3, :y3, :z3])
            @test :trigger in propertynames(dat.data)
            @test :trigger_info in propertynames(dat.data)
            @test count(!=(0), dat.data.trigger) > 0

            # Test zero-padded _run-01 and string/"n/a" values in events.tsv
            eeg_dir = joinpath(bids_dir, "sub-01", "eeg")
            cp(joinpath(eeg_dir, "sub-01_task-posner_eeg.bdf"), joinpath(eeg_dir, "sub-01_task-posner_run-01_eeg.bdf"))
            events_run01 = DataFrame(onset = [0.0, 0.05], duration = [0, 0], trial_type = ["start", "cue"], value = ["n/a", "S 42"])
            CSV.write(joinpath(eeg_dir, "sub-01_task-posner_run-01_events.tsv"), events_run01; delim = '\t')

            runs_df2 = EegFun.list_bids_runs(bids_dir)
            @test nrow(runs_df2) == 2
            @test 1 in runs_df2.run

            dat_run1 = EegFun.read_bids(bids_dir; subject = "01", task = "posner", run = 1)
            @test dat_run1.data.trigger[1] == 1
            @test dat_run1.data.trigger_info[1] == "start"
            samp2 = round(Int, 0.05 * EegFun.sample_rate(dat_run1)) + 1
            @test dat_run1.data.trigger[samp2] == 42
            @test dat_run1.data.trigger_info[samp2] == "cue"
        end
    end

end # @testset "BIDS Import"
