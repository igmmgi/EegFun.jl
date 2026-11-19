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
            orig_lay = EegFun.read_layout(layout_file)
            EegFun.polar_to_cartesian_xy!(orig_lay)
            @test maximum(abs.(abs.(dat.layout.data.inc) .- abs.(orig_lay.data.inc))) < 1e-3
            @test maximum(abs.(dat.layout.data.x2 .- orig_lay.data.x2)) < 1e-4
            @test maximum(abs.(dat.layout.data.y2 .- orig_lay.data.y2)) < 1e-4
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

            # Categorical-only trial_type without numeric value column assigns distinct integer codes
            cp(joinpath(eeg_dir, "sub-01_task-posner_eeg.bdf"), joinpath(eeg_dir, "sub-01_task-posner_acq-high_run-02_eeg.bdf"))
            events_run02 = DataFrame(onset = [0.0, 0.05, 0.10], duration = [0, 0, 0], trial_type = ["target", "standard", "target"])
            CSV.write(joinpath(eeg_dir, "sub-01_task-posner_acq-high_run-02_events.tsv"), events_run02; delim = '\t')
            dat_run2 = EegFun.read_bids(bids_dir; subject = "01", task = "posner", run = 2)
            s_t1 = 1
            s_std = round(Int, 0.05 * EegFun.sample_rate(dat_run2)) + 1
            s_t2 = round(Int, 0.10 * EegFun.sample_rate(dat_run2)) + 1
            @test dat_run2.data.trigger[s_std] == 1
            @test dat_run2.data.trigger[s_t1] == 2
            @test dat_run2.data.trigger[s_t2] == 2
        end

        # BrainVision multi-file (.vhdr + .eeg + .vmrk) export_bids -> read_bids round-trip
        bv_dir = joinpath(dirname(dirname(@__DIR__)), "resources", "data", "brainvision")
        if isfile(joinpath(bv_dir, "example1.vhdr"))
            mktempdir() do tmpbv
                bids_bv = joinpath(tmpbv, "bids_bv")
                EegFun.export_bids(task = "oddball", raw_dir = bv_dir, raw_pattern = "\\.vhdr\$", output_dir = bids_bv)
                @test isfile(joinpath(bids_bv, "sub-01", "eeg", "sub-01_task-oddball_eeg.vhdr"))
                @test isfile(joinpath(bids_bv, "sub-01", "eeg", "sub-01_task-oddball_eeg.eeg"))
                @test isfile(joinpath(bids_bv, "sub-01", "eeg", "sub-01_task-oddball_eeg.vmrk"))
                dat_bv = EegFun.read_bids(bids_bv; subject = "01", task = "oddball")
                @test dat_bv isa EegFun.ContinuousData
                @test EegFun.n_samples(dat_bv) > 0

                # Injecting fiducials (NAS/LPA/RPA) into electrodes.tsv does not pollute dat.layout or break repair_channels
                bv_eeg_dir = joinpath(bids_bv, "sub-01", "eeg")
                bv_chs = String.(EegFun.channel_labels(dat_bv))
                elec_names = vcat(["NAS", "LPA", "RPA"], bv_chs)
                n_el = length(elec_names)
                angles = range(0.1, 2π - 0.1, length = n_el)
                elec_df = DataFrame(name = elec_names, x = 0.5 .* cos.(angles), y = 0.5 .* sin.(angles), z = fill(0.8, n_el))
                CSV.write(joinpath(bv_eeg_dir, "sub-01_electrodes.tsv"), elec_df; delim = '\t')

                dat_bv_fid = EegFun.read_bids(bids_bv; subject = "01", task = "oddball")
                @test :NAS ∉ dat_bv_fid.layout.data.label
                @test dat_bv_fid.layout.data.label == Symbol.(bv_chs)

                # Even if a layout manually contains phantom fiducials, spherical & neighbor repair skip non-DataFrame rows
                dat_bv_fid.layout = EegFun._bids_electrodes_to_layout(elec_df)
                EegFun.get_neighbors_xyz!(dat_bv_fid.layout, 2.0)
                rep_sph = EegFun.repair_channels(dat_bv_fid, [Symbol(bv_chs[1])]; method = :spherical_spline)
                @test rep_sph isa EegFun.ContinuousData
                rep_nb = EegFun.repair_channels(dat_bv_fid, [Symbol(bv_chs[1])]; method = :neighbor_interpolation)
                @test rep_nb isa EegFun.ContinuousData
            end
        end
    end

end # @testset "BIDS Import"
