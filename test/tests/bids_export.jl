using DataFrames

@testset "BIDS Export" begin

    # --- Helpers ---
    function _make_test_continuous_data(filename; n_samples = 100, sample_rate = 256)
        labels = [:Fp1, :Fp2, :Cz, :Pz, :O1, :O2]
        n_channels = length(labels)
        time = collect(range(0, step = 1 / sample_rate, length = n_samples))
        data = randn(n_samples, n_channels)

        trigger = zeros(Int, n_samples)
        trigger[10] = 1
        trigger[30] = 2
        trigger[50] = 3
        trigger[70] = 4

        df = DataFrame(:time => time, :sample => collect(1:n_samples), :trigger => trigger)
        for (i, label) in enumerate(labels)
            df[!, label] = data[:, i]
        end

        layout_df = DataFrame(label = labels, inc = [20.0, 20.0, 0.0, 45.0, 80.0, 80.0], azi = [-45.0, 45.0, 0.0, 180.0, -60.0, 60.0])
        layout = EegFun.Layout(layout_df, nothing, nothing, nothing)

        return EegFun.ContinuousData(
            filename,
            df,
            layout,
            sample_rate,
            EegFun.AnalysisInfo(reference = :avg, hp_filter = 0.1, lp_filter = 0.0),
        )
    end


    @testset "Participant ID mapping" begin
        raw_files = ["example1.bdf", "example2.bdf", "subject10.bdf"]

        # Auto-extract
        pmap = EegFun._bids_build_participant_map(raw_files, nothing)
        @test pmap["example1"] == "sub-01"
        @test pmap["example2"] == "sub-02"
        @test pmap["subject10"] == "sub-10"

        # User-provided
        user_map = Dict("example1" => "sub-P01", "example2" => "sub-P02", "subject10" => "sub-P10")
        pmap2 = EegFun._bids_build_participant_map(raw_files, user_map)
        @test pmap2["example1"] == "sub-P01"
        @test pmap2["example2"] == "sub-P02"
    end


    @testset "Channel classification" begin
        @test EegFun._bids_classify_channel(:Fp1) == "EEG"
        @test EegFun._bids_classify_channel(:Cz) == "EEG"
        @test EegFun._bids_classify_channel(:vEOG) == "EOG"
        @test EegFun._bids_classify_channel(:hEOG) == "EOG"
        @test EegFun._bids_classify_channel(:IO1) == "EOG"
        @test EegFun._bids_classify_channel(:IO2) == "EOG"
        @test EegFun._bids_classify_channel(:M1) == "EEG"
        @test EegFun._bids_classify_channel(:EMG1) == "EMG"
        @test EegFun._bids_classify_channel(:ECG) == "ECG"
    end


    @testset "JSON serialiser" begin
        buf = IOBuffer()
        data = OrderedCollections.OrderedDict(
            "Name" => "Test",
            "BIDSVersion" => "1.9.0",
            "Values" => [1, 2, 3],
            "Nested" => OrderedCollections.OrderedDict("Key" => "Value"),
        )
        EegFun._json_print(buf, data, 0)
        json_str = String(take!(buf))
        @test contains(json_str, "\"Name\": \"Test\"")
        @test contains(json_str, "\"BIDSVersion\": \"1.9.0\"")
        @test contains(json_str, "[1, 2, 3]")
        @test contains(json_str, "\"Key\": \"Value\"")
    end


    @testset "Sidecar file generation" begin
        mktempdir() do tmpdir
            dat = _make_test_continuous_data("test_file")
            prefix = "sub-01_task-test"

            # EEG sidecar
            EegFun._bids_write_eeg_sidecar(tmpdir, prefix, dat, "test", 50, Dict{String,String}())
            eeg_json = joinpath(tmpdir, "$(prefix)_eeg.json")
            @test isfile(eeg_json)
            content = read(eeg_json, String)
            @test contains(content, "\"SamplingFrequency\": 256")
            @test contains(content, "\"TaskName\": \"test\"")
            @test contains(content, "\"PowerLineFrequency\": 50")
            @test contains(content, "\"EEGChannelCount\": 6")
            @test contains(content, "\"RecordingType\": \"continuous\"")

            # Channels TSV
            EegFun._bids_write_channels_tsv(tmpdir, prefix, dat)
            channels_tsv = joinpath(tmpdir, "$(prefix)_channels.tsv")
            @test isfile(channels_tsv)
            lines = readlines(channels_tsv)
            @test lines[1] == "name\ttype\tunits\tstatus"
            @test length(lines) == 7  # header + 6 channels

            # Events TSV
            EegFun._bids_write_events_tsv(tmpdir, prefix, dat)
            events_tsv = joinpath(tmpdir, "$(prefix)_events.tsv")
            @test isfile(events_tsv)
            events_lines = readlines(events_tsv)
            @test events_lines[1] == "onset\tduration\ttrial_type\tvalue"
            @test length(events_lines) == 5  # header + 4 triggers

            # Electrodes TSV
            EegFun._bids_write_electrodes_tsv(tmpdir, prefix, dat.layout)
            electrodes_tsv = joinpath(tmpdir, "$(prefix)_electrodes.tsv")
            @test isfile(electrodes_tsv)
            elec_lines = readlines(electrodes_tsv)
            @test elec_lines[1] == "name\tx\ty\tz"
            @test length(elec_lines) == 7  # header + 6 electrodes

            # Coordsystem JSON
            EegFun._bids_write_coordsystem_json(tmpdir, prefix)
            coord_json = joinpath(tmpdir, "$(prefix)_coordsystem.json")
            @test isfile(coord_json)
            @test contains(read(coord_json, String), "EEGCoordinateSystem")
        end
    end


    @testset "Dataset description" begin
        mktempdir() do tmpdir
            EegFun._bids_write_dataset_description(tmpdir, "My Study", ["Author A", "Author B"], "CC0")
            desc_path = joinpath(tmpdir, "dataset_description.json")
            @test isfile(desc_path)
            content = read(desc_path, String)
            @test contains(content, "\"Name\": \"My Study\"")
            @test contains(content, "\"BIDSVersion\": \"1.9.0\"")
            @test contains(content, "\"Author A\"")
        end
    end


    @testset "Participants TSV" begin
        mktempdir() do tmpdir
            pmap = Dict("file1" => "sub-01", "file2" => "sub-02")
            EegFun._bids_write_participants_tsv(tmpdir, pmap)

            tsv_path = joinpath(tmpdir, "participants.tsv")
            @test isfile(tsv_path)
            lines = readlines(tsv_path)
            @test lines[1] == "participant_id"
            @test "sub-01" in lines
            @test "sub-02" in lines

            json_path = joinpath(tmpdir, "participants.json")
            @test isfile(json_path)
        end
    end


    @testset "Derivatives renaming" begin
        mktempdir() do tmpdir
            # Create mock derivatives directory with JLD2 files
            deriv_dir = joinpath(tmpdir, "preprocessed")
            mkpath(deriv_dir)

            # Create dummy JLD2 files (just empty files for naming tests)
            for suffix in ["_epochs_uncorrected", "_erps_final", "_ica", "_artifact_info"]
                touch(joinpath(deriv_dir, "example1$(suffix).jld2"))
            end

            output_dir = joinpath(tmpdir, "bids_out")
            mkpath(output_dir)

            raw_files = ["example1.bdf"]
            pmap = Dict("example1" => "sub-01")

            EegFun._bids_copy_derivatives(output_dir, deriv_dir, raw_files, pmap, "posner", false)

            sub_dir = joinpath(output_dir, "derivatives", "EegFun", "sub-01", "eeg")
            @test isdir(sub_dir)

            # Check renamed files
            files = readdir(sub_dir)
            @test "sub-01_task-posner_desc-epochsUncorrected_eeg.jld2" in files
            @test "sub-01_task-posner_desc-erpsGood_eeg.jld2" in files
            @test "sub-01_task-posner_desc-ica_eeg.jld2" in files
            @test "sub-01_task-posner_desc-artifactInfo_eeg.jld2" in files

            # Check derivatives dataset_description.json
            deriv_desc = joinpath(output_dir, "derivatives", "EegFun", "dataset_description.json")
            @test isfile(deriv_desc)
            @test contains(read(deriv_desc, String), "derivative")
        end
    end


    @testset "Validation" begin
        # Missing task
        @test_throws Exception EegFun.export_bids(task = "")

        # Non-existent raw_dir
        @test_throws Exception EegFun.export_bids(task = "test", raw_dir = "/nonexistent_dir_xyz")

        # Overwrite protection
        mktempdir() do tmpdir
            output = joinpath(tmpdir, "exists")
            mkpath(output)
            @test_throws Exception EegFun.export_bids(task = "test", raw_dir = tmpdir, output_dir = output, overwrite = false)
        end
    end


    @testset "README generation" begin
        mktempdir() do tmpdir
            EegFun._bids_write_readme(tmpdir, "My Dataset")
            readme = joinpath(tmpdir, "README")
            @test isfile(readme)
            content = read(readme, String)
            @test contains(content, "My Dataset")
            @test contains(content, "EegFun.jl")
            @test contains(content, "BIDS")

            # Empty name gets a default
            EegFun._bids_write_readme(tmpdir, "")
            content2 = read(readme, String)
            @test contains(content2, "Untitled Dataset")
        end
    end


    @testset "Events JSON (no onset/duration redefinition)" begin
        mktempdir() do tmpdir
            prefix = "sub-01_task-test"
            EegFun._bids_write_events_json(tmpdir, prefix)
            content = read(joinpath(tmpdir, "$(prefix)_events.json"), String)
            @test contains(content, "value")
            @test !contains(content, "\"onset\"")      # must NOT redefine onset
            @test !contains(content, "\"duration\"")    # must NOT redefine duration
        end
    end


    @testset "Events TSV duration is 0" begin
        mktempdir() do tmpdir
            dat = _make_test_continuous_data("test_file")
            prefix = "sub-01_task-test"
            EegFun._bids_write_events_tsv(tmpdir, prefix, dat)
            lines = readlines(joinpath(tmpdir, "$(prefix)_events.tsv"))
            # Check a data line has 0 for duration (second column)
            fields = split(lines[2], '\t')
            @test fields[2] == "0"
        end
    end


    @testset "JSON Bool and Symbol serialisation" begin
        buf = IOBuffer()
        # Bool must serialize as true/false, not 1/0
        EegFun._json_print(buf, true, 0)
        @test String(take!(buf)) == "true"

        EegFun._json_print(buf, false, 0)
        @test String(take!(buf)) == "false"

        # Symbol gets quoted
        EegFun._json_print(buf, :test_symbol, 0)
        @test String(take!(buf)) == "\"test_symbol\""
    end


    @testset "Recommended metadata in EEG sidecar" begin
        mktempdir() do tmpdir
            dat = _make_test_continuous_data("test_file")
            prefix = "sub-01_task-test"
            recommended =
                Dict{String,String}("Manufacturer" => "BioSemi", "InstitutionName" => "Test University", "TaskDescription" => "A test task")
            EegFun._bids_write_eeg_sidecar(tmpdir, prefix, dat, "test", 50, recommended)
            content = read(joinpath(tmpdir, "$(prefix)_eeg.json"), String)
            @test contains(content, "\"Manufacturer\": \"BioSemi\"")
            @test contains(content, "\"InstitutionName\": \"Test University\"")
            @test contains(content, "\"TaskDescription\": \"A test task\"")
        end
    end


    @testset "get_files recursive" begin
        mktempdir() do tmpdir
            # Create a BIDS-like structure
            for sub in ["sub-01", "sub-02"]
                dir = joinpath(tmpdir, sub, "eeg")
                mkpath(dir)
                touch(joinpath(dir, "$(sub)_task-test_eeg.bdf"))
                touch(joinpath(dir, "$(sub)_task-test_eeg.json"))
            end

            # Non-recursive should find nothing (no BDFs in root)
            flat = EegFun.get_files(tmpdir, "\\.bdf")
            @test isempty(flat)

            # Recursive should find both BDFs
            recursive = EegFun.get_files(tmpdir, "\\.bdf", recursive = true)
            @test length(recursive) == 2
            @test all(endswith(".bdf"), recursive)
            @test all(isfile, recursive)

            # Should not match JSON files
            recursive_json = EegFun.get_files(tmpdir, "\\.bdf\$", recursive = true)
            @test all(endswith(".bdf"), recursive_json)
        end
    end

end
