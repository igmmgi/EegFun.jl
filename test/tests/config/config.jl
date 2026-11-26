# Tests for TOML configuration loading, default merging, validation, and template generation.

using Test
using TOML
using Dates
using Logging
using EegFun

@testset "Configuration System Tests" begin
    mktempdir() do test_dir
        load_toml_str(tag, content) = (p = joinpath(test_dir, "$tag.toml"); write(p, content); EegFun.read_config(p))

        @testset "ConfigParameter, constructors & PARAMETERS registry" begin
            p_f = EegFun.ConfigParameter{Float64}(description = "Test parameter", default = 1.0, min = 0.0, max = 2.0)
            @test (p_f.description, p_f.default, p_f.min, p_f.max, p_f.allowed) == ("Test parameter", 1.0, 0.0, 2.0, nothing)
            p_s = EegFun.ConfigParameter{String}(description = "Minimal")
            @test all(isnothing, (p_s.default, p_s.min, p_s.max, p_s.allowed))
            @test fieldtype(EegFun.ConfigParameter{String}, :allowed) == Union{Nothing,Vector{Any}}

            @test EegFun._string_param("s", "def", allowed = ["a", "b"]).allowed == ["a", "b"]
            @test EegFun._string_param("s").default == ""
            @test EegFun._string_param("s", "def", allowed = String[]).allowed == String[]
            @test EegFun._bool_param("b", true).default == true && EegFun._bool_param("b").default == false
            @test EegFun._number_param("n", 5.0, -10.0, 10.0).min == -10.0
            @test isnothing(EegFun._number_param("n", 5.0).min)
            @test EegFun._channel_groups_param("ch", [["Fp1"], ["Fp2"]]).default == [["Fp1"], ["Fp2"]]
            @test EegFun._param(Int, "i", 42; min = 0, max = 100) isa EegFun.ConfigParameter{Int}
            @test isnothing(EegFun._param(Union{String,Nothing}, "u", nothing).default)
            @test isnothing(EegFun._param(Bool, "b_none").default)

            fspec = EegFun._filter_param_spec("test.prefix", true, :hp, 1.0, 0.1, 10.0, 2, 1, 5)
            for sfx in ("apply", "type", "method", "func", "freq", "order")
                @test haskey(fspec, "test.prefix.$sfx")
            end
            @test (fspec["test.prefix.type"].default, fspec["test.prefix.freq"].default, fspec["test.prefix.order"].max) == (:hp, 1.0, 5)

            required_keys = [
                "files.input.directory",
                "files.input.raw_data_files",
                "files.input.layout_file",
                "files.output.directory",
                "files.output.save_ica_data",
                "files.output.save_epoch_data_raw",
                "files.output.save_epoch_data_corrected",
                "files.output.save_epoch_data",
                "files.output.save_erp_data_raw",
                "files.output.save_erp_data_corrected",
                "files.output.save_erp_data",
                "preprocess.filter.highpass.apply",
                "preprocess.filter.highpass.method",
                "preprocess.filter.highpass.func",
                "preprocess.filter.highpass.freq",
                "preprocess.filter.lowpass.apply",
                "preprocess.filter.lowpass.method",
                "preprocess.filter.lowpass.func",
                "preprocess.filter.lowpass.freq",
                "preprocess.filter.ica_highpass.apply",
                "preprocess.filter.ica_lowpass.apply",
                "preprocess.epoch_start",
                "preprocess.epoch_end",
                "preprocess.reference_channel",
                "preprocess.layout.neighbor_criterion",
                "preprocess.channel_repair.method",
                "preprocess.eog.vEOG_channels",
                "preprocess.eog.hEOG_channels",
                "preprocess.eog.vEOG_criterion",
                "preprocess.eog.hEOG_criterion",
                "preprocess.eeg.extreme_value_abs_criterion",
                "preprocess.eeg.artifact_value_abs_criterion",
                "preprocess.eeg.artifact_interval_end",
                "preprocess.ica.apply",
                "preprocess.ica.component_method",
            ]
            @test all(k -> haskey(EegFun.PARAMETERS, k), required_keys) && !haskey(EegFun.PARAMETERS, "nonexistent")
            @test all(p -> p isa EegFun.ConfigParameter && !isempty(p.description), values(EegFun.PARAMETERS))
            @test EegFun.PARAMETERS["files.input.directory"] isa EegFun.ConfigParameter{String}
            @test EegFun.PARAMETERS["preprocess.filter.highpass.freq"] isa EegFun.ConfigParameter{Real}
            @test EegFun.PARAMETERS["files.output.save_continuous_data_raw"] isa EegFun.ConfigParameter{Bool}
            @test EegFun.PARAMETERS["preprocess.channel_repair.method"].default == "spherical_spline"
            @test EegFun.PARAMETERS["preprocess.layout.neighbor_criterion"].default == 0.35
            for k in ("preprocess.eeg.extreme_value_abs_criterion", "preprocess.eeg.artifact_value_abs_criterion")
                @test contains(EegFun.PARAMETERS[k].description, "μV") && !contains(EegFun.PARAMETERS[k].description, "(mV)")
            end
        end

        @testset "read_config merging, coercion, boundaries & validation" begin
            cfg = load_toml_str(
                "valid",
                """
                [files.input]
                directory = "/path/with spaces/and/special@chars"
                raw_data_files = ["file1.bdf", "file2.bdf"]
                epoch_condition_file = ""

                [files.output]
                directory = "/custom/output"
                save_erp_data_raw = false
                save_continuous_data_raw = true
                save_ica_data = false

                [preprocess]
                epoch_start = -2.0
                epoch_end = 3.0

                [preprocess.filter.highpass]
                freq = 1
                order = 2.0
                apply = true
                """,
            )
            @test cfg["files"]["input"]["directory"] == "/path/with spaces/and/special@chars"
            @test cfg["files"]["input"]["raw_data_files"] == ["file1.bdf", "file2.bdf"]
            @test cfg["files"]["input"]["epoch_condition_file"] == ""
            @test cfg["files"]["output"]["directory"] == "/custom/output"
            @test cfg["files"]["output"]["save_erp_data_raw"] == false && cfg["files"]["output"]["save_ica_data"] == false
            @test cfg["preprocess"]["filter"]["highpass"]["freq"] == 1.0 && cfg["preprocess"]["filter"]["highpass"]["order"] == 2
            hp_typed = EegFun.PreprocessConfig(cfg["preprocess"]).filter.highpass
            @test hp_typed.freq === 1.0 && hp_typed.order === 2
            @test cfg["preprocess"]["filter"]["highpass"]["method"] == :iir
            @test cfg["preprocess"]["filter"]["lowpass"]["freq"] == 30
            @test cfg["preprocess"]["reference_channel"] == "avg"

            # Boundary min/max & empty config defaults
            c_min = load_toml_str("min_b", "[preprocess.filter.highpass]\nfreq = 0.01\norder = 1\n")
            c_max = load_toml_str("max_b", "[preprocess.filter.highpass]\nfreq = 20.0\norder = 4\n")
            @test (c_min["preprocess"]["filter"]["highpass"]["freq"], c_max["preprocess"]["filter"]["highpass"]["freq"]) == (0.01, 20.0)
            c_emp = load_toml_str("empty", "# empty\n")
            @test (c_emp["preprocess"]["epoch_start"], c_emp["preprocess"]["epoch_end"]) == (-1, 1)

            # Error cases
            @test_throws Exception EegFun.read_config(joinpath(test_dir, "missing_nope.toml"))
            for (i, bad_toml) in enumerate((
                "[section\nkey = value\n",
                "[preprocess.filter.highpass]\nfreq = -5.0\n",
                "[preprocess.filter.highpass]\nfreq = 0.005\n",
                "[preprocess.filter.highpass]\nfreq = \"not_a_number\"\n",
                "[preprocess.filter.highpass]\norder = \"not_a_number\"\n",
                "[preprocess.filter.highpass]\ntype = \"invalid_type\"\n",
                "[preprocess.filter.highpass]\nmethod = \"invalid\"\n",
                "[preprocess.layout]\nneighbour_criterion = 0.44\n",
            ))
                @test_throws Exception load_toml_str("bad_$i", bad_toml)
            end
        end

        @testset "Validation, section grouping, merging, template I/O & parameter display" begin
            @test EegFun.ValidationResult(success = true).success
            vr_err = EegFun.ValidationResult(success = false, error = "Err", key_path = "a.b")
            @test !vr_err.success && vr_err.error == "Err" && vr_err.key_path == "a.b"

            hp_dict(ov...) = Dict{String,Any}(
                "preprocess" => Dict{String,Any}(
                    "filter" => Dict{String,Any}(
                        "highpass" => merge(
                            Dict{String,Any}("apply" => true, "method" => :iir, "freq" => 0.1, "order" => 1),
                            Dict{String,Any}(ov...),
                        ),
                    ),
                ),
            )
            @test EegFun._validate_config(hp_dict()).success
            @test !EegFun._validate_config(hp_dict("apply" => "true")).success
            @test !EegFun._validate_config(hp_dict("freq" => -1)).success
            @test !EegFun._validate_config(hp_dict("method" => "invalid")).success

            p_num = EegFun.ConfigParameter{Real}(description = "Test", default = 5.0, min = 0.0, max = 10.0)
            @test EegFun._validate_parameter(5.0, p_num, "k").success
            @test !EegFun._validate_parameter(-1.0, p_num, "k").success
            @test !EegFun._validate_parameter(15.0, p_num, "k").success
            @test !EegFun._validate_parameter("str", p_num, "k").success
            @test EegFun._validate_parameter(5, EegFun.ConfigParameter{Float64}(description = "f", default = 5.0), "k").success
            p_all = EegFun.ConfigParameter{String}(description = "s", default = "a", allowed = ["a", "b"])
            @test EegFun._validate_parameter("b", p_all, "k").success && !EegFun._validate_parameter("d", p_all, "k").success

            # Grouping & merging helpers
            sections = EegFun._group_parameters_by_section()
            @test haskey(sections["files"], "input") && haskey(sections["preprocess"], "filter.highpass")
            @test EegFun._extract_subsection("files", "files.input.directory") == "input"
            @test EegFun._extract_subsection("preprocess.filter", "preprocess.filter.ica_highpass.freq") == "ica_highpass"
            @test EegFun._extract_subsection("preprocess.ica", "preprocess.ica.apply") == ""
            @test EegFun._extract_subsection("files", "preprocess.filter.highpass.freq") == ""
            grp = EegFun._group_params_by_subsection(
                "preprocess.filter",
                ["preprocess.filter.highpass.freq", "preprocess.filter.lowpass.freq"],
            )
            @test length(grp["highpass"]) == 1 && length(grp["lowpass"]) == 1

            @test isempty(EegFun._merge_configs(Dict(), Dict()))
            merged = EegFun._merge_configs(Dict("a" => Dict("b" => 1, "c" => 2)), Dict("a" => Dict("b" => 3, "d" => 4)))
            @test (merged["a"]["b"], merged["a"]["c"], merged["a"]["d"]) == (3, 2, 4)

            # Template & writer helpers
            tpl_path = joinpath(test_dir, "template.toml")
            EegFun.generate_config_template(filename = tpl_path)
            @test haskey(TOML.parsefile(tpl_path), "files") && contains(read(tpl_path, String), "[preprocess.ica]")

            fmt_val(v) = (io = IOBuffer(); EegFun._write_parameter_value(io, "p", v); String(take!(io)))
            @test contains(fmt_val("a\\b"), "a\\\\b")
            @test contains(fmt_val([["Fp1", "IO1"]]), "Fp1")
            @test fmt_val(String[]) == "p = []\n" && fmt_val(nothing) == "# p = \n"
            @test fmt_val(42) == "p = 42\n" && fmt_val(3.14) == "p = 3.14\n" && fmt_val(true) == "p = true\n"

            io_doc = IOBuffer()
            EegFun._write_parameter_docs(
                io_doc,
                EegFun.ConfigParameter{Real}(description = "D", default = 5.0, min = 0.0, max = 10.0, allowed = ["a", "b"]),
            )
            doc_str = String(take!(io_doc))
            @test all(s -> contains(doc_str, s), ("# D", "# Type: Real", "0.0 ≤ value ≤ 10.0", "Allowed values: a, b", "# Default: 5.0"))
            io_req = IOBuffer()
            EegFun._write_parameter_docs(io_req, EegFun.ConfigParameter{String}(description = "Req", default = nothing))
            @test contains(String(take!(io_req)), "# [REQUIRED]")

            # Execute display/info functions cleanly (replacing placeholder typeof checks)
            with_logger(NullLogger()) do
                redirect_stdout(devnull) do
                    EegFun._show_parameter_details("preprocess.filter.highpass.freq", EegFun.PARAMETERS["preprocess.filter.highpass.freq"])
                    EegFun._display_section("files", sections["files"])
                    EegFun._display_subsection("input", sections["files"]["input"])
                    EegFun._show_specific_parameter("preprocess.filter.highpass.freq")
                    EegFun._show_specific_parameter("filter")
                    EegFun._show_specific_parameter("nonexistent_param_xyz")
                    EegFun._show_all_parameters()
                    @test EegFun.show_parameter_info() === nothing
                    @test EegFun.show_parameter_info("filter") === nothing
                end
            end
        end

        @testset "Typed config structs, rereference, nullable fields & roundtrip" begin
            def_dict = EegFun._generate_default_config()
            p_def = EegFun.PreprocessConfig(def_dict["preprocess"])
            @test p_def.channel_repair_method == :spherical_spline

            cust_prep = deepcopy(def_dict["preprocess"])
            cust_prep["channel_repair"] = Dict("method" => "neighbor_interpolation")
            cust_prep["cleanline"]["line_frequencies"] = [50, 100]
            p_cust = EegFun.PreprocessConfig(cust_prep)
            @test p_cust.channel_repair_method == :neighbor_interpolation
            @test p_cust.cleanline.line_frequencies == [50.0, 100.0] && p_cust.cleanline.line_frequencies isa Vector{Float64}

            files_cfg = EegFun.FilesConfig(
                Dict(
                    "input" => Dict(
                        "directory" => "./data",
                        "raw_data_files" => "\\.bdf",
                        "recursive" => false,
                        "layout_file" => "biosemi72.csv",
                        "epoch_condition_file" => "epochs.toml",
                    ),
                    "output" => Dict("directory" => "./out", "save_continuous_data_raw" => true, "save_epoch_data_raw" => false),
                ),
            )
            @test (files_cfg.input.directory, files_cfg.input.layout_file, files_cfg.output.directory) ==
                  ("./data", "biosemi72.csv", "./out")
            @test files_cfg.output.save_continuous_data_raw && !files_cfg.output.save_epoch_data_raw

            p_full = EegFun.PipelineConfig(
                Dict(
                    "files" => Dict("input" => Dict("directory" => "./raw_dir")),
                    "preprocess" => Dict("channel_repair" => Dict("method" => "spherical_spline")),
                ),
            )
            @test p_full.files.input.directory == "./raw_dir" && p_full.preprocess.channel_repair_method == :spherical_spline

            @test EegFun.EegConfig().rereference == :avg
            @test EegFun.EegConfig(Dict("rereference" => ["M1", "M2"])).rereference == [:M1, :M2]
            ref_pcfg = EegFun.PipelineConfig(
                load_toml_str(
                    "linked_ref",
                    "[preprocess]\nreference_channel = [\"M1\", \"M2\"]\n[preprocess.ica]\ncomponent_method = \"correlation\"\n",
                ),
            )
            @test ref_pcfg.preprocess.reference_channel == [:M1, :M2] && ref_pcfg.preprocess.ica.component_method == :correlation

            @test EegFun.AutoRejectConfig(Dict("seed" => nothing)).seed === nothing
            eeg_none = EegFun.EegConfig(Dict("artifact_interval_start" => nothing, "artifact_interval_end" => nothing))
            @test isnothing(eeg_none.artifact_interval_start) && isnothing(eeg_none.artifact_interval_end)

            def_from_toml = EegFun.PipelineConfig(load_toml_str("defaults_align", "# empty\n"))
            def_struct = EegFun.PipelineConfig()
            @test def_from_toml.preprocess.neighbor_criterion == def_struct.preprocess.neighbor_criterion
            @test def_from_toml.preprocess.filter.highpass.order == def_struct.preprocess.filter.highpass.order
            @test def_from_toml.preprocess.filter.lowpass.freq == def_struct.preprocess.filter.lowpass.freq
            @test def_from_toml.files.output.save_ica_data == def_struct.files.output.save_ica_data

            rt_toml = joinpath(test_dir, "roundtrip_meta.toml")
            EegFun.print_config(def_dict, rt_toml)
            rt_loaded = EegFun.read_config(rt_toml)
            @test haskey(rt_loaded, "metadata") && EegFun.PipelineConfig(rt_loaded) isa EegFun.PipelineConfig
        end
    end
end
