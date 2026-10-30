using Test
using Makie
using EegFun
using DataFrames

@testset "plot_databrowser" begin
    dat = EegFun.create_test_continuous_data(n_channels = 3)

    @testset "basic functionality" begin
        fig = EegFun.plot_databrowser(dat; display_plot = false)
        @test fig.fig isa Figure
    end

    @testset "with ICA" begin
        ica_res = EegFun.run_ica(dat)
        fig = EegFun.plot_databrowser(dat, ica_res; display_plot = false)
        @test fig.fig isa Figure
    end

    @testset "epoched data navigation" begin
        epoch_dat, _ = EegFun.create_test_epoch_data_with_artifacts(n_epochs = 5, n_channels = 3)
        res = EegFun.plot_databrowser(epoch_dat; display_plot = false)
        @test res.fig isa Figure

        # Test in-place epoch stepping
        plot_kwargs = EegFun._merge_plot_kwargs(EegFun.PLOT_DATABROWSER_KWARGS, Dict{Symbol,Any}())
        fig = Figure()
        ax = Axis(fig[1, 1])
        state = EegFun._create_browser_state(epoch_dat, epoch_dat.layout.data.label, ax, nothing, plot_kwargs)
        EegFun._draw(ax, state)

        n_plots_start = length(ax.scene.plots)
        EegFun._step_epoch_forward(ax, state)
        @test state.data.current_epoch[] == 2
        @test length(ax.scene.plots) == n_plots_start

        # Check line coordinates match epoch 2
        ch1 = state.channels.labels[1]
        line_ch1 = state.channels.data_lines[ch1]
        y_pts = [p[2] for p in line_ch1[1][]]
        expected_y = epoch_dat.data[2][!, ch1] .* state.view.amplitude_scale[] .+ state.view.offset[1]
        @test maximum(abs.(y_pts .- expected_y)) < 1e-12
    end

    @testset "line points update once per scroll step" begin
        cdat = EegFun.create_test_continuous_data(n = 3000, n_channels = 3, seed = 1)
        plot_kwargs = EegFun._merge_plot_kwargs(EegFun.PLOT_DATABROWSER_KWARGS, Dict{Symbol,Any}())
        fig = Figure()
        ax = Axis(fig[1, 1])
        state = EegFun._create_browser_state(cdat, cdat.layout.data.label, ax, nothing, plot_kwargs)
        EegFun._draw(ax, state)
        state.view.xrange[] = 101:600

        ch = state.channels.labels[2]
        line = state.channels.data_lines[ch]
        n_updates = Ref(0)
        on(line[1]) do _
            n_updates[] += 1
        end
        for _ = 1:5
            state.view.xrange[] = state.view.xrange[] .+ 50
            line[1][]
        end
        @test n_updates[] == 5

        win = state.view.xrange[]
        pts = line[1][]
        @test pts isa Vector{Point2{Float64}}
        @test length(pts) == length(win)
        @test maximum(abs.([p[1] for p in pts] .- cdat.data.time[win])) < 1e-12
        expected_y = cdat.data[win, ch] .* state.view.amplitude_scale[] .+ state.view.offset[2]
        @test maximum(abs.([p[2] for p in pts] .- expected_y)) < 1e-12

        # amplitude scale changes also update the points
        state.view.amplitude_scale[] = 2.5
        expected_y = cdat.data[win, ch] .* 2.5 .+ state.view.offset[2]
        @test maximum(abs.([p[2] for p in line[1][]] .- expected_y)) < 1e-12

        # the point buffer is reused (mutated in place) while the window length is unchanged ...
        buf = line[1][]
        state.view.xrange[] = state.view.xrange[] .+ 25
        @test line[1][] === buf
        # ... and reallocated when it changes
        state.view.xrange[] = first(state.view.xrange[]):(last(state.view.xrange[])-100)
        @test length(line[1][]) == length(state.view.xrange[])
        @test line[1][] !== buf
    end

    @testset "ICA ghost lines" begin
        cdat = EegFun.create_test_continuous_data(n = 3000, n_channels = 4, seed = 2)
        ica_res = EegFun.run_ica(cdat)
        plot_kwargs = EegFun._merge_plot_kwargs(EegFun.PLOT_DATABROWSER_KWARGS, Dict{Symbol,Any}())
        fig = Figure()
        ax = Axis(fig[1, 1])
        state = EegFun._create_browser_state(cdat, cdat.layout.data.label, ax, ica_res, plot_kwargs)
        state.view.xrange[] = 101:500
        append!(state.ica.removed_components, [1, 2])
        EegFun._clear_and_save_limits!(ax, state)
        EegFun._draw(ax, state)

        # Independent dense reference (activations of the unfiltered original data, window only)
        function reference(win, ch_idx)
            labels = ica_res.layout.data.label
            current = state.data.current[].data
            restored = collect(Float64, current[win, labels[ch_idx]])
            for comp in state.ica.removed_components
                act = zeros(length(win))
                for (ci, lbl) in enumerate(labels)
                    act .+= ica_res.unmixing[comp, ci] .* ((cdat.data[win, lbl] .- ica_res.mean[ci]) ./ ica_res.scale)
                end
                restored .+= ica_res.mixing[ch_idx, comp] * ica_res.scale .* act
            end
            return restored, collect(Float64, current[win, labels[ch_idx]])
        end

        ch = state.channels.labels[3]
        orig_line = state.channels.original_lines[ch]
        sub_line = state.channels.subtracted_lines[ch]
        scale = state.view.amplitude_scale[]
        off = state.view.offset[3]

        # Hidden ghosts: NaN, and no activations are computed
        @test all(isnan, [p[2] for p in orig_line[1][]])
        @test all(isnan, [p[2] for p in sub_line[1][]])
        @test isempty(state.ica.activation_cache)
        state.view.xrange[] = 151:550
        orig_line[1][]
        @test isempty(state.ica.activation_cache)

        # Visible ghosts match the dense reference, also after scrolling
        state.view.show_original_ica[] = true
        state.view.show_subtracted_ica[] = true
        for shift in (0, 50, 120)
            state.view.xrange[] = (151:550) .+ shift
            win = state.view.xrange[]
            restored, current_raw = reference(win, 3)
            @test maximum(abs.([p[2] for p in orig_line[1][]] .- (restored .* scale .+ off))) < 1e-9
            @test maximum(abs.([p[2] for p in sub_line[1][]] .- ((restored .- current_raw) .* scale .+ off))) < 1e-9
        end

        # Cache is reused across scroll steps and dropped when removed components change
        entry = state.ica.activation_cache[1]
        state.view.xrange[] = state.view.xrange[] .+ 10
        orig_line[1][]
        @test state.ica.activation_cache[1] === entry
        push!(state.ica.removed_components, 3)
        EegFun._ica_activations!(state, 1, state.view.xrange[])
        @test state.ica.activation_key == [1, 2, 3]
        @test size(state.ica.activation_cache[1][2], 2) == 3

        # Toggling a ghost off returns it to NaN; toggling it back on recomputes it for the current window
        state.view.show_original_ica[] = false
        @test all(isnan, [p[2] for p in orig_line[1][]])
        state.view.xrange[] = state.view.xrange[] .+ 7   # scrolling while hidden leaves it NaN
        @test all(isnan, [p[2] for p in orig_line[1][]])
        state.view.show_original_ica[] = true
        @test !any(isnan, [p[2] for p in orig_line[1][]])
    end

    @testset "redraws do not stack observers" begin
        cdat = EegFun.create_test_continuous_data(n = 3000, n_channels = 4, seed = 4)
        plot_kwargs = EegFun._merge_plot_kwargs(EegFun.PLOT_DATABROWSER_KWARGS, Dict{Symbol,Any}())
        fig = Figure()
        ax = Axis(fig[1, 1])
        state = EegFun._create_browser_state(cdat, cdat.layout.data.label, ax, nothing, plot_kwargs)
        EegFun._draw(ax, state)
        n_xrange = length(state.view.xrange.listeners)
        n_current = length(state.data.current.listeners)
        n_amp = length(state.view.amplitude_scale.listeners)
        for _ = 1:4   # even number: back in the same butterfly mode (labels only exist outside butterfly mode)
            EegFun._butterfly_plot!(ax, state)
        end
        @test length(state.view.xrange.listeners) == n_xrange
        @test length(state.data.current.listeners) == n_current
        @test length(state.view.amplitude_scale.listeners) == n_amp

        # lines still follow the data after redraws
        state.view.xrange[] = 101:600
        state.view.xrange[] = 151:650
        ch = state.channels.labels[1]
        pts = state.channels.data_lines[ch][1][]
        expected_y = cdat.data[151:650, ch] .* state.view.amplitude_scale[] .+ state.view.offset[1]
        @test maximum(abs.([p[2] for p in pts] .- expected_y)) < 1e-12
    end

    @testset "selected channel is visually distinct with ICA" begin
        cdat = EegFun.create_test_continuous_data(n = 3000, n_channels = 4, seed = 5)
        ica_res = EegFun.run_ica(cdat)
        plot_kwargs = EegFun._merge_plot_kwargs(EegFun.PLOT_DATABROWSER_KWARGS, Dict{Symbol,Any}())
        fig = Figure()
        ax = Axis(fig[1, 1])
        state = EegFun._create_browser_state(cdat, cdat.layout.data.label, ax, ica_res, plot_kwargs)
        state.view.show_original_ica[] = true   # default in the ICA menu
        EegFun._clear_and_save_limits!(ax, state)
        EegFun._draw(ax, state)
        ch = state.channels.labels[2]
        line = state.channels.data_lines[ch]
        unselected_cmap = line.colormap[]
        unselected_width = line.linewidth[]
        EegFun._toggle_channel_visibility!(ax, state, 2)
        @test line.colormap[][1] == plot_kwargs[:selected_channel_color]
        @test line.colormap[][1] != unselected_cmap[1]
        @test line.linewidth[] > unselected_width
        EegFun._toggle_channel_visibility!(ax, state, 2)
        @test line.colormap[][1] == unselected_cmap[1]
        @test line.linewidth[] == unselected_width
    end

    @testset "_add_marker! has a docstring" begin
        @test !occursin("No documentation found", string(@doc EegFun._add_marker!))
    end

    @testset "path sanitization" begin
        @test EegFun._percent_decode("hello%20world") == "hello world"
        @test EegFun._percent_decode("test%C3%A4%C3%B6%C3%BC") == "testäöü"
        @test EegFun._sanitize_dropped_path("file:///home/ian/test%20file.jld2") == "/home/ian/test file.jld2"
        @test EegFun._sanitize_dropped_path("file://localhost/home/ian/test.jld2") == "/home/ian/test.jld2"
        @test EegFun._sanitize_dropped_path("/path/with%20percent.jld2") == "/path/with%20percent.jld2"
    end

    @testset "channel repair survives filter toggle" begin
        cdat = EegFun.create_test_continuous_data(n = 3000, n_channels = 5, seed = 42)
        plot_kwargs = EegFun._merge_plot_kwargs(EegFun.PLOT_DATABROWSER_KWARGS, Dict{Symbol,Any}())
        fig = Figure()
        ax = Axis(fig[1, 1])
        state = EegFun._create_browser_state(cdat, cdat.layout.data.label, ax, nothing, plot_kwargs)

        ch_bad = state.channels.labels[2]

        EegFun._repair_selected_channels!(state, [ch_bad], :spherical_spline, ax)
        @test length(state.channel_repair_history) == 1
        repaired_data = copy(state.data.current[].data[!, ch_bad])

        # Toggle HP filter on
        EegFun._apply_hp_filter!(state, true)
        @test state.data.filter_state.hp_active[] == true
        @test state.analysis_settings[].repaired_channels == [ch_bad]
        @test state.analysis_settings[].repair_method == :spherical_spline
        @test state.analysis_settings[].hp_filter > 0.0

        # Toggle HP filter off
        EegFun._apply_hp_filter!(state, false)
        @test state.data.filter_state.hp_active[] == false
        @test state.data.current[].data[!, ch_bad] ≈ repaired_data

        # Undo repair
        EegFun._undo_last_repair!(state, ax)
        @test isempty(state.channel_repair_history)
        @test isempty(state.analysis_settings[].repaired_channels)
    end

    @testset "ICA with channel subset does not throw KeyError in _draw" begin
        cdat = EegFun.create_test_continuous_data(n = 3000, n_channels = 4, seed = 7)
        # Run ICA on a subset of 3 channels
        ica_res = EegFun.run_ica(cdat; channel_selection = EegFun.channels(1:3))

        plot_kwargs = EegFun._merge_plot_kwargs(EegFun.PLOT_DATABROWSER_KWARGS, Dict{Symbol,Any}())
        fig = Figure()
        ax = Axis(fig[1, 1])
        # Browser on full data with 4 channels, but ICA on 3 channels
        state = EegFun._create_browser_state(cdat, cdat.layout.data.label, ax, ica_res, plot_kwargs)
        state.view.show_original_ica[] = true
        state.view.show_subtracted_ica[] = true
        push!(state.ica.removed_components, 1)

        # _draw must not throw KeyError for the 4th channel
        @test_nowarn EegFun._draw(ax, state)
    end

    @testset "ICA component removal with active filter matches filter after removal" begin
        cdat = EegFun.create_test_continuous_data(n = 3000, n_channels = 4, seed = 123)
        ica_res = EegFun.run_ica(cdat)
        plot_kwargs = EegFun._merge_plot_kwargs(EegFun.PLOT_DATABROWSER_KWARGS, Dict{Symbol,Any}())
        fig = Figure()
        ax = Axis(fig[1, 1])
        state = EegFun._create_browser_state(cdat, cdat.layout.data.label, ax, ica_res, plot_kwargs)

        EegFun._apply_hp_filter!(state, true)
        EegFun._toggle_ica_component!(state, 1)

        cdat_expected = copy(cdat)
        EegFun._apply_single_ica_removal!(cdat_expected, cdat, ica_res, 1)
        sym_method = Symbol(state.data.filter_state.hp_method[])
        sym_func = Symbol(state.data.filter_state.hp_func[])
        EegFun.highpass_filter!(
            cdat_expected,
            state.data.filter_state.hp_freq[];
            filter_method = sym_method,
            order = state.data.filter_state.hp_order[],
            filter_func = sym_func,
        )

        for col in state.channels.labels
            diff = maximum(abs.(state.data.current[].data[!, col] .- cdat_expected.data[!, col]))
            @test diff < 1e-12
        end
    end

    @testset "analysis_settings copies mutable state" begin
        cdat = EegFun.create_test_continuous_data(n = 2000, n_channels = 3, seed = 8)
        plot_kwargs = EegFun._merge_plot_kwargs(EegFun.PLOT_DATABROWSER_KWARGS, Dict{Symbol,Any}())
        fig = Figure()
        ax = Axis(fig[1, 1])
        state = EegFun._create_browser_state(cdat, cdat.layout.data.label, ax, nothing, plot_kwargs)

        push!(state.selection.selected_regions[], (0.1, 0.5))
        EegFun._update_analysis_settings!(state)
        settings = state.analysis_settings[]
        @test length(settings.selected_regions) == 1
        push!(state.selection.selected_regions[], (0.6, 0.8))
        @test length(settings.selected_regions) == 1

        EegFun._clear_all_selected_regions!(ax, state)
        @test isempty(state.analysis_settings[].selected_regions)
    end

    @testset "selection zero-width edge case" begin
        cdat = EegFun.create_test_continuous_data(n = 2000, n_channels = 3, seed = 9)
        plot_kwargs = EegFun._merge_plot_kwargs(EegFun.PLOT_DATABROWSER_KWARGS, Dict{Symbol,Any}())
        fig = Figure()
        ax = Axis(fig[1, 1])
        state = EegFun._create_browser_state(cdat, cdat.layout.data.label, ax, nothing, plot_kwargs)

        EegFun._start_selection!(ax, state, 0.5)
        @test state.selection.active[] == true
        EegFun._finish_selection!(ax, state, 0.5)
        @test state.selection.active[] == false
        @test isempty(state.selection.selected_regions[])
    end

    @testset "vector of datasets returns results" begin
        dats = [EegFun.create_test_continuous_data(n = 500, n_channels = 2, seed = i) for i = 1:2]
        res = EegFun.plot_databrowser(dats; display_plot = false)
        @test res isa Vector{<:NamedTuple}
        @test length(res) == 2
        @test haskey(res[1], :fig)
        @test haskey(res[1], :analysis_settings)
    end

    @testset "Data - ICA line color is green with ICA" begin
        cdat = EegFun.create_test_continuous_data(n = 1000, n_channels = 3, seed = 10)
        ica_res = EegFun.run_ica(cdat)
        plot_kwargs = EegFun._merge_plot_kwargs(EegFun.PLOT_DATABROWSER_KWARGS, Dict{Symbol,Any}())
        fig = Figure()
        ax = Axis(fig[1, 1])
        state = EegFun._create_browser_state(cdat, cdat.layout.data.label, ax, ica_res, plot_kwargs)
        EegFun._draw(ax, state)
        ch = state.channels.labels[1]
        line = state.channels.data_lines[ch]
        # Default with ICA is unselected channel color (only data shown)
        @test line.colormap[][1] == plot_kwargs[:unselected_channel_color]
        # When show_cleaned_ica is enabled, line colormap switches to green
        state.plot_kwargs[:show_cleaned_ica][] = true
        @test line.colormap[][1] == :green
    end

    @testset "plot_databrowser with ICA defaults to only data shown" begin
        cdat = EegFun.create_test_continuous_data(n = 1000, n_channels = 3, seed = 11)
        ica_res = EegFun.run_ica(cdat)
        plot_kwargs = EegFun._merge_plot_kwargs(EegFun.PLOT_DATABROWSER_KWARGS, Dict{Symbol,Any}())
        fig = Figure()
        ax = Axis(fig[1, 1])
        state = EegFun._create_browser_state(cdat, cdat.layout.data.label, ax, ica_res, plot_kwargs)
        EegFun._setup_ui(fig, ax, state, cdat, ica_res, plot_kwargs)
        EegFun._draw(ax, state)

        # Default state: only Data is shown, not Data - ICA, nor ICA activations
        @test state.view.show_original_ica[] == true
        @test state.plot_kwargs[:show_cleaned_ica][] == false
        @test state.view.show_subtracted_ica[] == false

        # Visual lines:
        ch = state.channels.labels[1]
        data_line = state.channels.data_lines[ch]
        orig_line = state.channels.original_lines[ch]
        sub_line = state.channels.subtracted_lines[ch]

        @test data_line.visible[] == true
        @test data_line.colormap[][1] == plot_kwargs[:unselected_channel_color]
        @test orig_line.visible[] == false
        @test sub_line.visible[] == false
    end

    @testset "topography scrubber on selected data" begin
        lay = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
        EegFun.polar_to_cartesian_xy!(lay)
        cdat = EegFun.create_test_continuous_data(n = 200, n_channels = 72, seed = 42)
        cdat.layout = lay
        names_dict = Dict(Symbol("Ch$i") => Symbol(lay.data.label[i]) for i = 1:72)
        DataFrames.rename!(cdat.data, names_dict)

        sub = EegFun.subset(cdat, sample_selection = x -> 0.0 .<= x.time .<= 0.5)
        # 2D scrubber
        res_2d = EegFun.plot_topography(sub; scrubber = true, display_plot = false)
        @test res_2d.fig isa Figure
        @test length(res_2d.axes) == 1

        # 3D scrubber
        res_3d_scrubber = EegFun.plot_topography_3d(sub; scrubber = true, display_plot = false)
        @test res_3d_scrubber.fig isa Figure
        @test length(res_3d_scrubber.axes) == 1

        # 3D static (default)
        res_3d_static = EegFun.plot_topography_3d(sub; scrubber = false, display_plot = false)
        @test res_3d_static.fig isa Figure
        @test length(res_3d_static.axes) == 1
    end

    @testset "region trigger detection and epoch plotting" begin
        cdat = EegFun.create_test_continuous_data(n = 2000, n_channels = 3, seed = 55)
        # Add triggers: code 10 at sample 500 and 1500, code 20 at sample 1800
        cdat.data.trigger[500] = 10
        cdat.data.trigger[1500] = 10
        cdat.data.trigger[1800] = 20

        plot_kwargs = EegFun._merge_plot_kwargs(EegFun.PLOT_DATABROWSER_KWARGS, Dict{Symbol,Any}())
        fig = Figure()
        ax = Axis(fig[1, 1])
        state = EegFun._create_browser_state(cdat, cdat.layout.data.label, ax, nothing, plot_kwargs)

        # 1. No region selected -> empty trigger events
        @test isempty(EegFun._find_region_triggers(state))

        # 2. Select region containing no triggers (samples 1 to 200)
        t_start = cdat.data.time[1]
        t_end = cdat.data.time[200]
        push!(state.selection.selected_regions[], (t_start, t_end))
        @test isempty(EegFun._find_region_triggers(state))

        # 3. Select region containing trigger 10 (around sample 500)
        t_trig = cdat.data.time[500]
        x_min = t_trig - 0.2
        x_max = t_trig + 0.4
        push!(state.selection.selected_regions[], (x_min, x_max))

        trig_events = EegFun._find_region_triggers(state)
        @test length(trig_events) == 1
        @test trig_events[1].sample == 500
        @test trig_events[1].code == 10
        @test abs(trig_events[1].time - t_trig) < 1e-12

        # 4. Extract single epoch
        single_res = EegFun._plot_single_trigger_epoch(state, trig_events[1], x_min, x_max; display_plot = false)
        @test haskey(single_res, :fig)
        @test single_res.fig isa Figure

        # 5. Extract all epochs for trigger 10
        all_res = EegFun._plot_all_trigger_epochs(state, 10, x_min, x_max, t_trig; display_plot = false)
        @test haskey(all_res, :fig)
        @test all_res.fig isa Figure

        # 6. Select region containing multiple triggers (around 1500 and 1800)
        t_trig2 = cdat.data.time[1500]
        t_trig3 = cdat.data.time[1800]
        push!(state.selection.selected_regions[], (t_trig2 - 0.1, t_trig3 + 0.1))
        multi_events = EegFun._find_region_triggers(state)
        @test length(multi_events) == 2
        @test multi_events[1].code == 10
        @test multi_events[2].code == 20
    end
end
