using Test
using Makie
using EegFun

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
        plot_kwargs = EegFun._merge_plot_kwargs(EegFun.PLOT_DATABROWSER_KWARGS, Dict{Symbol, Any}())
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
        plot_kwargs = EegFun._merge_plot_kwargs(EegFun.PLOT_DATABROWSER_KWARGS, Dict{Symbol, Any}())
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
        for _ in 1:5
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
        state.view.xrange[] = first(state.view.xrange[]):(last(state.view.xrange[]) - 100)
        @test length(line[1][]) == length(state.view.xrange[])
        @test line[1][] !== buf
    end

    @testset "ICA ghost lines" begin
        cdat = EegFun.create_test_continuous_data(n = 3000, n_channels = 4, seed = 2)
        ica_res = EegFun.run_ica(cdat)
        plot_kwargs = EegFun._merge_plot_kwargs(EegFun.PLOT_DATABROWSER_KWARGS, Dict{Symbol, Any}())
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
        plot_kwargs = EegFun._merge_plot_kwargs(EegFun.PLOT_DATABROWSER_KWARGS, Dict{Symbol, Any}())
        fig = Figure()
        ax = Axis(fig[1, 1])
        state = EegFun._create_browser_state(cdat, cdat.layout.data.label, ax, nothing, plot_kwargs)
        EegFun._draw(ax, state)
        n_xrange = length(state.view.xrange.listeners)
        n_current = length(state.data.current.listeners)
        n_amp = length(state.view.amplitude_scale.listeners)
        for _ in 1:4   # even number: back in the same butterfly mode (labels only exist outside butterfly mode)
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
        plot_kwargs = EegFun._merge_plot_kwargs(EegFun.PLOT_DATABROWSER_KWARGS, Dict{Symbol, Any}())
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
end
