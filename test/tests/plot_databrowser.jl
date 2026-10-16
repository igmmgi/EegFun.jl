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
end
