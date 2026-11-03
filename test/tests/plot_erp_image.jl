using Test
using Makie
using EegFun

@testset "plot_erp_image" begin
    dat = EegFun.create_test_epoch_data(n_epochs = 20, n_channels = 3)

    @testset "basic functionality" begin
        result = EegFun.plot_erp_image(dat; display_plot = false)
        @test result.fig isa Figure
        @test length(result.axes) > 0
    end

    @testset "channel selection" begin
        result = EegFun.plot_erp_image(dat; channel_selection = EegFun.channels([1]), display_plot = false)
        @test result.fig isa Figure
    end

    @testset "grid layout" begin
        result = EegFun.plot_erp_image(dat; layout = :grid, display_plot = false)
        @test result.fig isa Figure
        @test length(result.axes) >= 3
    end

    @testset "topo layout" begin
        result = EegFun.plot_erp_image(dat; layout = :topo, display_plot = false)
        @test result.fig isa Figure
        @test length(result.axes) >= 3
    end

    @testset "boxcar average" begin
        result = EegFun.plot_erp_image(dat; boxcar_average = 3, display_plot = false)
        @test result.fig isa Figure
    end

    @testset "channel_plot_order and symbol channel_selection" begin
        result = EegFun.plot_erp_image(dat; channel_selection = :Ch1, channel_plot_order = ["Ch1"], display_plot = false)
        @test result.fig isa Figure
    end

    @testset "symmetric diverging colorrange, non-mutating multi-channel colorrange, and epoch ylim preservation" begin
        dat_scale = EegFun.create_test_epoch_data(n_epochs = 12, n_channels = 2)
        for df in dat_scale.data
            df[!, :Ch1] .= 0.5 .* df[!, :Ch1] .+ 2.0
            df[!, :Ch2] .= 50.0 .* df[!, :Ch2]
        end

        res_multi = EegFun.plot_erp_image(dat_scale; layout = :grid, colormap = :coolwarm, ylim = (-3.0, 3.0), display_plot = false)
        @test res_multi.fig isa Figure
        hm1 = only(filter(p -> p isa Makie.Heatmap, res_multi.axes[1].scene.plots))
        hm2 = only(filter(p -> p isa Makie.Heatmap, res_multi.axes[2].scene.plots))

        # Diverging :coolwarm auto colorrange is symmetric around 0
        cr1 = hm1.colorrange[]
        cr2 = hm2.colorrange[]
        @test cr1[1] ≈ -cr1[2]
        @test cr2[1] ≈ -cr2[2]
        # Channel 2's colorrange was NOT locked to Channel 1's smaller colorrange
        @test cr2[2] > cr1[2] * 5

        # Heatmap axis preserves epoch index limits (1, 12) despite waveform ylim = (-3.0, 3.0)
        hm_ylo = res_multi.axes[1].finallimits[].origin[2]
        hm_yhi = hm_ylo + res_multi.axes[1].finallimits[].widths[2]
        @test hm_ylo <= 1.0
        @test hm_yhi >= 12.0
    end
end
