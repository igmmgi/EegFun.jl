using Test
using DataFrames
using Makie
using EegFun

@testset "Makie Recipes" begin
    @testset "ErpData recipes" begin
        erp = EegFun.create_test_erp_data(n_channels = 3)
        fig = Figure()
        ax1 = Axis(fig[1, 1])
        plt_line_pos = lines!(ax1, erp, :Ch1)
        @test plt_line_pos isa Makie.Lines

        ax2 = Axis(fig[1, 2])
        plt_series = series!(ax2, erp)
        @test plt_series isa Makie.Series
    end

    @testset "SpectrumData recipes use freq column" begin
        cdat = EegFun.create_test_continuous_data(n = 1024, fs = 256, n_channels = 3, seed = 1)
        spec = EegFun.freq_spectrum(cdat; max_freq = 50.0)

        pts_tuple = Makie.convert_arguments(Makie.PointBased(), spec, :Ch1)
        pts = pts_tuple[1]
        @test length(pts) == nrow(spec.data)
        @test [p[1] for p in pts] ≈ spec.data.freq
        @test [p[2] for p in pts] ≈ spec.data.Ch1

        fig = Figure()
        ax1 = Axis(fig[1, 1])
        @test lines!(ax1, spec, :Ch1) isa Makie.Lines

        ax2 = Axis(fig[1, 2])
        @test series!(ax2, spec) isa Makie.Series
    end

    @testset "TimeFreqData throws clear error directing to plot_tf" begin
        ep = EegFun.create_test_epoch_data(n_epochs = 4, n_channels = 2)
        tf = EegFun.tf_morlet(ep; frequencies = [8.0, 12.0], cycles = 3)
        @test_throws EegFun.EegFunError Makie.convert_arguments(Makie.PointBased(), tf)
        @test_throws EegFun.EegFunError Makie.convert_arguments(Makie.PointBased(), tf, :Ch1)
        @test_throws EegFun.EegFunError Makie.convert_arguments(Makie.Series, tf)
    end

    @testset "Layout recipe computes 2D coordinates without mutating input" begin
        raw_df = DataFrame(label = [:Fz, :Cz, :Pz], inc = [30.0, 0.0, 30.0], azi = [0.0, 0.0, 180.0])
        raw_layout = EegFun.Layout(raw_df, nothing, nothing, nothing)
        @test !hasproperty(raw_layout.data, :x2)

        pts_tuple = Makie.convert_arguments(Makie.PointBased(), raw_layout)
        pts = pts_tuple[1]
        @test length(pts) == 3
        # Input layout was not mutated
        @test !hasproperty(raw_layout.data, :x2)

        fig = Figure()
        ax = Axis(fig[1, 1])
        @test scatter!(ax, raw_layout) isa Makie.Scatter
    end
end
