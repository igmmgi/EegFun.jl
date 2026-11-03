using Test
using Makie
using EegFun

@testset "plot_power_spectrum" begin
    dat = EegFun.create_test_continuous_data(n_channels = 3)

    @testset "plot_channel_spectrum" begin
        result = EegFun.plot_channel_spectrum(dat; display_plot = false)
        @test result.fig isa Figure
        @test length(result.axes) > 0
    end

    @testset "plot_channel_spectrum with kwargs" begin
        result = EegFun.plot_channel_spectrum(
            dat;
            max_freq = 50.0,
            xscale = :linear,
            yscale = :linear,
            unit = :dB,
            show_legend = true,
            display_plot = false,
        )
        @test result.fig isa Figure
        @test result.axes isa Vector{Axis}

        # Backwards-compatible x_scale / y_scale aliases
        result_alias = EegFun.plot_channel_spectrum(dat; max_freq = 50.0, x_scale = :log10, y_scale = :linear, display_plot = false)
        @test result_alias.fig isa Figure
    end

    @testset "OrderedDict freq_bands and ICA spectrum centering + labels" begin
        res_bands = EegFun.plot_channel_spectrum(dat; max_freq = 120.0, show_freq_bands = true, display_plot = false)
        all_axes = filter(b -> b isa Axis, res_bands.fig.content)
        @test length(all_axes) == 2
        band_texts = reduce(
            vcat,
            [
                p.text[] isa AbstractVector ? String.(p.text[]) : [String(p.text[])] for
                p in filter(p -> p isa Makie.Text, all_axes[2].scene.plots)
            ],
        )
        @test band_texts == ["δ", "θ", "α", "β", "γ"]

        ica_res = EegFun.run_ica(dat)
        dat_shifted = copy(dat)
        ica_shifted = deepcopy(ica_res)
        ch1 = ica_shifted.layout.data.label[1]
        dat_shifted.data[!, ch1] .+= 250.0
        ica_shifted.mean[1] += 250.0

        res_orig = EegFun.plot_ica_component_spectrum(dat, ica_res; component_selection = EegFun.components(1:2), display_plot = false)
        res_shift =
            EegFun.plot_ica_component_spectrum(dat_shifted, ica_shifted; component_selection = EegFun.components(1:2), display_plot = false)
        @test res_orig.fig isa Figure
        line_orig = first(filter(p -> p isa Makie.Lines, res_orig.axes[1].scene.plots))
        line_shift = first(filter(p -> p isa Makie.Lines, res_shift.axes[1].scene.plots))
        @test [pt[2] for pt in line_orig[1][]] ≈ [pt[2] for pt in line_shift[1][]]
        @test occursin("%)", line_orig.label[])
    end
end
