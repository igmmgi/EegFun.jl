using Test
using Makie
using EegFun

@testset "plot_decoding" begin
    dat = EegFun.create_test_epoch_data_vector(conditions = 1:2, n_epochs = 10, n_channels = 3)

    res = EegFun.decode_libsvm(dat, n_iterations = 2)

    @testset "basic functionality" begin
        result = EegFun.plot_decoding(res; color = nothing, linestyle = nothing, display_plot = false)
        @test result.fig isa Figure
        @test result.axes isa Vector{Axis}
        @test length(result.axes) > 0
    end

    stats_res = EegFun.test_against_chance([res, res])

    @testset "with stats and significance kwargs" begin
        result = EegFun.plot_decoding(res, stats_res; show_significance = true, significance_position = :top, display_plot = false)
        @test result.fig isa Figure
        @test result.axes isa Vector{Axis}
        @test length(result.axes) == 1

        # Alias :sig_bar_position and show_significance = false
        result_no_sig = EegFun.plot_decoding(res, stats_res; show_significance = false, sig_bar_position = :bottom, display_plot = false)
        @test result_no_sig.fig isa Figure
    end

    @testset "multi-subject shared limits" begin
        res_high = deepcopy(res)
        res_high.average_score .= 0.95
        res_low = deepcopy(res)
        res_low.average_score .= 0.52
        multi = EegFun.plot_decoding([res_low, res_high]; display_plot = false)
        @test multi.fig isa Figure
        @test length(multi.axes) == 2
        yl1 = multi.axes[1].finallimits[].origin[2], multi.axes[1].finallimits[].origin[2] + multi.axes[1].finallimits[].widths[2]
        yl2 = multi.axes[2].finallimits[].origin[2], multi.axes[2].finallimits[].origin[2] + multi.axes[2].finallimits[].widths[2]
        @test yl1[1] ≈ yl2[1]
        @test yl1[2] ≈ yl2[2]
    end

    @testset "asymmetric confusion matrix orientation" begin
        res_cm = deepcopy(res)
        # Row 1 = True Class 1 ([0.9, 0.1]), Row 2 = True Class 2 ([0.4, 0.6])
        asym = [0.9 0.1; 0.4 0.6]
        for t in axes(res_cm.confusion_matrix, 1)
            res_cm.confusion_matrix[t, :, :] .= asym
        end
        cm_out = EegFun.plot_confusion_matrix(res_cm; display_plot = false)
        @test cm_out.fig isa Figure
        @test cm_out.axes isa Vector{Axis}
        ax_cm = cm_out.axes[1]
        @test ax_cm.yreversed[] == true

        hm_plot = only(filter(p -> p isa Makie.Heatmap, ax_cm.scene.plots))
        z_mat = hm_plot[3][]
        # In Makie heatmap(1:n, 1:n, z), z[col_pred, row_true]:
        # Pred=2, True=1 -> z_mat[2, 1] == 0.1; Pred=1, True=2 -> z_mat[1, 2] == 0.4
        @test z_mat[1, 1] ≈ 0.9
        @test z_mat[2, 1] ≈ 0.1
        @test z_mat[1, 2] ≈ 0.4
        @test z_mat[2, 2] ≈ 0.6
    end
end
