using Test
using Makie
using EegFun

@testset "plot_epochs_rejection" begin
    dat, _ = EegFun.create_test_epoch_data_with_artifacts(n_epochs = 16, n_channels = 3)
    EegFun.is_extreme_value!(dat, 100)
    auto_info = EegFun.detect_bad_epochs_automatic(dat)
    # Ensure both page 1 (epochs 1..6 for dims=(2,3)) and page 3 (epochs 13..16) have auto-flagged epochs
    push!(auto_info.rejected, EegFun.Rejection(:Ch1, 2))
    push!(auto_info.rejected, EegFun.Rejection(:Ch2, 15))

    @testset "colormap = nothing, return shape, and unvisited page seeding + manual uncheck persistence" begin
        state =
            EegFun.detect_bad_epochs_interactive(dat; dims = (2, 3), artifact_info = auto_info, colormap = nothing, display_plot = false)
        @test state isa EegFun.EpochRejectionState
        @test state.fig isa Figure
        @test state.axes isa Vector{Axis}
        @test length(state.axes) == 6
        # .epoch_axes property alias and tuple destructuring
        @test state.epoch_axes === state.axes
        fig_d, axes_d = state
        @test fig_d === state.fig
        @test axes_d === state.axes

        # Unvisited page 3 epoch (15) is included in rejected_epochs before visiting page 3
        @test 2 in EegFun.get_rejected_epochs(state)
        @test 15 in EegFun.get_rejected_epochs(state)

        # User manually unchecks epoch 2 on page 1
        state.checkboxes[2].active[] = false
        @test !(2 in EegFun.get_rejected_epochs(state))

        # Navigate to page 2, then page 3, then back to page 1
        state.current_page[] = 2
        state.current_page[] = 3
        @test 15 in EegFun.get_rejected_epochs(state)
        state.current_page[] = 1

        # Manual uncheck on epoch 2 persists across page navigation
        @test !(2 in EegFun.get_rejected_epochs(state))
        @test 15 in EegFun.get_rejected_epochs(state)

        cleaned = EegFun.reject_epochs(dat, state)
        @test EegFun.n_epochs(cleaned) == 16 - length(EegFun.get_rejected_epochs(state))
    end
end
