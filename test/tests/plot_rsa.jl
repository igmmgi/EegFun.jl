using Test
using Makie
using EegFun

@testset "plot_rsa" begin
    dat = EegFun.create_test_epoch_data_vector(conditions = 1:3, n_epochs = 5, n_channels = 3)

    model = EegFun.create_rdm_from_categorical([1, 1, 2])

    res = EegFun.rsa(dat)
    res = EegFun.compare_models(res, [model])

    @testset "basic functionality" begin
        result = EegFun.plot_rsa(res; display_plot = false)
        @test result isa Figure
    end
end
