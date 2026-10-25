using Test
using EegFun
using LinearAlgebra

@testset "Current Source Density (CSD)" begin
    # 1. ContinuousData
    c_dat = EegFun.create_test_continuous_data(n = 500, fs = 256, n_channels = 4, seed = 1)
    c_orig = copy(c_dat)

    # Out-of-place
    csd_c = EegFun.compute_csd(c_dat; m = 4, lambda = 1e-5)
    @test csd_c isa EegFun.ContinuousData
    # Input not modified
    @test c_dat.data == c_orig.data
    # Output has same dimensions
    @test size(csd_c.data) == size(c_orig.data)
    # Output data modified by spatial filter
    @test csd_c.data != c_orig.data
    @test all(isfinite.(Matrix(csd_c.data[:, EegFun.channel_labels(csd_c)])))

    # In-place
    c_mut = copy(c_dat)
    res_c = EegFun.compute_csd!(c_mut; m = 4, lambda = 1e-5)
    @test res_c === c_mut
    @test c_mut.data == csd_c.data

    # 2. EpochData
    ep_dat = EegFun.create_test_epoch_data(n = 256, fs = 256, n_epochs = 4, n_channels = 4, seed = 1)
    ep_orig = copy(ep_dat)

    csd_ep = EegFun.compute_csd(ep_dat)
    @test csd_ep isa EegFun.EpochData
    @test length(csd_ep.data) == 4
    @test ep_dat.data[1] == ep_orig.data[1]
    @test csd_ep.data[1] != ep_orig.data[1]
    @test size(csd_ep.data[1]) == size(ep_orig.data[1])

    ep_mut = copy(ep_dat)
    EegFun.compute_csd!(ep_mut)
    @test ep_mut.data[1] == csd_ep.data[1]

    # 3. ErpData
    erp_dat = EegFun.create_test_erp_data(n_channels = 4, seed = 1)
    erp_orig = copy(erp_dat)

    csd_erp = EegFun.compute_csd(erp_dat)
    @test csd_erp isa EegFun.ErpData
    @test erp_dat.data == erp_orig.data
    @test csd_erp.data != erp_orig.data
    @test size(csd_erp.data) == size(erp_orig.data)

    erp_mut = copy(erp_dat)
    EegFun.compute_csd!(erp_mut)
    @test erp_mut.data == csd_erp.data

    # Vector dispatch
    erp_vec = [erp_dat, copy(erp_dat)]
    csd_vec = EegFun.compute_csd(erp_vec)
    @test length(csd_vec) == 2
    @test csd_vec[1].data == csd_erp.data

    erp_vec_mut = [copy(erp_dat), copy(erp_dat)]
    EegFun.compute_csd!(erp_vec_mut)
    @test erp_vec_mut[1].data == csd_erp.data
    @test erp_vec_mut[2].data == csd_erp.data

    # 4. Legendre evaluation helper
    cosang = [1.0 0.5; 0.5 1.0]
    H = EegFun._calc_h(cosang, 4, 50)
    @test size(H) == (2, 2)
    @test all(isfinite, H)
    @test H[1, 1] == H[2, 2] # Symmetry on diagonal
    @test H[1, 2] == H[2, 1] # Symmetry off diagonal
end
