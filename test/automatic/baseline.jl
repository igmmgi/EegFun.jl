using Test
using DataFrames
using Statistics
using eegfun

@testset "baseline" begin

    # 1) Baseline over first sample for Ch1, Ch2, and Ch3
    dat = create_test_data(n=6)
    eegfun.baseline!(dat, eegfun.IntervalIdx(1, 1))

    @test isapprox(mean(dat.data.Ch1[1]), 0.0; atol = 1e-9)  
    @test isapprox(mean(dat.data.Ch2[1]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch3[1]), 0.0; atol = 1e-9)

    dat = create_test_data(n=6)
    dat = eegfun.baseline(dat, eegfun.IntervalIdx(1, 1))

    @test isapprox(mean(dat.data.Ch1[1]), 0.0; atol = 1e-9) 
    @test isapprox(mean(dat.data.Ch2[1]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch3[1]), 0.0; atol = 1e-9)

    # 2) Baseline over entire range for all channels
    dat = create_test_data(n=6)
    eegfun.baseline!(dat)

    @test isapprox(mean(dat.data.Ch1), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch2), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch3), 0.0; atol = 1e-9)

    dat = create_test_data(n=6)
    dat = eegfun.baseline(dat)

    @test isapprox(mean(dat.data.Ch1), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch2), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch3), 0.0; atol = 1e-9)


    # 3) EpochData: each epoch baselined independently
    dat = create_test_epoch_data(n=3)
    eegfun.baseline!(dat, eegfun.IntervalIdx(1, 1))

    @test isapprox(mean(dat.data[1].Ch1[1]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data[2].Ch2[1]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data[3].Ch3[1]), 0.0; atol = 1e-9)

    dat = create_test_epoch_data(n=3)
    eegfun.baseline!(dat, eegfun.IntervalIdx(1, 3))
    
    @test isapprox(mean(dat.data[1].Ch1[1:3]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data[2].Ch2[1:3]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data[3].Ch3[1:3]), 0.0; atol = 1e-9)

    # 4) IntervalTime converted to indices correctly
    dat = create_test_data(n=6)
    eegfun.baseline!(dat, eegfun.IntervalTime(0.0, 0.0))

    @test isapprox(mean(dat.data.Ch1[1]), 0.0; atol = 1e-9)  
    @test isapprox(mean(dat.data.Ch2[1]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch3[1]), 0.0; atol = 1e-9)

    # 4) IntervalTime converted to indices correctly
    dat = create_test_data(n=6)
    dat = eegfun.baseline(dat, eegfun.IntervalTime(0.0, 0.0))

    @test isapprox(mean(dat.data.Ch1[1]), 0.0; atol = 1e-9)  
    @test isapprox(mean(dat.data.Ch2[1]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch3[1]), 0.0; atol = 1e-9)


    # 4) IntervalTime converted to indices correctly
    dat = create_test_data(n=6)
    eegfun.baseline!(dat, eegfun.IntervalTime(0.003, 0.003))

    @test isapprox(mean(dat.data.Ch1[4]), 0.0; atol = 1e-9)  
    @test isapprox(mean(dat.data.Ch2[4]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch3[4]), 0.0; atol = 1e-9)

    dat = create_test_data(n=6)
    dat = eegfun.baseline(dat, eegfun.IntervalTime(0.003, 0.003))

    @test isapprox(mean(dat.data.Ch1[4]), 0.0; atol = 1e-9)  
    @test isapprox(mean(dat.data.Ch2[4]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data.Ch3[4]), 0.0; atol = 1e-9)


    dat = create_test_epoch_data(n=6)
    eegfun.baseline!(dat, eegfun.IntervalTime(0.003, 0.003))

    @test isapprox(mean(dat.data[1].Ch1[4]), 0.0; atol = 1e-9)  
    @test isapprox(mean(dat.data[2].Ch2[4]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data[3].Ch3[4]), 0.0; atol = 1e-9)

    dat = create_test_epoch_data(n=6)
    dat = eegfun.baseline(dat, eegfun.IntervalTime(0.003, 0.003))

    @test isapprox(mean(dat.data[1].Ch1[4]), 0.0; atol = 1e-9) 
    @test isapprox(mean(dat.data[2].Ch2[4]), 0.0; atol = 1e-9)
    @test isapprox(mean(dat.data[3].Ch3[4]), 0.0; atol = 1e-9)

end
