using Test
using DataFrames
using Statistics
using eegfun

@testset "baseline" begin
    # Build layout and ContinuousData
    layout_df = DataFrame(label = [:A, :B, :C], inc = [0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0])
    layout = eegfun.Layout(layout_df, nothing, nothing)
    n = 6
    t = collect(0:0.001:(n-1)*0.001)
    df = DataFrame(time = t, sample = 1:n)
    # A increases, B constant offset, C decreasing
    df.A = [1, 2, 3, 4, 5, 6]
    df.B = [10, 10, 10, 10, 10, 10]
    df.C = [6, 5, 4, 3, 2, 1]
    dat = eegfun.ContinuousData(copy(df, copycols=true), layout, 1000, eegfun.AnalysisInfo())

    # 1) Baseline over first 3 samples for A and B only
    dat1 = copy(dat)
    eegfun.baseline!(dat1, eegfun.IntervalIdx(1, 3); channel_selection = eegfun.channels([:A, :B]))
    @test !isapprox(mean(df.A[1:3]), 0.0; atol=1e-9)  # sanity on original
    @test isapprox(mean(dat1.data.A[1:3]), 0.0; atol=1e-9)
    @test isapprox(mean(dat1.data.B[1:3]), 0.0; atol=1e-9)
    @test all(dat1.data.C .== df.C)  # C unchanged

    # 2) Baseline over entire range for all channels
    dat2 = copy(dat)
    eegfun.baseline!(dat2)
    @test isapprox(mean(dat2.data.A), 0.0; atol=1e-9)
    @test isapprox(mean(dat2.data.B), 0.0; atol=1e-9)
    @test isapprox(mean(dat2.data.C), 0.0; atol=1e-9)

    # 3) EpochData: each epoch baselined independently
    df1 = copy(df); df2 = copy(df)
    df1.epoch = fill(1, n); df2.epoch = fill(2, n)
    ep = eegfun.EpochData([df1, df2], layout, 1000, eegfun.AnalysisInfo())
    eegfun.baseline!(ep, eegfun.IntervalIdx(1, 3); channel_selection = eegfun.channels([:A, :B, :C]))
    @test isapprox(mean(ep.data[1].A[1:3]), 0.0; atol=1e-9)
    @test isapprox(mean(ep.data[2].B[1:3]), 0.0; atol=1e-9)

    # 4) IntervalTime converted to indices correctly
    dat3 = copy(dat)
    eegfun.baseline!(dat3, eegfun.IntervalTime(0.0, 0.002); channel_selection = eegfun.channels([:A]))
    @test isapprox(mean(dat3.data.A[1:3]), 0.0; atol=1e-9)

    # 5) Non-mutating: original unchanged
    dat4 = eegfun.baseline(dat, eegfun.IntervalIdx(1, 2); channel_selection = eegfun.channels([:A]))
    @test !isapprox(mean(dat.data.A[1:2]), 0.0; atol=1e-9)
    @test isapprox(mean(dat4.data.A[1:2]), 0.0; atol=1e-9)
end


