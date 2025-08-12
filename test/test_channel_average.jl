using Test
using DataFrames
using eegfun

@testset "channel_average" begin
    # Build a tiny fake layout with minimal polar (inc/azi)
    layout_df = DataFrame(
        label = [:A, :B, :C],
        inc = [30.0, 40.0, 50.0],
        azi = [0.0, 90.0, 180.0],
    )
    layout = eegfun.Layout(layout_df, nothing, nothing)

    # Build tiny ContinuousData with 3 channels + meta
    n = 5
    time = collect(range(0, step=0.001, length=n))
    df = DataFrame(time = time, sample = 1:n)
    df[!, :A] = 1.0 .* collect(1:n)
    df[!, :B] = 3.0 .* collect(1:n)
    df[!, :C] = 5.0 .* collect(1:n)
    dat = eegfun.ContinuousData(df, layout, 1000, eegfun.AnalysisInfo())

    # 1) Append averaged columns only (Symbols input)
    eegfun.channel_average!(dat, channel_selections = [eegfun.channels([:A, :B])])
    @test :A_B ∈ propertynames(dat.data)
    @test all(dat.data.A_B .== (df.A .+ df.B) ./ 2)

    # 1a) Predicate equals explicit groups
    dat_pred = copy(dat)
    eegfun.channel_average!(dat_pred, channel_selections = [eegfun.channels([:A, :B])])
    @test all(dat_pred.data.A_B .== dat.data.A_B)

    # 2) Reduce to only averages and create averaged layout
    dat2 = eegfun.channel_average(dat, channel_selections = [eegfun.channels([:A, :C])]; reduce = true)
    @test all(propertynames(dat2.data) .== [:time, :sample, :A_C])
    @test size(dat2.layout.data, 1) == 1
    @test :inc ∈ propertynames(dat2.layout.data)
    @test :azi ∈ propertynames(dat2.layout.data)

    # 3) Append averaged channels to layout
    dat3 = eegfun.channel_average(dat, channel_selections = [eegfun.channels([:B, :C])])
    @test :B_C ∈ propertynames(dat3.data)
    @test any(dat3.layout.data.label .== :B_C)

    # 4) Mixed Symbol input
    dat4 = copy(dat)
    eegfun.channel_average!(dat4, channel_selections = [eegfun.channels([:A, :C])])
    @test :A_C ∈ propertynames(dat4.data)

    # 5) Auto-label :avg for all channels
    dat5 = eegfun.channel_average(dat, channel_selections = [eegfun.channels([:A,:B,:C])]; reduce = true)
    @test :avg ∈ propertynames(dat5.data)

    # 6) Custom output_labels applied and length mismatch errors
    dat6 = copy(dat)
    eegfun.channel_average!(dat6, channel_selections = [eegfun.channels([:A,:B])]; output_labels = [:AB_label])
    @test :AB_label ∈ propertynames(dat6.data)
    @test_throws Any eegfun.channel_average!(copy(dat), channel_selections = [eegfun.channels([:A,:B]), eegfun.channels([:A,:C])]; output_labels = [:x])

    # 7) Duplicate labels in layout (should accumulate)
    # First add B_C, then add again to verify rows accumulate
    dat7 = eegfun.channel_average(dat, channel_selections = [eegfun.channels([:B,:C])])
    n1 = sum(dat7.layout.data.label .== :B_C)
    dat7b = eegfun.channel_average(dat7, channel_selections = [eegfun.channels([:B,:C])])
    n2 = sum(dat7b.layout.data.label .== :B_C)
    @test n1 == 1 && n2 == 2  # Should accumulate rows

    # 8) Non-mutating variant: original unchanged (use a fresh base)
    df_base = DataFrame(time = time, sample = 1:n)
    df_base[!, :A] = 1.0 .* collect(1:n)
    df_base[!, :B] = 3.0 .* collect(1:n)
    df_base[!, :C] = 5.0 .* collect(1:n)
    dat8 = eegfun.ContinuousData(copy(df_base, copycols=true), layout, 1000, eegfun.AnalysisInfo())
    dat8b = eegfun.channel_average(dat8, channel_selections = [eegfun.channels([:A,:B])])
    @test :A_B ∈ propertynames(dat8b.data)
    @test :A_B ∉ propertynames(dat8.data)

    # 9) EpochData append and reduce
    # Create simple EpochData (2 epochs) from dat
    df1 = copy(df); df2 = copy(df)
    df1.epoch = fill(1, n); df2.epoch = fill(2, n)
    ep = eegfun.EpochData([df1, df2], layout, 1000, eegfun.AnalysisInfo())
    ep_app = eegfun.channel_average(ep, channel_selections = [eegfun.channels([:A,:B])])
    @test :A_B ∈ propertynames(ep_app.data[1]) && :A_B ∈ propertynames(ep_app.data[2])
    ep_red = eegfun.channel_average(ep, channel_selections = [eegfun.channels([:A,:B])]; reduce = true)
    # Expect meta columns (leading) + A_B only
    cols_ep1 = propertynames(ep_red.data[1])
    @test cols_ep1[end] == :A_B
    @test :A ∉ cols_ep1 && :B ∉ cols_ep1 && :C ∉ cols_ep1

    # 10) ErpData reduce path
    # Build a tiny ErpData from dat by averaging across time (re-using ContinuousData schema)
    erp_df = select(df, [:time, :A, :B, :C])
    erp = eegfun.ErpData(erp_df, layout, 1000, eegfun.AnalysisInfo(), 10)
    erp_red = eegfun.channel_average(erp, channel_selections = [eegfun.channels([:A,:B])]; reduce = true)
    @test all(propertynames(erp_red.data) .== [:time, :A_B])

    # 11) Test default behavior (average all channels)
    dat_default = copy(dat)
    eegfun.channel_average!(dat_default)  # Should use default channels() from second function
    @test :avg ∈ propertynames(dat_default.data)
    @test all(dat_default.data.avg .== (df.A .+ df.B .+ df.C) ./ 3)

    # 12) Test single channel selection with custom label
    dat_single = copy(dat)
    eegfun.channel_average!(dat_single, channel_selections = [eegfun.channels([:A, :B])], output_labels = [:AB_custom])
    @test :AB_custom ∈ propertynames(dat_single.data)
    @test all(dat_single.data.AB_custom .== (df.A .+ df.B) ./ 2)

    # 13) Test empty channel selections (should do nothing)
    dat_empty = copy(dat)
    original_cols = propertynames(dat_empty.data)
    eegfun.channel_average!(dat_empty, channel_selections = [])
    @test propertynames(dat_empty.data) == original_cols  # No new columns added
end


