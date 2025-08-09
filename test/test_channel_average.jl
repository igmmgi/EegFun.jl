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

    # 1) Append averaged columns only (Strings input)
    eegfun.channel_average!(dat, [["A", "B"]])
    @test :A_B ∈ propertynames(dat.data)
    @test all(dat.data.A_B .== (df.A .+ df.B) ./ 2)

    # 1a) Predicate equals explicit groups
    dat_pred = copy(dat)
    eegfun.channel_average!(dat_pred, [eegfun.channels([:A, :B])])
    @test all(dat_pred.data.A_B .== dat.data.A_B)

    # 2) Reduce to only averages and create averaged layout
    dat2 = eegfun.channel_average(dat, [["A", "C"]]; reduce = true)
    @test all(propertynames(dat2.data) .== [:time, :sample, :A_C])
    @test size(dat2.layout.data, 1) == 1
    @test :inc ∈ propertynames(dat2.layout.data)
    @test :azi ∈ propertynames(dat2.layout.data)

    # 3) Update layout with an additional averaged label
    dat3 = eegfun.channel_average(dat, [["B", "C"]]; update_layout = true)
    @test :B_C ∈ propertynames(dat3.data)
    @test any(dat3.layout.data.label .== :B_C)

    # 4) Mixed Symbol/String input
    dat4 = copy(dat)
    eegfun.channel_average!(dat4, [["A", :C]])
    @test :A_C ∈ propertynames(dat4.data)

    # 5) Auto-label :avg for all channels
    dat5 = eegfun.channel_average(dat, [["A","B","C"]]; reduce = true)
    @test :avg ∈ propertynames(dat5.data)

    # 6) Custom output_labels applied and length mismatch errors
    dat6 = copy(dat)
    eegfun.channel_average!(dat6, [["A","B"]]; output_labels = [:AB_label])
    @test :AB_label ∈ propertynames(dat6.data)
    @test_throws Any eegfun.channel_average!(copy(dat), [["A","B"],["A","C"]]; output_labels = [:x])

    # 7) Duplicate label in layout replaced, not duplicated
    # First add B_C, then add again with update_layout to verify single row remains
    dat7 = eegfun.channel_average(dat, [["B","C"]]; update_layout = true)
    n1 = sum(dat7.layout.data.label .== :B_C)
    dat7b = eegfun.channel_average(dat7, [["B","C"]]; update_layout = true)
    n2 = sum(dat7b.layout.data.label .== :B_C)
    @test n1 == 1 && n2 == 1

    # 8) Non-mutating variant: original unchanged (use a fresh base)
    df_base = DataFrame(time = time, sample = 1:n)
    df_base[!, :A] = 1.0 .* collect(1:n)
    df_base[!, :B] = 3.0 .* collect(1:n)
    df_base[!, :C] = 5.0 .* collect(1:n)
    dat8 = eegfun.ContinuousData(copy(df_base, copycols=true), layout, 1000, eegfun.AnalysisInfo())
    dat8b = eegfun.channel_average(dat8, [["A","B"]])
    @test :A_B ∈ propertynames(dat8b.data)
    @test :A_B ∉ propertynames(dat8.data)

    # 9) EpochData append and reduce
    # Create simple EpochData (2 epochs) from dat
    df1 = copy(df); df2 = copy(df)
    df1.epoch = fill(1, n); df2.epoch = fill(2, n)
    ep = eegfun.EpochData([df1, df2], layout, 1000, eegfun.AnalysisInfo())
    ep_app = eegfun.channel_average(ep, [["A","B"]])
    @test :A_B ∈ propertynames(ep_app.data[1]) && :A_B ∈ propertynames(ep_app.data[2])
    ep_red = eegfun.channel_average(ep, [["A","B"]]; reduce = true)
    # Expect meta columns (leading) + A_B only
    cols_ep1 = propertynames(ep_red.data[1])
    @test cols_ep1[end] == :A_B
    @test :A ∉ cols_ep1 && :B ∉ cols_ep1 && :C ∉ cols_ep1

    # 10) ErpData reduce path
    # Build a tiny ErpData from dat by averaging across time (re-using ContinuousData schema)
    erp_df = select(df, [:time, :A, :B, :C])
    erp = eegfun.ErpData(erp_df, layout, 1000, eegfun.AnalysisInfo(), 10)
    erp_red = eegfun.channel_average(erp, [["A","B"]]; reduce = true)
    @test all(propertynames(erp_red.data) .== [:time, :A_B])
end


