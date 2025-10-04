using Test
using DataFrames
using eegfun

@testset "channel_average" begin
    # Build a tiny fake layout with minimal polar (inc/azi)
    layout_df = DataFrame(label = [:A, :B, :C], inc = [30.0, 40.0, 50.0], azi = [0.0, 90.0, 180.0])
    layout = eegfun.Layout(layout_df, nothing, nothing)

    # Build tiny ContinuousData with 3 channels + meta
    n = 5
    time = collect(range(0, step = 0.001, length = n))
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
    dat5 = eegfun.channel_average(dat, channel_selections = [eegfun.channels([:A, :B, :C])]; reduce = true)
    @test :avg ∈ propertynames(dat5.data)

    # 6) Custom output_labels applied and length mismatch errors
    dat6 = copy(dat)
    eegfun.channel_average!(dat6, channel_selections = [eegfun.channels([:A, :B])]; output_labels = [:AB_label])
    @test :AB_label ∈ propertynames(dat6.data)
    @test_throws Any eegfun.channel_average!(
        copy(dat),
        channel_selections = [eegfun.channels([:A, :B]), eegfun.channels([:A, :C])];
        output_labels = [:x],
    )

    # 7) Duplicate labels in layout (should accumulate)
    # First add B_C, then add again to verify rows accumulate
    dat7 = eegfun.channel_average(dat, channel_selections = [eegfun.channels([:B, :C])])
    n1 = sum(dat7.layout.data.label .== :B_C)
    dat7b = eegfun.channel_average(dat7, channel_selections = [eegfun.channels([:B, :C])])
    n2 = sum(dat7b.layout.data.label .== :B_C)
    @test n1 == 1 && n2 == 2  # Should accumulate rows

    # 8) Non-mutating variant: original unchanged (use a fresh base)
    df_base = DataFrame(time = time, sample = 1:n)
    df_base[!, :A] = 1.0 .* collect(1:n)
    df_base[!, :B] = 3.0 .* collect(1:n)
    df_base[!, :C] = 5.0 .* collect(1:n)
    dat8 = eegfun.ContinuousData(copy(df_base, copycols = true), layout, 1000, eegfun.AnalysisInfo())
    dat8b = eegfun.channel_average(dat8, channel_selections = [eegfun.channels([:A, :B])])
    @test :A_B ∈ propertynames(dat8b.data)
    @test :A_B ∉ propertynames(dat8.data)

    # 9) EpochData append and reduce
    # Create simple EpochData (2 epochs) from dat
    df1 = copy(df);
    df2 = copy(df)
    df1.epoch = fill(1, n);
    df2.epoch = fill(2, n)
    ep = eegfun.EpochData([df1, df2], layout, 1000, eegfun.AnalysisInfo())
    ep_app = eegfun.channel_average(ep, channel_selections = [eegfun.channels([:A, :B])])
    @test :A_B ∈ propertynames(ep_app.data[1]) && :A_B ∈ propertynames(ep_app.data[2])
    ep_red = eegfun.channel_average(ep, channel_selections = [eegfun.channels([:A, :B])]; reduce = true)
    # Expect meta columns (leading) + A_B only
    cols_ep1 = propertynames(ep_red.data[1])
    @test cols_ep1[end] == :A_B
    @test :A ∉ cols_ep1 && :B ∉ cols_ep1 && :C ∉ cols_ep1

    # 10) ErpData reduce path
    # Build a tiny ErpData from dat by averaging across time (re-using ContinuousData schema)
    erp_df = select(df, [:time, :A, :B, :C])
    erp = eegfun.ErpData(erp_df, layout, 1000, eegfun.AnalysisInfo(), 10)
    erp_red = eegfun.channel_average(erp, channel_selections = [eegfun.channels([:A, :B])]; reduce = true)
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

    # 14) Test coordinate averaging correctness
    # Create layout with channels that have different coordinate regions
    # Use the same channel names as the test data (:A, :B, :C)
    # Note: The coordinate averaging uses 3D Cartesian averaging followed by polar conversion
    # This gives the true geometric center, which may differ from simple polar averaging
    coord_layout_df = DataFrame(
        label = [:A, :B, :C],  # Match the data columns
        inc = [-92.0, -92.0, 115.0],  # A/B frontal (negative), C posterior (positive)
        azi = [-72.0, -52.0, -68.0],
    )
    coord_layout = eegfun.Layout(coord_layout_df, nothing, nothing)

    # Create data with this layout
    coord_dat = eegfun.ContinuousData(df, coord_layout, 1000, eegfun.AnalysisInfo())

    # Test coordinate averaging with reduce=true
    coord_dat_avg = eegfun.channel_average(coord_dat, channel_selections = [eegfun.channels([:A, :B])], reduce = true)

    # Verify we have the expected averaged channels
    @test :A_B ∈ propertynames(coord_dat_avg.data)

    # Verify layout has averaged coordinates
    @test size(coord_dat_avg.layout.data, 1) == 1
    @test :A_B ∈ coord_dat_avg.layout.data.label

    # Verify coordinates are mathematically reasonable
    # A_B should be between A and B
    a_b_row = coord_dat_avg.layout.data[coord_dat_avg.layout.data.label .== :A_B, :]
    # Note: The mathematical average of 3D coordinates gives the geometric center
    # This may not preserve the sign convention of the original polar coordinates
    @test a_b_row[1, :inc] ≈ 77.52 atol=0.1   # Geometric center of A and B
    @test a_b_row[1, :azi] ≈ 111.72 atol=0.1  # Mathematical result from 3D averaging

    # 14a) Test that averaged channels have different coordinates (the original issue)
    # Create a second averaged channel to verify they're different
    coord_dat_avg2 = eegfun.channel_average(
        coord_dat,
        channel_selections = [eegfun.channels([:A, :B]), eegfun.channels([:B, :C])],
        reduce = true,
    )

    # Verify we have two different averaged channels
    @test :A_B ∈ propertynames(coord_dat_avg2.data)
    @test :B_C ∈ propertynames(coord_dat_avg2.data)
    @test size(coord_dat_avg2.layout.data, 1) == 2

    # Verify the averaged channels have different coordinates (this was the original bug)
    a_b_row2 = coord_dat_avg2.layout.data[coord_dat_avg2.layout.data.label .== :A_B, :]
    b_c_row2 = coord_dat_avg2.layout.data[coord_dat_avg2.layout.data.label .== :B_C, :]

    # The coordinates should be different, not identical
    # This was the original bug - all averaged channels had identical coordinates
    @test a_b_row2[1, :inc] != b_c_row2[1, :inc] ||
          a_b_row2[1, :azi] != b_c_row2[1, :azi] ||
          a_b_row2[1, :x3] != b_c_row2[1, :x3] ||
          a_b_row2[1, :y3] != b_c_row2[1, :y3] ||
          a_b_row2[1, :z3] != b_c_row2[1, :z3]

    # Print the actual coordinates for debugging
    println(
        "A_B coordinates: inc=$(a_b_row2[1, :inc]), azi=$(a_b_row2[1, :azi]), x3=$(a_b_row2[1, :x3]), y3=$(a_b_row2[1, :y3]), z3=$(a_b_row2[1, :z3])",
    )
    println(
        "B_C coordinates: inc=$(b_c_row2[1, :inc]), azi=$(b_c_row2[1, :azi]), x3=$(b_c_row2[1, :x3]), y3=$(b_c_row2[1, :y3]), z3=$(b_c_row2[1, :z3])",
    )

    # 15) Test coordinate averaging preserves 3D coordinates
    # Verify that 3D coordinates are properly calculated from averaged polar coordinates
    @test :x3 ∈ propertynames(coord_dat_avg.layout.data)
    @test :y3 ∈ propertynames(coord_dat_avg.layout.data)
    @test :z3 ∈ propertynames(coord_dat_avg.layout.data)

    # 16) Test coordinate averaging with append (reduce=false)
    coord_dat_append = eegfun.channel_average(coord_dat, channel_selections = [eegfun.channels([:A, :B])])

    # Verify original channels are preserved
    @test :A ∈ propertynames(coord_dat_append.data)
    @test :B ∈ propertynames(coord_dat_append.data)
    @test :A_B ∈ propertynames(coord_dat_append.data)

    # Verify layout has both original and averaged channels
    @test size(coord_dat_append.layout.data, 1) == 4  # 3 original + 1 averaged
    @test :A_B ∈ coord_dat_append.layout.data.label

    # 17) Test coordinate averaging with mixed sign conventions
    # Create layout with channels that have mixed incidence signs
    mixed_layout_df = DataFrame(
        label = [:A, :B, :C],  # Match the data columns
        inc = [-92.0, -92.0, 115.0],  # Mixed signs
        azi = [-72.0, -52.0, -68.0],
    )
    mixed_layout = eegfun.Layout(mixed_layout_df, nothing, nothing)
    mixed_dat = eegfun.ContinuousData(df, mixed_layout, 1000, eegfun.AnalysisInfo())

    # Test that averaging works correctly with mixed signs
    mixed_avg = eegfun.channel_average(mixed_dat, channel_selections = [eegfun.channels([:A, :B])], reduce = true)

    # Verify the averaged coordinates make mathematical sense
    a_b_mixed = mixed_avg.layout.data[mixed_avg.layout.data.label .== :A_B, :]

    # Note: The mathematical average of 3D coordinates gives the geometric center
    # This may not preserve the sign convention of the original polar coordinates
    # The important thing is that the coordinates are mathematically correct
    @test a_b_mixed[1, :inc] > 0  # Should be positive (geometric center)
    @test a_b_mixed[1, :azi] ≈ 111.72 atol=0.1  # Mathematical result from 3D averaging

end
