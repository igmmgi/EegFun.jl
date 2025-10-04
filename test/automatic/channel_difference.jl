using Test
using DataFrames
using eegfun

@testset "channel_difference" begin
    # Layout with three channels
    layout_df = DataFrame(label = [:A, :B, :C], inc = [30.0, 40.0, 50.0], azi = [0.0, 90.0, 180.0])
    layout = eegfun.Layout(layout_df, nothing, nothing)

    # Build ContinuousData with known values
    n = 5
    time = collect(range(0, step = 0.001, length = n))
    df = DataFrame(time = time, sample = 1:n)
    df[!, :A] = collect(1.0:n)
    df[!, :B] = 2 .* collect(1.0:n)
    df[!, :C] = 4 .* collect(1.0:n)
    dat = eegfun.ContinuousData(copy(df, copycols = true), layout, 1000, eegfun.AnalysisInfo())

    # 1) Basic difference A - B
    eegfun.channel_difference!(
        dat;
        channel_selection1 = eegfun.channels([:A]),
        channel_selection2 = eegfun.channels([:B]),
        channel_out = :A_minus_B,
    )
    @test :A_minus_B ∈ propertynames(dat.data)
    @test all(dat.data.A_minus_B .== (df.A .- df.B))

    # 2) Group average difference mean(A,B) - C
    eegfun.channel_difference!(
        dat;
        channel_selection1 = eegfun.channels([:A, :B]),
        channel_selection2 = eegfun.channels([:C]),
        channel_out = :AB_minus_C,
    )
    @test :AB_minus_C ∈ propertynames(dat.data)
    @test all(dat.data.AB_minus_C .== ((df.A .+ df.B) ./ 2 .- df.C))

    # 3) Predicate selection (same as explicit)
    dat_pred = copy(dat)
    eegfun.channel_difference!(
        dat_pred;
        channel_selection1 = eegfun.channels([:A, :B]),
        channel_selection2 = eegfun.channels([:C]),
        channel_out = :P_AB_minus_C,
    )
    @test all(dat_pred.data.P_AB_minus_C .== dat.data.AB_minus_C)

    # 4) Non-mutating version returns a new object; original unchanged
    dat_nm = eegfun.channel_difference(
        dat;
        channel_selection1 = eegfun.channels([:A]),
        channel_selection2 = eegfun.channels([:C]),
        channel_out = :A_minus_C,
    )
    @test :A_minus_C ∈ propertynames(dat_nm.data)
    @test :A_minus_C ∉ propertynames(dat.data)
    @test all(dat_nm.data.A_minus_C .== (df.A .- df.C))

    # 5) Overwrite behavior: write then overwrite with different groups
    dat_ow = copy(dat)
    dat_ow.data[!, :X] = zeros(n)
    eegfun.channel_difference!(
        dat_ow;
        channel_selection1 = eegfun.channels([:B]),
        channel_selection2 = eegfun.channels([:A]),
        channel_out = :X,
    )
    @test all(dat_ow.data.X .== (df.B .- df.A))

    # 6) EpochData: append to each epoch
    df1 = copy(df);
    df2 = copy(df)
    df1.epoch = fill(1, n);
    df2.epoch = fill(2, n)
    ep = eegfun.EpochData([df1, df2], layout, 1000, eegfun.AnalysisInfo())
    eegfun.channel_difference!(
        ep;
        channel_selection1 = eegfun.channels([:A]),
        channel_selection2 = eegfun.channels([:B]),
        channel_out = :A_minus_B,
    )
    @test :A_minus_B ∈ propertynames(ep.data[1]) && :A_minus_B ∈ propertynames(ep.data[2])
    @test all(ep.data[1].A_minus_B .== (df1.A .- df1.B))
    @test all(ep.data[2].A_minus_B .== (df2.A .- df2.B))

    # 7) ErpData (SingleDataFrameEeg): append
    erp_df = select(df, [:time, :A, :B, :C])
    erp = eegfun.ErpData(copy(erp_df, copycols = true), layout, 1000, eegfun.AnalysisInfo(), 10)
    eegfun.channel_difference!(
        erp;
        channel_selection1 = eegfun.channels([:C]),
        channel_selection2 = eegfun.channels([:B]),
        channel_out = :C_minus_B,
    )
    @test :C_minus_B ∈ propertynames(erp.data)
    @test all(erp.data.C_minus_B .== (erp_df.C .- erp_df.B))

    # 8) Missing/nonexistent in group1 resolves to empty selection -> result is NaN vector
    dat_miss = copy(dat)
    eegfun.channel_difference!(
        dat_miss;
        channel_selection1 = eegfun.channels([:Z]),
        channel_selection2 = eegfun.channels([:A]),
        channel_out = :Z_minus_A,
    )
    @test :Z_minus_A ∈ propertynames(dat_miss.data)
    @test all(isnan, dat_miss.data.Z_minus_A)

    # 9) Empty selection (both groups empty) -> NaN vector (0/0)
    dat_empty = copy(dat)
    eegfun.channel_difference!(
        dat_empty;
        channel_selection1 = x -> falses(length(x)),
        channel_selection2 = x -> falses(length(x)),
        channel_out = :none,
    )
    @test :none ∈ propertynames(dat_empty.data)
    @test all(isnan, dat_empty.data.none)

    # 10) Commutativity sanity: A-B == -(B-A)
    dat_ab = copy(dat)
    eegfun.channel_difference!(
        dat_ab;
        channel_selection1 = eegfun.channels([:A]),
        channel_selection2 = eegfun.channels([:B]),
        channel_out = :AB,
    )
    dat_ba = copy(dat)
    eegfun.channel_difference!(
        dat_ba;
        channel_selection1 = eegfun.channels([:B]),
        channel_selection2 = eegfun.channels([:A]),
        channel_out = :BA,
    )
    @test all(dat_ab.data.AB .== .-(dat_ba.data.BA))

    # 11) Layout invariance: no change to layout size or labels
    nrows_before = nrow(dat.layout.data)
    labels_before = copy(dat.layout.data.label)
    dat_layout = copy(dat)
    eegfun.channel_difference!(
        dat_layout;
        channel_selection1 = eegfun.channels([:A]),
        channel_selection2 = eegfun.channels([:B]),
        channel_out = :tmp,
    )
    @test nrow(dat_layout.layout.data) == nrows_before
    @test all(dat_layout.layout.data.label .== labels_before)

    # 12) Non-mutating invariants
    dat_orig = copy(dat)
    dat_new = eegfun.channel_difference(
        dat_orig;
        channel_selection1 = eegfun.channels([:A]),
        channel_selection2 = eegfun.channels([:B]),
        channel_out = :A_minus_B2,
    )
    @test dat_new.sample_rate == dat_orig.sample_rate
    @test eegfun.reference(dat_new.analysis_info) == eegfun.reference(dat_orig.analysis_info)
    @test :A_minus_B2 ∉ propertynames(dat_orig.data)
end
