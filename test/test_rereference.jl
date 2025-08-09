using Test
using DataFrames
using eegfun

@testset "rereference" begin
    # Layout and data
    layout_df = DataFrame(label = [:A, :B, :M1, :M2], inc = [0.0, 0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0, 0.0])
    layout = eegfun.Layout(layout_df, nothing, nothing)
    n = 5
    t = collect(0:0.001:(n-1)*0.001)
    df = DataFrame(time = t, sample = 1:n)
    df.A = [1, 2, 3, 4, 5]
    df.B = [2, 4, 6, 8, 10]
    df.M1 = [1, 1, 1, 1, 1]
    df.M2 = [3, 3, 3, 3, 3]
    dat = eegfun.ContinuousData(copy(df, copycols=true), layout, 1000, eegfun.AnalysisInfo())

    # 1) Average reference (:avg) subtracts mean of A,B,M1,M2 from each selected channel
    dat_avg = copy(dat)
    eegfun.rereference!(dat_avg, :avg, eegfun.channels([:A, :B]))
    avg_ref = (df.A .+ df.B .+ df.M1 .+ df.M2) ./ 4
    @test all(dat_avg.data.A .== df.A .- avg_ref)
    @test all(dat_avg.data.B .== df.B .- avg_ref)
    @test eegfun.reference(dat_avg.analysis_info) == :avg

    # 2) Mastoid reference uses M1 and M2
    dat_mast = copy(dat)
    eegfun.rereference!(dat_mast, :mastoid, eegfun.channels([:A, :B]))
    mast_ref = (df.M1 .+ df.M2) ./ 2
    @test all(dat_mast.data.A .== df.A .- mast_ref)
    @test all(dat_mast.data.B .== df.B .- mast_ref)
    @test eegfun.reference(dat_mast.analysis_info) == :mastoid

    # 3) Specific reference [A] applied to B only
    dat_spec = copy(dat)
    eegfun.rereference!(dat_spec, [:A], eegfun.channels([:B]))
    @test all(dat_spec.data.B .== df.B .- df.A)
    @test eegfun.reference(dat_spec.analysis_info) == Symbol("A")

    # 4) EpochData: per-epoch reference computed independently
    df1 = copy(df); df2 = copy(df)
    df1o = copy(df1); df2o = copy(df2)
    df1.epoch = fill(1, n); df2.epoch = fill(2, n)
    ep = eegfun.EpochData([df1, df2], layout, 1000, eegfun.AnalysisInfo())
    eegfun.rereference!(ep, :mastoid, eegfun.channels([:A, :B]))
    @test all(ep.data[1].A .== df1o.A .- ((df1o.M1 .+ df1o.M2) ./ 2))
    @test all(ep.data[2].B .== df2o.B .- ((df2o.M1 .+ df2o.M2) ./ 2))

    # 5) Non-mutating version returns a new object; original unchanged
    dat_nm = eegfun.rereference(dat, :mastoid, eegfun.channels([:A]))
    @test :A ∈ propertynames(dat_nm.data) && :A ∈ propertynames(dat.data)
    @test !all(dat_nm.data.A .== dat.data.A)

    # 6) Only selected channels modified (B, M1, M2 unchanged)
    dat_sel = copy(dat)
    B0, M10, M20 = copy(df.B), copy(df.M1), copy(df.M2)
    eegfun.rereference!(dat_sel, :mastoid, eegfun.channels([:A]))
    mast_ref2 = (df.M1 .+ df.M2) ./ 2
    @test all(dat_sel.data.A .== df.A .- mast_ref2)
    @test all(dat_sel.data.B .== B0)
    @test all(dat_sel.data.M1 .== M10)
    @test all(dat_sel.data.M2 .== M20)

    # 7) Channel included in both reference and selection becomes zero
    dat_zero = copy(dat)
    eegfun.rereference!(dat_zero, [:A], eegfun.channels([:A]))
    @test all(dat_zero.data.A .== 0.0)

    # 8) ErpData mutate path
    erp_df = select(df, [:time, :A, :B, :M1, :M2])
    erp2 = eegfun.ErpData(copy(erp_df, copycols=true), layout, 1000, eegfun.AnalysisInfo(), 10)
    eegfun.rereference!(erp2, :mastoid, eegfun.channels([:A]))
    @test all(erp2.data.A .== erp_df.A .- ((erp_df.M1 .+ erp_df.M2) ./ 2))

    # 9) Missing reference channel should throw
    @test_throws Any eegfun.rereference!(copy(dat), [:Z], eegfun.channels([:A]))

    # 10) Reference info set for explicit vector
    dat_vec = copy(dat)
    eegfun.rereference!(dat_vec, [:M1, :M2], eegfun.channels([:A]))
    @test eegfun.reference(dat_vec.analysis_info) == Symbol("M1_M2")
end


