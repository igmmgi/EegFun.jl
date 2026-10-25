using Test
using DataFrames
using Random
using Statistics
using LinearAlgebra

@testset "ica" begin

    # Helper: create synthetic continuous data with EOG and masks
    function create_synthetic_continuous(; n::Int = 2000, fs::Int = 500, nch::Int = 5)
        Random.seed!(1234)
        t = collect(0:(n-1)) ./ fs
        # Base sources
        s1 = sin.(2π .* 8 .* t)
        s2 = 0.7 .* sin.(2π .* 12 .* t .+ 0.3)
        s3 = randn(n)
        s4 = 0.2 .* sin.(2π .* 50 .* t)  # line frequency present
        s5 = 0.5 .* randn(n)
        S = hcat(s1, s2, s3, s4, s5)
        # Mixing
        mixing = randn(nch, size(S, 2))
        X = S * mixing'  # n x nch
        # Channels
        cols = Dict{Symbol,Vector{Float64}}()
        for i = 1:nch
            cols[Symbol("ch", i)] = X[:, i]
        end
        # EOG channels roughly correlated to first two sources
        vEOG = s1 .+ 0.1 .* randn(n)
        hEOG = s2 .+ 0.1 .* randn(n)
        df = DataFrame(time = t, sample = 1:n, trigger = zeros(Int, n))
        for (k, v) in cols
            df[!, k] = v
        end
        df[!, :vEOG] = vEOG
        df[!, :hEOG] = hEOG
        # Keep mask: first half only
        df[!, :keepmask] = [i <= n ÷ 2 for i = 1:n]
        layout =
            EegFun.Layout(DataFrame(label = Symbol.(collect(keys(cols))), inc = zeros(nch), azi = zeros(nch)), nothing, nothing, nothing)
        dat = EegFun.ContinuousData("test_data", copy(df, copycols = true), layout, fs, EegFun.AnalysisInfo())
        return dat
    end

    @testset "run_ica basic and parameters" begin
        dat = create_synthetic_continuous()
        ica_res = EegFun.run_ica(dat)
        @test ica_res isa EegFun.InfoIca
        @test size(ica_res.unmixing, 1) == EegFun.n_channels(dat) - 1
        @test length(ica_res.ica_label) == size(ica_res.unmixing, 1)
        @test isapprox(sum(ica_res.variance), 1.0; atol = 1e-6)

        # Custom number of components
        ica_res2 = EegFun.run_ica(dat; n_components = 3)
        @test size(ica_res2.unmixing, 1) == 3
        # n_components larger than channels -> adjusted to channels-1
        ica_res2b = EegFun.run_ica(dat; n_components = 999)
        @test size(ica_res2b.unmixing, 1) == EegFun.n_channels(dat) - 1

        # Channel selection subset
        ica_res3 = EegFun.run_ica(dat; channel_selection = EegFun.channels([:ch1, :ch2, :ch3]))
        @test size(ica_res3.unmixing, 1) == 2

        # Sample selection subset
        ica_res4 = EegFun.run_ica(dat; sample_selection = EegFun.samples(:keepmask))
        @test ica_res4 isa EegFun.InfoIca

        # Preprocessing should be external; run_ica has no filter flags

        # Include extra channels (vEOG, hEOG) in ICA
        nch = EegFun.n_channels(dat)
        sel = vcat(dat.layout.data.label, [:vEOG, :hEOG])
        ica_extra = EegFun.run_ica(dat; include_extra = true, channel_selection = EegFun.channels(sel))
        @test size(ica_extra.unmixing, 1) == (nch + 2 - 1)

        # Error when no channels selected
        @test_throws Any EegFun.run_ica(dat; channel_selection = EegFun.channels(Symbol[]))
        # Error when no samples selected
        @test_throws Any EegFun.run_ica(dat; sample_selection = x -> falses(nrow(x)))

    end

    @testset "_create_ica_data_matrix" begin
        dat = create_synthetic_continuous(; n = 200, fs = 200, nch = 4)
        # Duplicate samples to test unique by :sample
        df = vcat(dat.data, dat.data)
        mat = EegFun._create_ica_data_matrix(df, [:ch1, :ch2, :ch3, :ch4], 1:nrow(df))
        @test size(mat, 1) == 4
        @test size(mat, 2) == nrow(df)
        # Non-existent channels get dropped by intersect
        mat2 = EegFun._create_ica_data_matrix(df, [:ch1, :chX], 1:nrow(df))
        @test size(mat2, 1) == 1
    end

    @testset "remove and restore components" begin
        dat = create_synthetic_continuous()
        ica_res = EegFun.run_ica(dat; n_components = 3)
        # Non-mutating removal
        cleaned_df, ica_updated = EegFun.subtract_ica_components(dat.data, ica_res; component_selection = EegFun.components([1]))
        @test !isempty(ica_updated.removed_activations)
        @test cleaned_df isa DataFrame
        # Restore non-mutating
        restored_df, ica_restored = EegFun.add_ica_components(cleaned_df, ica_updated; component_selection = EegFun.components([1]))
        @test isempty(ica_restored.removed_activations)

        # Mutating on ContinuousData
        dat2 = copy(dat)
        ica2 = copy(ica_res)
        EegFun.subtract_ica_components!(dat2, ica2; component_selection = EegFun.components([1, 2]))
        @test length(keys(ica2.removed_activations)) == 2
        EegFun.add_ica_components!(dat2, ica2; component_selection = EegFun.components([1]))
        @test length(keys(ica2.removed_activations)) == 1

        # Invalid component index (no-op when selection yields no components)
        before = length(keys(ica2.removed_activations))
        EegFun.subtract_ica_components!(dat2, ica2; component_selection = EegFun.components([999]))
        @test length(keys(ica2.removed_activations)) == before
        EegFun.add_ica_components!(dat2.data, ica2; component_selection = EegFun.components([999]))
        @test length(keys(ica2.removed_activations)) == before

        # Restoring a valid component that was not removed should throw
        dat3 = create_synthetic_continuous()
        ica3 = EegFun.run_ica(dat3; n_components = 3)
        @test_throws EegFun.EegFunError EegFun.add_ica_components!(dat3.data, ica3; component_selection = EegFun.components([1]))

        # Roundtrip (mutating): remove then restore yields original data
        dat4 = create_synthetic_continuous()
        # Use full-rank ICA (components == channels)
        nfull = EegFun.n_channels(dat4)
        ica4 = EegFun.run_ica(dat4; n_components = nfull)
        orig_ch = select(copy(dat4.data, copycols = true), dat4.layout.data.label)
        EegFun.subtract_ica_components!(dat4, ica4; component_selection = EegFun.components([1, 2]))
        EegFun.add_ica_components!(dat4, ica4; component_selection = EegFun.components([1, 2]))
        roundtrip_ch = select(dat4.data, dat4.layout.data.label)
        @test all(isapprox.(Matrix(roundtrip_ch), Matrix(orig_ch); rtol = 1e-6, atol = 1e-8))

        # 1. Bare subtract without component_selection throws error
        @test_throws EegFun.EegFunError EegFun.subtract_ica_components!(dat4, ica4)

        # 2. Round-trip when n_components < n_channels (PCA reduction)
        dat_pca = create_synthetic_continuous(; nch = 6, n = 500)
        ica_pca = EegFun.run_ica(dat_pca; n_components = 3)
        orig_pca = Matrix{Float64}(dat_pca.data[!, dat_pca.layout.data.label])
        EegFun.subtract_ica_components!(dat_pca, ica_pca; component_selection = EegFun.components([1, 3]))
        @test dat_pca.analysis_info.n_ica_components_removed == 2
        # Verify activation of removed components is ~0
        X_sub = Matrix{Float64}(dat_pca.data[!, dat_pca.layout.data.label])
        acts_rem = (ica_pca.unmixing[[1, 3], :] * transpose(X_sub) .- ica_pca.unmixing[[1, 3], :] * ica_pca.mean) ./ ica_pca.scale
        @test maximum(abs.(acts_rem)) < 1e-10
        # Default add restores all removed components
        EegFun.add_ica_components!(dat_pca, ica_pca)
        @test dat_pca.analysis_info.n_ica_components_removed == 0
        restored_pca = Matrix{Float64}(dat_pca.data[!, dat_pca.layout.data.label])
        @test maximum(abs.(restored_pca - orig_pca)) < 1e-10

        # 3. Duplicate indices: components([1, 1]) removes once and increments counter by 1
        dat_dup = create_synthetic_continuous(; nch = 4, n = 200)
        ica_dup1 = EegFun.run_ica(dat_dup; n_components = 3)
        ica_dup2 = copy(ica_dup1)
        dat_dup1 = copy(dat_dup)
        dat_dup2 = copy(dat_dup)
        EegFun.subtract_ica_components!(dat_dup1, ica_dup1; component_selection = EegFun.components([1]))
        EegFun.subtract_ica_components!(dat_dup2, ica_dup2; component_selection = EegFun.components([1, 1]))
        @test dat_dup1.analysis_info.n_ica_components_removed == 1
        @test dat_dup2.analysis_info.n_ica_components_removed == 1
        @test Matrix(dat_dup1.data[!, dat_dup1.layout.data.label]) == Matrix(dat_dup2.data[!, dat_dup2.layout.data.label])

        # 4. Sequential subtract: component 2, then [2, 3], then restore [2]
        dat_seq = create_synthetic_continuous(; nch = 5, n = 300)
        ica_seq = EegFun.run_ica(dat_seq; n_components = 4)
        orig_seq = Matrix{Float64}(dat_seq.data[!, dat_seq.layout.data.label])
        EegFun.subtract_ica_components!(dat_seq, ica_seq; component_selection = EegFun.components([2]))
        @test dat_seq.analysis_info.n_ica_components_removed == 1
        # Subtracting [2, 3] skips 2 with warning and removes 3
        EegFun.subtract_ica_components!(dat_seq, ica_seq; component_selection = EegFun.components([2, 3]))
        @test dat_seq.analysis_info.n_ica_components_removed == 2
        @test sort(collect(keys(ica_seq.removed_activations))) == [2, 3]
        # Restoring 2 leaves 3 removed
        EegFun.add_ica_components!(dat_seq, ica_seq; component_selection = EegFun.components([2]))
        @test dat_seq.analysis_info.n_ica_components_removed == 1
        @test collect(keys(ica_seq.removed_activations)) == [3]
        # Now restoring 3 recovers original data completely
        EegFun.add_ica_components!(dat_seq, ica_seq; component_selection = EegFun.components([3]))
        @test dat_seq.analysis_info.n_ica_components_removed == 0
        @test maximum(abs.(Matrix{Float64}(dat_seq.data[!, dat_seq.layout.data.label]) - orig_seq)) < 1e-10

        # 5. Float32 DataFrame support
        df32 = DataFrame(:time => 1:100, :ch1 => Float32.(sin.(1:100)), :ch2 => Float32.(cos.(1:100)))
        lay32 = EegFun.Layout(DataFrame(:label => [:ch1, :ch2], :x => [0.0, 1.0], :y => [0.0, 0.0]), nothing, nothing, nothing)
        ica32 = EegFun.InfoIca(
            "test32",
            [1.0 0.0; 0.0 1.0],
            [1.0 0.0; 0.0 1.0],
            [1.0 0.0; 0.0 1.0],
            [1.0, 1.0],
            1.0,
            [0.0, 0.0],
            [:comp1, :comp2],
            Dict{Int,Matrix{Float64}}(),
            lay32,
            [false, false],
        )
        orig_ch1 = copy(df32.ch1)
        EegFun.subtract_ica_components!(df32, ica32; component_selection = EegFun.components([1]))
        @test eltype(df32.ch1) == Float32
        EegFun.add_ica_components!(df32, ica32)
        @test isapprox(df32.ch1, orig_ch1; atol = 1e-6)

        # 6. Epoched data (Vector{DataFrame} and EpochData)
        ep_dfs = [DataFrame(:ch1 => randn(50), :ch2 => randn(50)), DataFrame(:ch1 => randn(50), :ch2 => randn(50))]
        ep_dat = EegFun.EpochData("test_ep", 1, "c1", ep_dfs, lay32, 100, EegFun.AnalysisInfo())
        orig_ep_mat = Matrix(ep_dat.data[1][!, [:ch1, :ch2]])
        EegFun.subtract_ica_components!(ep_dat, ica32; component_selection = EegFun.components([1]))
        @test ep_dat.analysis_info.n_ica_components_removed == 1
        EegFun.add_ica_components!(ep_dat, ica32)
        @test ep_dat.analysis_info.n_ica_components_removed == 0
        @test isapprox(Matrix(ep_dat.data[1][!, [:ch1, :ch2]]), orig_ep_mat; atol = 1e-10)

        # Unequal epoch lengths throw EegFunError
        bad_dfs = [DataFrame(:ch1 => randn(50), :ch2 => randn(50)), DataFrame(:ch1 => randn(40), :ch2 => randn(40))]
        @test_throws EegFun.EegFunError EegFun.subtract_ica_components!(bad_dfs, ica32; component_selection = EegFun.components([1]))
    end

    @testset "artifact identification helpers" begin
        dat = create_synthetic_continuous()
        ica_res = EegFun.run_ica(dat)

        # EOG components (presence of vEOG and hEOG columns)
        eog_dict, eog_df = EegFun.identify_eog_components(dat, ica_res)
        @test eog_dict isa Dict
        @test all(haskey(eog_dict, k) for k in [:vEOG, :hEOG])
        @test eog_df isa DataFrame
        @test all(c in propertynames(eog_df) for c in [:Component, :vEOG_corr, :vEOG_zscore, :hEOG_corr, :hEOG_zscore])

        # Missing EOG columns -> returns empty dict with Int[] instead of crashing pipeline
        dat_no_eog = copy(dat)
        select!(dat_no_eog.data, Not([:vEOG, :hEOG]))
        result = EegFun.identify_eog_components(dat_no_eog, ica_res)
        @test result[1] == Dict(:vEOG => Int[], :hEOG => Int[])
        @test result[2] isa DataFrame

        # No samples selected -> returns empty dict and typed empty df
        empty_eog, empty_eog_df = EegFun.identify_eog_components(dat, ica_res; sample_selection = x -> falses(EegFun.n_samples(x)))
        @test empty_eog == Dict(:vEOG => Int[], :hEOG => Int[])
        @test nrow(empty_eog_df) == 0

        # ECG components
        ecg_vec, ecg_df = EegFun.identify_ecg_components(dat, ica_res)
        @test ecg_vec isa Vector{Int}
        @test ecg_df isa DataFrame
        @test all(
            c in propertynames(ecg_df) for
            c in [:Component, :num_peaks, :num_valid_ibis, :mean_ibi_s, :std_ibi_s, :peak_ratio, :heart_rate_bpm]
        )
        # No samples selected -> empty results with typed schema
        ecg_vec0, ecg_df0 = EegFun.identify_ecg_components(dat, ica_res; sample_selection = x -> falses(EegFun.n_samples(x)))
        @test isempty(ecg_vec0) && nrow(ecg_df0) == 0

        # Spatial kurtosis components
        sk_vec, sk_df = EegFun.identify_spatial_kurtosis_components(ica_res)
        @test sk_vec isa Vector{Int}
        @test sk_df isa DataFrame
        @test all(c in propertynames(sk_df) for c in [:Component, :spatial_kurtosis, :z_spatial_kurtosis])

        # Line noise components
        ln_vec, ln_df = EegFun.identify_line_noise_components(dat, ica_res)
        @test ln_vec isa Vector{Int}
        @test ln_df isa DataFrame
        @test all(
            c in propertynames(ln_df) for
            c in [:Component, :line_power, :flanking_power, :power_ratio, :harmonic_ratio, :power_ratio_zscore]
        )
        # No samples selected -> empty results with typed schema
        ln_vec0, ln_df0 = EegFun.identify_line_noise_components(dat, ica_res; sample_selection = x -> falses(EegFun.n_samples(x)))
        @test isempty(ln_vec0) && nrow(ln_df0) == 0

        # Unified identify_components with top-level kwargs forwarding
        artifacts_unified, metrics_unified = EegFun.identify_components(dat, ica_res; z_threshold = 2.5)
        @test artifacts_unified isa EegFun.ArtifactComponents
        @test all(haskey(metrics_unified, k) for k in [:eog_metrics, :ecg_metrics, :line_noise_metrics, :channel_noise_metrics])

        # Unified call on dataset without EOG channels does not crash
        artifacts_no_eog, _ = EegFun.identify_components(dat_no_eog, ica_res)
        @test artifacts_no_eog isa EegFun.ArtifactComponents
        @test isempty(artifacts_no_eog.eog[:vEOG])
        @test isempty(artifacts_no_eog.eog[:hEOG])
    end

    @testset "artifact combination helpers" begin
        artifacts = EegFun.combine_artifact_components(Dict(:vEOG => [1, 2], :hEOG => [3]), [4], [5, 6], [2, 7])
        @test artifacts isa EegFun.ArtifactComponents
        allc = EegFun.get_all_ica_components(artifacts)
        @test sort(allc) == sort([1, 2, 3, 4, 5, 6, 7])
        # get_all_components alias
        @test EegFun.get_all_components(artifacts) == allc

        # show single-line and MIME text/plain should not error
        io = IOBuffer()
        show(io, artifacts)
        @test startswith(String(take!(io)), "ArtifactComponents(")

        io_mime = IOBuffer()
        show(io_mime, MIME("text/plain"), artifacts)
        s_mime = String(take!(io_mime))
        @test occursin("ICA Artifact Components", s_mime)
        @test occursin("vEOG: [1, 2]", s_mime)
    end

    @testset "infomax and extended infomax correctness" begin
        # 1. Non-mutation of inputs
        rng = MersenneTwister(42)
        X = randn(rng, 4, 1000)
        X_orig = copy(X)
        prms = EegFun.IcaPrms(l_rate = 0.002)
        prms_orig = deepcopy(prms)
        layout = EegFun.Layout(DataFrame(label = [:C1, :C2, :C3, :C4], inc = zeros(4), azi = zeros(4)), nothing, nothing, nothing)

        res = EegFun.infomax_ica(X, layout, "test"; n_components = 3, params = prms, rng = MersenneTwister(1))
        @test X == X_orig
        @test prms.l_rate == prms_orig.l_rate

        # 2. Reproducibility test with rng
        res1 = EegFun.infomax_ica(X, layout, "test"; n_components = 3, rng = MersenneTwister(123))
        res2 = EegFun.infomax_ica(X, layout, "test"; n_components = 3, rng = MersenneTwister(123))
        @test isapprox(res1.unmixing, res2.unmixing; atol = 1e-12)

        # 3. Rank deficiency test (average referenced data has rank n_channels - 1)
        X_avg = X .- mean(X, dims = 1)
        @test_throws EegFun.EegFunError EegFun.infomax_ica(X_avg, layout, "test"; n_components = 4)

        # 4. Dimension guards
        @test_throws EegFun.EegFunError EegFun.infomax_ica(X, layout, "test"; n_components = 0)
        @test_throws EegFun.EegFunError EegFun.infomax_ica(X, layout, "test"; n_components = 5)

        # 5. Extended Infomax on sub-Gaussian sources (Lee et al. 1999 rule)
        N = 5000
        s1 = rand(MersenneTwister(1), N) .- 0.5
        s2 = rand(MersenneTwister(2), N) .- 0.5
        S = Matrix(transpose(hcat(s1, s2)))
        A = [0.8 0.6; 0.3 0.9]
        X_sub = A * S
        layout2 = EegFun.Layout(DataFrame(label = [:C1, :C2], inc = zeros(2), azi = zeros(2)), nothing, nothing, nothing)
        res_ext = EegFun.infomax_extended_ica(X_sub, layout2, "test"; n_components = 2, rng = MersenneTwister(42))
        S_rec = res_ext.unmixing * ((X_sub .- res_ext.mean) ./ res_ext.scale)
        c11 = abs(cor(S_rec[1, :], s1))
        c12 = abs(cor(S_rec[1, :], s2))
        c21 = abs(cor(S_rec[2, :], s1))
        c22 = abs(cor(S_rec[2, :], s2))
        @test max(c11, c12) > 0.95
        @test max(c21, c22) > 0.95
        @test count(res_ext.is_sub_gaussian) == 2

        # 6. Mixed super-Gaussian (Laplace) and sub-Gaussian (Uniform) sources
        u = rand(MersenneTwister(3), N) .- 0.5
        s_laplace = -sign.(u) .* log.(1.0 .- 2.0 .* abs.(u))
        s_uniform = rand(MersenneTwister(4), N) .- 0.5
        S_mix = Matrix(transpose(hcat(s_laplace, s_uniform)))
        X_mix = A * S_mix
        res_mix = EegFun.infomax_extended_ica(X_mix, layout2, "test"; n_components = 2, rng = MersenneTwister(99))
        S_rec_mix = res_mix.unmixing * ((X_mix .- res_mix.mean) ./ res_mix.scale)
        cm11 = abs(cor(S_rec_mix[1, :], s_laplace))
        cm12 = abs(cor(S_rec_mix[1, :], s_uniform))
        cm21 = abs(cor(S_rec_mix[2, :], s_laplace))
        cm22 = abs(cor(S_rec_mix[2, :], s_uniform))
        @test max(cm11, cm21) > 0.95
        @test max(cm12, cm22) > 0.95
        @test count(res_mix.is_sub_gaussian) == 1
    end

    @testset "picard and extended picard correctness" begin
        # 1. Non-mutation of inputs
        rng = MersenneTwister(42)
        X = randn(rng, 4, 1000)
        X_orig = copy(X)
        prms = EegFun.IcaPrms()
        prms_orig = deepcopy(prms)
        layout = EegFun.Layout(DataFrame(label = [:C1, :C2, :C3, :C4], inc = zeros(4), azi = zeros(4)), nothing, nothing, nothing)

        res_pic = EegFun.picard_ica(X, layout, "test"; n_components = 3, params = prms)
        @test X == X_orig
        @test prms.default_stop == prms_orig.default_stop

        # 2. Rank deficiency test (average referenced data has rank n_channels - 1)
        X_avg = X .- mean(X, dims = 1)
        @test_throws EegFun.EegFunError EegFun.picard_ica(X_avg, layout, "test"; n_components = 4)

        # 3. Dimension guards
        @test_throws EegFun.EegFunError EegFun.picard_ica(X, layout, "test"; n_components = 0)
        @test_throws EegFun.EegFunError EegFun.picard_ica(X, layout, "test"; n_components = 5)

        # 4. Extended Picard on mixed super-Gaussian (Laplace) and sub-Gaussian (Uniform) sources
        N = 5000
        u = rand(MersenneTwister(3), N) .- 0.5
        s_laplace = -sign.(u) .* log.(1.0 .- 2.0 .* abs.(u))
        s_uniform = rand(MersenneTwister(4), N) .- 0.5
        S_mix = Matrix(transpose(hcat(s_laplace, s_uniform)))
        A = [0.8 0.6; 0.3 0.9]
        X_mix = A * S_mix
        layout2 = EegFun.Layout(DataFrame(label = [:C1, :C2], inc = zeros(2), azi = zeros(2)), nothing, nothing, nothing)

        res_ext_pic = EegFun.picard_ica(X_mix, layout2, "test"; n_components = 2, extended = true)
        S_rec_mix = res_ext_pic.unmixing * ((X_mix .- res_ext_pic.mean) ./ res_ext_pic.scale)
        c11 = abs(cor(S_rec_mix[1, :], s_laplace))
        c12 = abs(cor(S_rec_mix[1, :], s_uniform))
        c21 = abs(cor(S_rec_mix[2, :], s_laplace))
        c22 = abs(cor(S_rec_mix[2, :], s_uniform))
        @test max(c11, c21) > 0.95
        @test max(c12, c22) > 0.95
        # Verify sub-Gaussian flags are correctly preserved
        @test count(res_ext_pic.is_sub_gaussian) == 1

        # 5. Standard Picard (extended=false) marks all sources as super-Gaussian
        res_std_pic = EegFun.picard_ica(X_mix, layout2, "test"; n_components = 2, extended = false)
        @test all(.!res_std_pic.is_sub_gaussian)
    end

    @testset "amica correctness" begin
        # 1. Non-mutation of inputs
        rng = MersenneTwister(42)
        X = randn(rng, 4, 1000)
        X_orig = copy(X)
        prms = EegFun.IcaPrms(amica_lrate = 0.05)
        prms_orig = deepcopy(prms)
        layout = EegFun.Layout(DataFrame(label = [:C1, :C2, :C3, :C4], inc = zeros(4), azi = zeros(4)), nothing, nothing, nothing)

        res_amica = EegFun.amica_ica(X, layout, "test"; n_components = 3, params = prms, rng = MersenneTwister(1))
        @test X == X_orig
        @test prms.amica_lrate == prms_orig.amica_lrate

        # 2. Reproducibility test with rng
        res1 = EegFun.amica_ica(X, layout, "test"; n_components = 3, rng = MersenneTwister(123))
        res2 = EegFun.amica_ica(X, layout, "test"; n_components = 3, rng = MersenneTwister(123))
        @test isapprox(res1.unmixing, res2.unmixing; atol = 1e-12)

        # 3. Rank deficiency test (average referenced data has rank n_channels - 1)
        X_avg = X .- mean(X, dims = 1)
        @test_throws EegFun.EegFunError EegFun.amica_ica(X_avg, layout, "test"; n_components = 4)

        # 4. Dimension guards
        @test_throws EegFun.EegFunError EegFun.amica_ica(X, layout, "test"; n_components = 0)
        @test_throws EegFun.EegFunError EegFun.amica_ica(X, layout, "test"; n_components = 5)

        # 5. Separation of mixed super-Gaussian (Laplace) and sub-Gaussian (Uniform) sources
        N = 5000
        u = rand(MersenneTwister(3), N) .- 0.5
        s_laplace = -sign.(u) .* log.(1.0 .- 2.0 .* abs.(u))
        s_uniform = rand(MersenneTwister(4), N) .- 0.5
        S_mix = Matrix(transpose(hcat(s_laplace, s_uniform)))
        A = [0.8 0.6; 0.3 0.9]
        X_mix = A * S_mix
        layout2 = EegFun.Layout(DataFrame(label = [:C1, :C2], inc = zeros(2), azi = zeros(2)), nothing, nothing, nothing)

        res_mix = EegFun.amica_ica(X_mix, layout2, "test"; n_components = 2, rng = MersenneTwister(42))
        S_rec_mix = res_mix.unmixing * ((X_mix .- res_mix.mean) ./ res_mix.scale)
        cm11 = abs(cor(S_rec_mix[1, :], s_laplace))
        cm12 = abs(cor(S_rec_mix[1, :], s_uniform))
        cm21 = abs(cor(S_rec_mix[2, :], s_laplace))
        cm22 = abs(cor(S_rec_mix[2, :], s_uniform))
        @test max(cm11, cm21) > 0.90
        @test max(cm12, cm22) > 0.90
        # Verify sub-Gaussian flag derived from kurtosis
        @test count(res_mix.is_sub_gaussian) == 1

        # 6. Rollback safeguard on extreme learning rate
        prms_wild = EegFun.IcaPrms(amica_lrate = 1e8, amica_max_iter = 5)
        res_wild = EegFun.amica_ica(X_mix, layout2, "test"; n_components = 2, params = prms_wild, rng = MersenneTwister(1))
        @test all(isfinite, res_wild.unmixing)
        @test all(isfinite, res_wild.mixing)
    end

    @testset "artifact identification advanced and robustness" begin
        dat = create_synthetic_continuous()
        ica_res = EegFun.run_ica(dat)

        # 1. two_step = false skips lag calculations and keeps spatial/temporal columns as NaN
        eog_1step, df_1step = EegFun.identify_eog_components(dat, ica_res; two_step = false)
        @test all(isnan, df_1step.vEOG_spatial_corr)
        @test all(isnan, df_1step.vEOG_temporal_corr)

        # 2. Robust z-scoring avoids outlier masking
        corrs = vcat([0.95, 0.92], 0.1 .+ 0.01 .* rand(MersenneTwister(42), 18))
        rob_z = EegFun._robust_z(corrs)
        @test rob_z[1] > 3.0 && rob_z[2] > 3.0

        # 3. ECG polarity invariance: unidirectional peak finder detects peaks in both signs
        ts = zeros(2000)
        for idx = 100:250:1900
            ts[idx] = 8.0
        end
        pk_pos = EegFun._find_unidirectional_peaks(ts; min_prominence_std = 3.0)
        pk_neg = EegFun._find_unidirectional_peaks(-ts; min_prominence_std = 3.0)
        @test length(pk_pos) == 8
        @test length(pk_neg) == 0

        # When inverted, negative peak finder detects them
        pk_inv_neg = EegFun._find_unidirectional_peaks(-(-ts); min_prominence_std = 3.0)
        @test length(pk_inv_neg) == 8

        # 4. Non-contiguous sample coordinates in ECG IBI calculation
        selected = vcat(1:950, 2050:3000)
        peaks_local = [500, 900, 1000]
        orig_samples = selected[peaks_local] # [500, 900, 2099]
        n_valid, _, _, _ = EegFun._calculate_ibi_metrics(orig_samples, 250.0, 0.5, 2.0)
        @test n_valid == 1 # only 500->900 is valid, seam gap is cleanly rejected

        # 5. Line noise with Nyquist guard and unobservable harmonic
        t_128 = (0:2559) ./ 128.0
        noise_50hz = sin.(2π * 50.0 .* t_128)
        df_128 = DataFrame(:C1 => noise_50hz, :C2 => randn(MersenneTwister(1), 2560))
        layout_128 = EegFun.Layout(DataFrame(label = [:C1, :C2], inc = zeros(2), azi = zeros(2)), nothing, nothing, nothing)
        cont_128 = EegFun.ContinuousData("test_128", df_128, layout_128, 128.0, EegFun.AnalysisInfo())
        ica_128 = EegFun.InfoIca(
            "test",
            [1.0 0.0; 0.0 1.0],
            [1.0 0.0; 0.0 1.0],
            [1.0 0.0; 0.0 1.0],
            [0.5, 0.5],
            1.0,
            zeros(2),
            [:IC1, :IC2],
            Dict{Int,Matrix{Float64}}(),
            layout_128,
            falses(2),
        )
        ln_comps, ln_df = EegFun.identify_line_noise_components(cont_128, ica_128; line_freq = 50.0, z_threshold = 0.5)
        @test 1 in ln_comps # Found even though 2nd harmonic (100 Hz) is above Nyquist (64 Hz)
    end

end
