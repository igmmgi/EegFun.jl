using Test
using DataFrames
using Random
using Statistics
using LinearAlgebra
using eegfun

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
        df = DataFrame(time = t, sample = 1:n, triggers = zeros(Int, n))
        for (k, v) in cols
            df[!, k] = v
        end
        df[!, :vEOG] = vEOG
        df[!, :hEOG] = hEOG
        # Keep mask: first half only
        df[!, :keepmask] = [i <= n ÷ 2 for i = 1:n]
        layout = eegfun.Layout(
            DataFrame(label = Symbol.(collect(keys(cols))), inc = zeros(nch), azi = zeros(nch)),
            nothing,
            nothing,
        )
        dat = eegfun.ContinuousData("test_data", copy(df, copycols = true), layout, fs, eegfun.AnalysisInfo())
        return dat
    end

    @testset "run_ica basic and parameters" begin
        dat = create_synthetic_continuous()
        ica_res = eegfun.run_ica(dat)
        @test ica_res isa eegfun.InfoIca
        @test size(ica_res.unmixing, 1) == eegfun.n_channels(dat) - 1
        @test length(ica_res.ica_label) == size(ica_res.unmixing, 1)
        @test isapprox(sum(ica_res.variance), 1.0; atol = 1e-6)

        # Custom number of components
        ica_res2 = eegfun.run_ica(dat; n_components = 3)
        @test size(ica_res2.unmixing, 1) == 3
        # n_components larger than channels -> adjusted to channels-1
        ica_res2b = eegfun.run_ica(dat; n_components = 999)
        @test size(ica_res2b.unmixing, 1) == eegfun.n_channels(dat) - 1

        # Channel selection subset
        ica_res3 = eegfun.run_ica(dat; channel_selection = eegfun.channels([:ch1, :ch2, :ch3]))
        @test size(ica_res3.unmixing, 1) == 2

        # Sample selection subset
        ica_res4 = eegfun.run_ica(dat; sample_selection = eegfun.samples(:keepmask))
        @test ica_res4 isa eegfun.InfoIca

        # Preprocessing should be external; run_ica has no filter flags

        # Include extra channels (vEOG, hEOG) in ICA
        nch = eegfun.n_channels(dat)
        sel = vcat(dat.layout.data.label, [:vEOG, :hEOG])
        ica_extra = eegfun.run_ica(dat; include_extra = true, channel_selection = eegfun.channels(sel))
        @test size(ica_extra.unmixing, 1) == (nch + 2 - 1)

        # Error when no channels selected
        @test_throws Any eegfun.run_ica(dat; channel_selection = eegfun.channels(Symbol[]))
        # Error when no samples selected
        @test_throws Any eegfun.run_ica(dat; sample_selection = x -> falses(nrow(x)))

    end

    @testset "create_ica_data_matrix" begin
        dat = create_synthetic_continuous(; n = 200, fs = 200, nch = 4)
        # Duplicate samples to test unique by :sample
        df = vcat(dat.data, dat.data)
        mat = eegfun.create_ica_data_matrix(df, [:ch1, :ch2, :ch3, :ch4], 1:nrow(df))
        @test size(mat, 1) == 4
        @test size(mat, 2) == nrow(df)
        # Non-existent channels get dropped by intersect
        mat2 = eegfun.create_ica_data_matrix(df, [:ch1, :chX], 1:nrow(df))
        @test size(mat2, 1) == 1
    end

    @testset "remove and restore components" begin
        dat = create_synthetic_continuous()
        ica_res = eegfun.run_ica(dat; n_components = 3)
        # Non-mutating removal
        cleaned_df, ica_updated =
            eegfun.remove_ica_components(dat.data, ica_res; component_selection = eegfun.components([1]))
        @test !isempty(ica_updated.removed_activations)
        @test cleaned_df isa DataFrame
        # Restore non-mutating
        restored_df, ica_restored =
            eegfun.restore_ica_components(cleaned_df, ica_updated; component_selection = eegfun.components([1]))
        @test isempty(ica_restored.removed_activations)

        # Mutating on ContinuousData
        dat2 = copy(dat);
        ica2 = copy(ica_res)
        eegfun.remove_ica_components!(dat2, ica2; component_selection = eegfun.components([1, 2]))
        @test length(keys(ica2.removed_activations)) == 2
        eegfun.restore_ica_components!(dat2, ica2; component_selection = eegfun.components([1]))
        @test length(keys(ica2.removed_activations)) == 1

        # Invalid component index (no-op when selection yields no components)
        before = length(keys(ica2.removed_activations))
        eegfun.remove_ica_components!(dat2, ica2; component_selection = eegfun.components([999]))
        @test length(keys(ica2.removed_activations)) == before
        eegfun.restore_ica_components!(dat2.data, ica2; component_selection = eegfun.components([999]))
        @test length(keys(ica2.removed_activations)) == before

        # Restoring a valid component that was not removed should throw
        dat3 = create_synthetic_continuous();
        ica3 = eegfun.run_ica(dat3; n_components = 3)
        @test_throws ArgumentError eegfun.restore_ica_components!(
            dat3.data,
            ica3;
            component_selection = eegfun.components([1]),
        )

        # Roundtrip (mutating): remove then restore yields original data
        dat4 = create_synthetic_continuous()
        # Use full-rank ICA (components == channels) to make the transform invertible
        nfull = eegfun.n_channels(dat4)
        ica4 = eegfun.run_ica(dat4; n_components = nfull)
        orig_ch = select(copy(dat4.data, copycols = true), dat4.layout.data.label)
        eegfun.remove_ica_components!(dat4, ica4; component_selection = eegfun.components([1, 2]))
        eegfun.restore_ica_components!(dat4, ica4; component_selection = eegfun.components([1, 2]))
        roundtrip_ch = select(dat4.data, dat4.layout.data.label)
        @test all(isapprox.(Matrix(roundtrip_ch), Matrix(orig_ch); rtol = 1e-6, atol = 1e-8))

        # Note: Non-mutating roundtrip is sensitive to numeric conditioning; covered by mutating check above.
    end

    @testset "artifact identification helpers" begin
        dat = create_synthetic_continuous()
        ica_res = eegfun.run_ica(dat)

        # EOG components (presence of vEOG and hEOG columns)
        eog_dict, eog_df = eegfun.identify_eog_components(dat, ica_res)
        @test eog_dict isa Dict
        @test all(haskey(eog_dict, k) for k in [:vEOG, :hEOG])
        @test eog_df isa DataFrame
        @test all(c in propertynames(eog_df) for c in [:Component, :vEOG_corr, :vEOG_zscore, :hEOG_corr, :hEOG_zscore])
        # Missing EOG columns -> returns nothing
        dat_no_eog = copy(dat)
        select!(dat_no_eog.data, Not([:vEOG, :hEOG]))
        @test eegfun.identify_eog_components(dat_no_eog, ica_res) === nothing
        # No samples selected -> returns empty dict and empty df
        empty_eog, empty_eog_df = eegfun.identify_eog_components(dat, ica_res; sample_selection = x -> falses(nrow(x)))
        @test empty_eog == Dict(:vEOG => Int[], :hEOG => Int[])
        @test size(empty_eog_df) == (0, 0)

        # ECG components
        ecg_vec, ecg_df = eegfun.identify_ecg_components(dat, ica_res)
        @test ecg_vec isa Vector{Int}
        @test ecg_df isa DataFrame
        @test all(
            c in propertynames(ecg_df) for c in [
                :Component,
                :num_peaks,
                :num_valid_ibis,
                :mean_ibi_s,
                :std_ibi_s,
                :peak_ratio,
                :heart_rate_bpm,
                :is_ecg_artifact,
            ]
        )
        # No samples selected -> empty results
        ecg_vec0, ecg_df0 = eegfun.identify_ecg_components(dat, ica_res; sample_selection = x -> falses(nrow(x)))
        @test isempty(ecg_vec0) && size(ecg_df0) == (0, 0)

        # Spatial kurtosis components
        sk_vec, sk_df = eegfun.identify_spatial_kurtosis_components(ica_res)
        @test sk_vec isa Vector{Int}
        @test sk_df isa DataFrame
        @test all(c in propertynames(sk_df) for c in [:Component, :SpatialKurtosis, :SpatialKurtosisZScore])

        # Line noise components
        ln_vec, ln_df = eegfun.identify_line_noise_components(dat, ica_res)
        @test ln_vec isa Vector{Int}
        @test ln_df isa DataFrame
        @test all(
            c in propertynames(ln_df) for c in [
                :Component,
                :LinePower,
                :SurroundingPower,
                :PowerRatio,
                :Harmonic2Ratio,
                :Harmonic3Ratio,
                :PowerRatioZScore,
            ]
        )
        # No samples selected -> empty results
        ln_vec0, ln_df0 = eegfun.identify_line_noise_components(dat, ica_res; sample_selection = x -> falses(nrow(x)))
        @test isempty(ln_vec0) && size(ln_df0) == (0, 0)
    end

    @testset "artifact combination helpers" begin
        artifacts = eegfun.combine_artifact_components(Dict(:vEOG => [1, 2], :hEOG => [3]), [4], [5, 6], [2, 7])
        @test artifacts isa eegfun.ArtifactComponents
        allc = eegfun.get_all_ica_components(artifacts)
        @test sort(allc) == sort([1, 2, 3, 4, 5, 6, 7])
        # show should not error
        io = IOBuffer()
        show(io, artifacts)
        @test length(String(take!(io))) > 0
    end

end
