using Test
using Random
using Statistics
using StatsBase
using LinearAlgebra
using EegFun

@testset "RSA" begin

    # ────────────────────────────────────────────────────────────
    # 1. Dissimilarity measures
    # ────────────────────────────────────────────────────────────
    @testset "dissimilarity measures" begin
        Random.seed!(500)
        a = [1.0, 2.0, 3.0, 4.0, 5.0]
        b = [5.0, 4.0, 3.0, 2.0, 1.0]
        c = copy(a)

        # Correlation: identical → 0, anti-correlated → 2
        @test EegFun._compute_dissimilarity(a, c, :correlation) ≈ 0.0 atol = 1e-10
        @test EegFun._compute_dissimilarity(a, b, :correlation) ≈ 2.0 atol = 1e-10

        # Spearman
        @test EegFun._compute_dissimilarity(a, c, :spearman) ≈ 0.0 atol = 1e-10
        @test EegFun._compute_dissimilarity(a, b, :spearman) ≈ 2.0 atol = 1e-10

        # Euclidean
        @test EegFun._compute_dissimilarity(a, c, :euclidean) ≈ 0.0 atol = 1e-10
        @test EegFun._compute_dissimilarity(a, b, :euclidean) > 0.0

        # Mahalanobis with identity → same as Euclidean
        inv_I = Matrix{Float64}(LinearAlgebra.I, 5, 5)
        mahal = EegFun._compute_dissimilarity(a, b, :mahalanobis; inv_covariance_matrix = inv_I)
        eucl = EegFun._compute_dissimilarity(a, b, :euclidean)
        @test mahal ≈ eucl atol = 1e-10

        # Mahalanobis without covariance matrix → error
        @test_throws EegFun.EegFunError EegFun._compute_dissimilarity(a, b, :mahalanobis)

        # Unknown measure → error
        @test_throws Exception EegFun._compute_dissimilarity(a, b, :unknown)
    end

    # ────────────────────────────────────────────────────────────
    # 2. RDM computation
    # ────────────────────────────────────────────────────────────
    @testset "RDM computation" begin
        Random.seed!(501)
        patterns = [randn(10) for _ in 1:4]

        rdm = EegFun._compute_rdm(patterns, :correlation)

        @test size(rdm) == (4, 4)
        # Diagonal should be zero
        for i in 1:4
            @test rdm[i, i] ≈ 0.0
        end
        # Symmetric
        for i in 1:4, j in 1:4
            @test rdm[i, j] ≈ rdm[j, i]
        end
        # All finite
        @test all(isfinite, rdm)
    end

    @testset "RDM in-place" begin
        Random.seed!(502)
        patterns = [randn(8) for _ in 1:3]
        rdm = zeros(Float64, 3, 3)

        result = EegFun._compute_rdm!(rdm, patterns, :euclidean)

        @test result === rdm  # same object
        @test all(diag(rdm) .≈ 0.0)
        @test issymmetric(rdm)
    end

    # ────────────────────────────────────────────────────────────
    # 3. Upper triangular extraction
    # ────────────────────────────────────────────────────────────
    @testset "upper triangular extraction" begin
        rdm = [0.0 1.0 2.0; 1.0 0.0 3.0; 2.0 3.0 0.0]
        vec = EegFun._extract_upper_triangular(rdm)

        @test length(vec) == 3  # n*(n-1)/2
        @test vec == [1.0, 2.0, 3.0]
    end

    # ────────────────────────────────────────────────────────────
    # 4. RDM normalization
    # ────────────────────────────────────────────────────────────
    @testset "normalize_rdm" begin
        rdm = [0.0 1.0 3.0; 1.0 0.0 2.0; 3.0 2.0 0.0]

        # :none returns copy
        rdm_none = EegFun.normalize_rdm(rdm; method = :none)
        @test rdm_none == rdm
        @test rdm_none !== rdm  # different object

        # :zscore
        rdm_z = EegFun.normalize_rdm(rdm; method = :zscore)
        upper = EegFun._extract_upper_triangular(rdm_z)
        @test mean(upper) ≈ 0.0 atol = 1e-10
        @test std(upper) ≈ 1.0 atol = 1e-10
        @test all(diag(rdm_z) .≈ 0.0)

        # :minmax
        rdm_mm = EegFun.normalize_rdm(rdm; method = :minmax)
        upper_mm = EegFun._extract_upper_triangular(rdm_mm)
        @test minimum(upper_mm) ≈ 0.0
        @test maximum(upper_mm) ≈ 1.0

        # :rank
        rdm_rank = EegFun.normalize_rdm(rdm; method = :rank)
        @test all(diag(rdm_rank) .≈ 0.0)

        # Unknown method → error
        @test_throws Exception EegFun.normalize_rdm(rdm; method = :invalid)
    end

    # ────────────────────────────────────────────────────────────
    # 5. Model RDM helpers
    # ────────────────────────────────────────────────────────────
    @testset "create_rdm_from_reaction_times" begin
        rts = [0.3, 0.5, 0.4]
        rdm = EegFun.create_rdm_from_reaction_times(rts)

        @test size(rdm) == (3, 3)
        @test all(diag(rdm) .≈ 0.0)
        @test rdm[1, 2] ≈ 0.2
        @test rdm[1, 3] ≈ 0.1
        @test rdm[2, 3] ≈ 0.1
        @test issymmetric(rdm)
    end

    @testset "create_rdm_from_categorical" begin
        cats = [1, 1, 2, 2, 3]
        rdm = EegFun.create_rdm_from_categorical(cats)

        @test size(rdm) == (5, 5)
        @test rdm[1, 2] ≈ 0.0  # same category
        @test rdm[1, 3] ≈ 1.0  # different category
        @test rdm[3, 4] ≈ 0.0  # same category
        @test issymmetric(rdm)
    end

    @testset "create_rdm_from_vectors" begin
        Random.seed!(503)
        vecs = [randn(20) for _ in 1:3]

        rdm = EegFun.create_rdm_from_vectors(vecs; dissimilarity_measure = :euclidean)

        @test size(rdm) == (3, 3)
        @test all(diag(rdm) .≈ 0.0)
        @test issymmetric(rdm)
        @test all(rdm .>= 0.0)
    end

    @testset "create_rdm_from_matrix" begin
        Random.seed!(504)
        mat = randn(3, 10)  # 3 conditions, 10 features

        rdm = EegFun.create_rdm_from_matrix(mat; dissimilarity_measure = :correlation)

        @test size(rdm) == (3, 3)
        @test all(diag(rdm) .≈ 0.0)
        @test issymmetric(rdm)
    end

    @testset "create_rdm_from_distances" begin
        dists = [0.2, 0.1, 0.1]  # (1,2), (1,3), (2,3)
        rdm = EegFun.create_rdm_from_distances(dists, 3)

        @test size(rdm) == (3, 3)
        @test rdm[1, 2] ≈ 0.2
        @test rdm[1, 3] ≈ 0.1
        @test rdm[2, 3] ≈ 0.1

        # Wrong length → error
        @test_throws Exception EegFun.create_rdm_from_distances([0.1, 0.2], 3)
    end

    @testset "create_rdm_from_similarity_ratings" begin
        sim = [1.0 2.0 5.0; 2.0 1.0 4.0; 5.0 4.0 1.0]

        rdm = EegFun.create_rdm_from_similarity_ratings(sim; convert_to_dissimilarity = true)

        @test size(rdm) == (3, 3)
        @test all(diag(rdm) .≈ 0.0)
        # Most similar pair (1,3) should have lowest dissimilarity
        @test rdm[1, 3] < rdm[1, 2]
    end

    # ────────────────────────────────────────────────────────────
    # 6. RSA on epoch data
    # ────────────────────────────────────────────────────────────
    @testset "rsa - basic" begin
        Random.seed!(510)
        cond1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 12, n_channels = 3, n = 50, fs = 1000)
        cond2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 12, n_channels = 3, n = 50, fs = 1000)
        cond3 = EegFun.create_test_epoch_data(condition = 3, n_epochs = 12, n_channels = 3, n = 50, fs = 1000)

        result = EegFun.rsa([cond1, cond2, cond3])

        @test result isa EegFun.RsaData
        @test length(result.condition_names) == 3
        @test length(result.times) == 50
        @test size(result.rdm) == (50, 3, 3)
        @test result.dissimilarity_measure == :correlation
        @test length(result.channels) == 3

        # RDMs should be symmetric with zero diagonal
        for t in 1:length(result.times)
            rdm_t = result.rdm[t, :, :]
            @test all(diag(rdm_t) .≈ 0.0)
            @test issymmetric(rdm_t)
        end
    end

    @testset "rsa - euclidean" begin
        Random.seed!(511)
        cond1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 10, n_channels = 3, n = 30, fs = 1000)
        cond2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 10, n_channels = 3, n = 30, fs = 1000)

        result = EegFun.rsa([cond1, cond2]; dissimilarity_measure = :euclidean)

        @test result.dissimilarity_measure == :euclidean
        @test all(result.rdm .>= 0.0)
    end

    @testset "rsa - no trial averaging" begin
        Random.seed!(512)
        cond1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 10, n_channels = 3, n = 30, fs = 1000)
        cond2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 10, n_channels = 3, n = 30, fs = 1000)

        result = EegFun.rsa([cond1, cond2]; average_trials = false)

        @test result isa EegFun.RsaData
        @test size(result.rdm, 1) == 30
    end

    @testset "rsa - channel selection" begin
        Random.seed!(513)
        cond1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 10, n_channels = 3, n = 30, fs = 1000)
        cond2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 10, n_channels = 3, n = 30, fs = 1000)

        result = EegFun.rsa([cond1, cond2]; channel_selection = EegFun.channels([:Ch1, :Ch2]))

        @test length(result.channels) == 2
        @test result.channels == [:Ch1, :Ch2]
    end

    @testset "rsa - error cases" begin
        Random.seed!(514)
        cond1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 10, n_channels = 3, n = 30, fs = 1000)

        # Need at least 2 conditions
        @test_throws Exception EegFun.rsa([cond1])

        # Empty input
        @test_throws Exception EegFun.rsa(EegFun.EpochData[])
    end

    @testset "rsa - batch (vector of participants)" begin
        Random.seed!(515)
        all_participants = Vector{Vector{EegFun.EpochData}}()
        for _ in 1:3
            c1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 10, n_channels = 3, n = 30, fs = 1000)
            c2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 10, n_channels = 3, n = 30, fs = 1000)
            push!(all_participants, [c1, c2])
        end

        results = EegFun.rsa(all_participants)

        @test length(results) == 3
        @test all(r -> r isa EegFun.RsaData, results)
    end

    # ────────────────────────────────────────────────────────────
    # 7. Model comparison
    # ────────────────────────────────────────────────────────────
    @testset "compare_models - static" begin
        Random.seed!(520)
        cond1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 10, n_channels = 3, n = 30, fs = 1000)
        cond2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 10, n_channels = 3, n = 30, fs = 1000)
        cond3 = EegFun.create_test_epoch_data(condition = 3, n_epochs = 10, n_channels = 3, n = 30, fs = 1000)

        rsa_result = EegFun.rsa([cond1, cond2, cond3])

        # Create a model RDM
        model_rdm = EegFun.create_rdm_from_reaction_times([0.3, 0.5, 0.4])

        result = EegFun.compare_models(rsa_result, [model_rdm]; model_names = ["RT"], n_permutations = 0)

        @test result isa EegFun.RsaData
        @test !isnothing(result.model_correlations)
        @test size(result.model_correlations) == (30, 1)
        @test result.model_names == ["RT"]
        @test all(isfinite, result.model_correlations)
    end

    @testset "compare_models - with permutations" begin
        Random.seed!(521)
        cond1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 10, n_channels = 3, n = 20, fs = 1000)
        cond2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 10, n_channels = 3, n = 20, fs = 1000)
        cond3 = EegFun.create_test_epoch_data(condition = 3, n_epochs = 10, n_channels = 3, n = 20, fs = 1000)

        rsa_result = EegFun.rsa([cond1, cond2, cond3])
        model_rdm = EegFun.create_rdm_from_reaction_times([0.3, 0.5, 0.4])

        result = EegFun.compare_models(rsa_result, [model_rdm]; model_names = ["RT"], n_permutations = 20)

        @test !isnothing(result.p_values)
        @test size(result.p_values) == (20, 1)
        @test all(0 .<= result.p_values .<= 1)
    end

    # ────────────────────────────────────────────────────────────
    # 8. Cross-validated RSA
    # ────────────────────────────────────────────────────────────
    @testset "rsa_crossvalidated - splithalf" begin
        Random.seed!(530)
        cond1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 20, n_channels = 3, n = 20, fs = 1000)
        cond2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 20, n_channels = 3, n = 20, fs = 1000)

        result = EegFun.rsa_crossvalidated([cond1, cond2]; cv_method = :splithalf, n_iterations = 5)

        @test result isa EegFun.RsaData
        @test size(result.rdm, 1) == 20
    end

    @testset "rsa_crossvalidated - leaveoneout" begin
        Random.seed!(531)
        cond1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 10, n_channels = 3, n = 20, fs = 1000)
        cond2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 10, n_channels = 3, n = 20, fs = 1000)

        result = EegFun.rsa_crossvalidated([cond1, cond2]; cv_method = :leaveoneout)

        @test result isa EegFun.RsaData
    end

    @testset "rsa_crossvalidated - kfold" begin
        Random.seed!(532)
        cond1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 15, n_channels = 3, n = 20, fs = 1000)
        cond2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 15, n_channels = 3, n = 20, fs = 1000)

        result = EegFun.rsa_crossvalidated([cond1, cond2]; cv_method = :kfold, n_folds = 3)

        @test result isa EegFun.RsaData
    end

    @testset "rsa_crossvalidated - error cases" begin
        Random.seed!(533)
        cond1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 3, n_channels = 3, n = 20, fs = 1000)
        cond2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 3, n_channels = 3, n = 20, fs = 1000)

        # kfold with too few trials → error
        @test_throws Exception EegFun.rsa_crossvalidated([cond1, cond2]; cv_method = :kfold, n_folds = 5)

        # Unknown method → error
        @test_throws Exception EegFun.rsa_crossvalidated([cond1, cond2]; cv_method = :invalid)
    end

    # ────────────────────────────────────────────────────────────
    # 9. Noise ceiling
    # ────────────────────────────────────────────────────────────
    @testset "compute_noise_ceiling" begin
        Random.seed!(540)
        rsa_list = EegFun.RsaData[]
        for p in 1:4
            c1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 10, n_channels = 3, n = 20, fs = 1000)
            c2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 10, n_channels = 3, n = 20, fs = 1000)
            c3 = EegFun.create_test_epoch_data(condition = 3, n_epochs = 10, n_channels = 3, n = 20, fs = 1000)
            push!(rsa_list, EegFun.rsa([c1, c2, c3]))
        end

        nc = EegFun.compute_noise_ceiling(rsa_list)

        @test nc isa EegFun.NoiseCeiling
        @test nc.n_participants == 4
        @test length(nc.lower_bound) == 20
        @test length(nc.upper_bound) == 20
        # Upper bound should generally be >= lower bound
        @test all(isfinite, nc.lower_bound)
        @test all(isfinite, nc.upper_bound)
    end

    @testset "compute_noise_ceiling - too few participants" begin
        Random.seed!(541)
        c1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 10, n_channels = 3, n = 20, fs = 1000)
        c2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 10, n_channels = 3, n = 20, fs = 1000)
        single = [EegFun.rsa([c1, c2])]

        @test_throws Exception EegFun.compute_noise_ceiling(single)
    end

    # ────────────────────────────────────────────────────────────
    # 10. Grand average
    # ────────────────────────────────────────────────────────────
    @testset "grand_average (RsaData)" begin
        Random.seed!(550)
        rsa_list = EegFun.RsaData[]
        for p in 1:3
            c1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 10, n_channels = 3, n = 20, fs = 1000)
            c2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 10, n_channels = 3, n = 20, fs = 1000)
            push!(rsa_list, EegFun.rsa([c1, c2]))
        end

        ga = EegFun.grand_average(rsa_list)

        @test ga isa EegFun.RsaData
        @test ga.file == "Grand Average"
        @test size(ga.rdm) == size(rsa_list[1].rdm)
        @test !isnothing(ga.noise_ceiling)  # auto-computed with >= 2 participants

        # Grand average RDM should be mean of individual RDMs
        expected = mean([r.rdm for r in rsa_list])
        @test ga.rdm ≈ expected

        # Single participant returns itself
        single_ga = EegFun.grand_average([rsa_list[1]])
        @test single_ga === rsa_list[1]
    end

    # ────────────────────────────────────────────────────────────
    # 11. Temporal resampling
    # ────────────────────────────────────────────────────────────
    @testset "resample_temporal_data" begin
        Random.seed!(560)
        # 3 conditions, 2 features, 10 time points
        data = randn(3, 2, 10)
        original_times = collect(0.0:0.1:0.9)
        target_times = collect(0.0:0.05:0.9)  # double the rate

        resampled, _ = EegFun.resample_temporal_data(data, original_times, target_times; method = :linear)

        @test size(resampled) == (3, 2, length(target_times))
        @test all(isfinite, resampled)

        # At original time points, values should match (or be very close)
        @test resampled[1, 1, 1] ≈ data[1, 1, 1] atol = 1e-10
    end

    @testset "resample_temporal_data - nearest" begin
        Random.seed!(561)
        data = randn(2, 1, 5)
        original_times = collect(0.0:0.1:0.4)
        target_times = collect(0.0:0.05:0.4)

        resampled, _ = EegFun.resample_temporal_data(data, original_times, target_times; method = :nearest)

        @test size(resampled) == (2, 1, length(target_times))
    end

    # ────────────────────────────────────────────────────────────
    # 12. Timeseries RDM creation
    # ────────────────────────────────────────────────────────────
    @testset "create_rdm_from_timeseries" begin
        Random.seed!(570)
        # 3 conditions, 4 features, 10 timepoints
        data = randn(3, 4, 10)
        times = collect(0.0:0.1:0.9)

        rdms = EegFun.create_rdm_from_timeseries(data, times; dissimilarity_measure = :euclidean)

        @test size(rdms) == (10, 3, 3)
        # Diagonal should be zero at each time
        for t in 1:10
            @test all(diag(rdms[t, :, :]) .≈ 0.0)
        end
    end

    # ────────────────────────────────────────────────────────────
    # 13. create_model_rdms dispatch
    # ────────────────────────────────────────────────────────────
    @testset "create_model_rdms" begin
        model_data = Dict{String,Any}(
            "RTs" => [0.3, 0.5, 0.4],
            "Categories" => [1, 1, 2],
        )

        rdms, names = EegFun.create_model_rdms(model_data)

        @test length(rdms) == 2
        @test length(names) == 2
        @test all(rdm -> size(rdm) == (3, 3), rdms)
    end

    # ────────────────────────────────────────────────────────────
    # 14. Show methods
    # ────────────────────────────────────────────────────────────
    @testset "RsaData show" begin
        Random.seed!(580)
        c1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 10, n_channels = 3, n = 20, fs = 1000)
        c2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 10, n_channels = 3, n = 20, fs = 1000)

        result = EegFun.rsa([c1, c2])

        io = IOBuffer()
        show(io, result)
        output = String(take!(io))

        @test contains(output, "RsaData")
        @test contains(output, "correlation")
    end

    @testset "NoiseCeiling show" begin
        nc = EegFun.NoiseCeiling([0.5, 0.6], [0.8, 0.9], 4)

        io = IOBuffer()
        show(io, nc)
        output = String(take!(io))

        @test contains(output, "NoiseCeiling")
        @test contains(output, "4")
    end

end # @testset "RSA"
