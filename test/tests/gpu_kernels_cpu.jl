using Test
using EegFun
using DataFrames
using Random

@testset "GPU Kernel Execution via KernelAbstractions CPU Backend" begin
    # 1. Verify graceful CPU fallback when use_gpu = true and no GPU backend is active
    c_dat = EegFun.create_test_continuous_data(n = 200, fs = 200, n_channels = 3, seed = 42)
    ep_dat = EegFun.create_test_epoch_data(n = 64, fs = 128, n_epochs = 4, n_channels = 2, seed = 42)

    if !EegFun.is_gpu_available()
        ica_fb =
            EegFun.run_ica(c_dat; algorithm = :infomax, use_gpu = true, params = EegFun.IcaPrms(max_iter = 5), rng = MersenneTwister(1))
        @test ica_fb isa EegFun.InfoIca

        tf_fb = EegFun.tf_morlet(ep_dat; frequencies = [10.0, 20.0], cycles = 3, use_gpu = true)
        @test tf_fb isa EegFun.TimeFreqData
    end

    # 2. Temporarily register KernelAbstractions.CPU() to exercise all GPU kernels & CPU-vs-GPU equivalence
    prev_state = (
        EegFun.GLOBAL_GPU_STATE.active,
        EegFun.GLOBAL_GPU_STATE.backend_name,
        EegFun.GLOBAL_GPU_STATE.device_name,
        EegFun.GLOBAL_GPU_STATE.backend_obj,
        EegFun.GLOBAL_GPU_STATE.converter,
    )

    try
        EegFun.register_gpu_backend!(:CPU_KA, "KernelAbstractions CPU", EegFun.KernelAbstractions.CPU(), identity)
        @test EegFun.is_gpu_available() == true
        @test EegFun.gpu_device_name() == "KernelAbstractions CPU"
        @test EegFun.gpu_backend() == EegFun.KernelAbstractions.CPU()
        arr = Float32[1.0, 2.0, 3.0]
        @test EegFun.gpu_array(arr) === arr

        # 2a. Time-Frequency CPU vs GPU numerical equivalence (Morlet, Multitaper, STFT)
        @testset "TF GPU kernels (Morlet, Multitaper, STFT)" begin
            for (fn, kw) in [
                (EegFun.tf_morlet, (; frequencies = [8.0, 16.0], cycles = 3)),
                (EegFun.tf_multitaper, (; frequencies = [8.0, 16.0], cycles = 4)),
                (EegFun.tf_stft, (; frequencies = [8.0, 16.0], window_length = 0.25)),
            ]
                tf_cpu = fn(ep_dat; use_gpu = false, pad = :both, kw...)
                tf_gpu = fn(ep_dat; use_gpu = true, pad = :both, kw...)
                @test isapprox(tf_cpu.data_power.Ch1, tf_gpu.data_power.Ch1; rtol = 1e-5, nans = true)
                @test isapprox(tf_cpu.data_phase.Ch1, tf_gpu.data_phase.Ch1; rtol = 1e-5, nans = true)

                # Single-trial output mode (return_trials = true)
                tf_trials_cpu = fn(ep_dat; use_gpu = false, return_trials = true, kw...)
                tf_trials_gpu = fn(ep_dat; use_gpu = true, return_trials = true, kw...)
                @test isapprox(tf_trials_cpu.data_power[1].Ch1, tf_trials_gpu.data_power[1].Ch1; rtol = 1e-5, nans = true)
            end
        end

        # 2b. ICA GPU kernels (Infomax, Extended Infomax, Picard, Extended Picard, AMICA)
        @testset "ICA GPU kernels" begin
            prms = EegFun.IcaPrms(max_iter = 12, amica_max_iter = 15, default_stop = 1e-12)
            for alg in (:infomax, :infomax_extended, :picard, :picard_extended, :amica)
                res = EegFun.run_ica(c_dat; algorithm = alg, use_gpu = true, params = prms, rng = MersenneTwister(7))
                @test res isa EegFun.InfoIca
                @test size(res.unmixing, 1) == 3
                @test all(isfinite, res.unmixing)
            end
        end

        # 2c. ERP Permutation Statistics GPU kernels (paired & independent)
        @testset "Permutation Statistics GPU kernels" begin
            layout = EegFun.create_test_layout(n_channels = 3, layout_type = :topo)
            chs = layout.data.label[1:3]
            layout = EegFun.subset_layout(layout; channel_selection = EegFun.channels(chs))
            EegFun.get_neighbors_xy!(layout, 1.5)
            t_vec = collect(-0.1:0.02:0.3)
            n_t = length(t_vec)

            # Paired design
            paired_erps = EegFun.ErpData[]
            for p = 1:5, (cond, amp) in ((1, 3.0), (2, 0.0))
                df = DataFrame(
                    time = t_vec,
                    sample = 1:n_t,
                    condition = fill(cond, n_t),
                    condition_name = fill("c$cond", n_t),
                    participant = fill(p, n_t),
                )
                rng_p = MersenneTwister(100 * p + cond)
                for ch in chs
                    df[!, ch] = amp .* exp.(-((t_vec .- 0.1) ./ 0.05) .^ 2) .+ 0.2 .* randn(rng_p, n_t)
                end
                push!(paired_erps, EegFun.ErpData("sub$p", cond, "c$cond", df, copy(layout), 50, EegFun.AnalysisInfo(), 20))
            end

            prep_paired = EegFun.prepare_stats(paired_erps; design = :paired)
            res_cpu = EegFun.permutation_test(
                prep_paired;
                n_permutations = 12,
                cluster_type = :spatiotemporal,
                use_gpu = false,
                show_progress = false,
                rng = MersenneTwister(42),
            )
            res_gpu = EegFun.permutation_test(
                prep_paired;
                n_permutations = 12,
                cluster_type = :spatiotemporal,
                use_gpu = true,
                show_progress = false,
                rng = MersenneTwister(42),
            )
            @test isapprox(res_cpu.stat_matrix.t, res_gpu.stat_matrix.t; rtol = 1e-4)
            @test isapprox(res_cpu.permutation_distribution.positive, res_gpu.permutation_distribution.positive; rtol = 1e-3)

            # Independent design + nonparametric_common threshold on GPU t-matrices
            indep_erps = EegFun.ErpData[]
            for p = 1:8
                cond = p <= 4 ? 1 : 2
                amp = cond == 1 ? 3.0 : 0.0
                df = DataFrame(
                    time = t_vec,
                    sample = 1:n_t,
                    condition = fill(cond, n_t),
                    condition_name = fill("c$cond", n_t),
                    participant = fill(p, n_t),
                )
                rng_p = MersenneTwister(200 + p)
                for ch in chs
                    df[!, ch] = amp .* sin.(2π * 5 .* t_vec) .+ 0.2 .* randn(rng_p, n_t)
                end
                push!(indep_erps, EegFun.ErpData("sub$p", cond, "c$cond", df, copy(layout), 50, EegFun.AnalysisInfo(), 20))
            end

            prep_indep = EegFun.prepare_stats(indep_erps; design = :independent)
            res_indep_gpu = EegFun.permutation_test(
                prep_indep;
                n_permutations = 12,
                threshold_method = :nonparametric_common,
                cluster_type = :spatiotemporal,
                use_gpu = true,
                show_progress = false,
                rng = MersenneTwister(99),
            )
            @test res_indep_gpu isa EegFun.PermutationResult
            @test length(res_indep_gpu.permutation_distribution.positive) == 12
        end
    finally
        EegFun.GLOBAL_GPU_STATE.active = prev_state[1]
        EegFun.GLOBAL_GPU_STATE.backend_name = prev_state[2]
        EegFun.GLOBAL_GPU_STATE.device_name = prev_state[3]
        EegFun.GLOBAL_GPU_STATE.backend_obj = prev_state[4]
        EegFun.GLOBAL_GPU_STATE.converter = prev_state[5]
    end
end
