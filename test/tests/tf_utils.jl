using Test
using DataFrames

@testset "TF Utils" begin

    @testset "_tf_build_row_index" begin
        df = DataFrame(freq = [1.0, 2.0, 1.0, 2.0], time = [0.0, 0.0, 0.1, 0.1], Ch1 = [10.0, 20.0, 30.0, 40.0])
        idx = EegFun._tf_build_row_index(df)
        @test length(idx) == 4
        @test idx[(1.0, 0.0)] == 1
        @test idx[(2.0, 0.0)] == 2
        @test idx[(1.0, 0.1)] == 3
        @test idx[(2.0, 0.1)] == 4
    end

    @testset "_tf_df_to_matrix" begin
        # Create a simple TF DataFrame: 3 freqs × 4 times
        freqs = [5.0, 10.0, 15.0]
        times = [0.0, 0.1, 0.2, 0.3]
        n_f = length(freqs)
        n_t = length(times)

        freq_col = repeat(freqs, inner = n_t)
        time_col = repeat(times, outer = n_f)
        ch1_vals = collect(1.0:(n_f*n_t))

        df = DataFrame(freq = freq_col, time = time_col, Ch1 = ch1_vals)

        mat = EegFun._tf_df_to_matrix(df, :Ch1, freqs, times)

        @test size(mat) == (3, 4)

        # Verify all values roundtrip correctly
        row_index = EegFun._tf_build_row_index(df)
        for (fi, f) in enumerate(freqs)
            for (ti, t) in enumerate(times)
                row = row_index[(round(f, digits = 6), round(t, digits = 6))]
                @test mat[fi, ti] == df[row, :Ch1]
            end
        end

        # Test with subset of frequencies
        mat_sub = EegFun._tf_df_to_matrix(df, :Ch1, [5.0, 15.0], times)
        @test size(mat_sub) == (2, 4)
        @test mat_sub[1, :] == mat[1, :]  # 5 Hz row
        @test mat_sub[2, :] == mat[3, :]  # 15 Hz row

        # Test with missing freq/time → NaN
        mat_miss = EegFun._tf_df_to_matrix(df, :Ch1, [5.0, 99.0], times)
        @test size(mat_miss) == (2, 4)
        @test !any(isnan, mat_miss[1, :])  # 5 Hz exists
        @test all(isnan, mat_miss[2, :])   # 99 Hz doesn't exist
    end

    @testset "_tf_matrix_to_df!" begin
        freqs = [5.0, 10.0]
        times = [0.0, 0.1, 0.2]
        n_f = length(freqs)
        n_t = length(times)

        freq_col = repeat(freqs, inner = n_t)
        time_col = repeat(times, outer = n_f)

        df = DataFrame(freq = freq_col, time = time_col, Ch1 = zeros(n_f * n_t))

        mat = [1.0 2.0 3.0; 4.0 5.0 6.0]  # [2 freqs × 3 times]
        EegFun._tf_matrix_to_df!(df, :Ch1, mat, freqs, times)

        # Read back and verify
        mat_back = EegFun._tf_df_to_matrix(df, :Ch1, freqs, times)
        @test mat_back == mat
    end

    @testset "Roundtrip with float tolerance" begin
        # Test that slight float variations still match
        freqs = [1.0 / 3.0, 2.0 / 3.0]  # imprecise floats
        times = [0.0, 0.1]

        freq_col = repeat(freqs, inner = 2)
        time_col = repeat(times, outer = 2)
        df = DataFrame(freq = freq_col, time = time_col, Ch1 = [10.0, 20.0, 30.0, 40.0])

        mat = EegFun._tf_df_to_matrix(df, :Ch1, freqs, times)
        @test !any(isnan, mat)  # All values should be found
        @test mat[1, 1] == 10.0
        @test mat[1, 2] == 20.0
        @test mat[2, 1] == 30.0
        @test mat[2, 2] == 40.0
    end

    @testset "Row index reuse across channels" begin
        freqs = [1.0, 2.0]
        times = [0.0, 0.1]
        n = 4

        df = DataFrame(
            freq = repeat(freqs, inner = 2),
            time = repeat(times, outer = 2),
            Ch1 = [1.0, 2.0, 3.0, 4.0],
            Ch2 = [10.0, 20.0, 30.0, 40.0],
        )

        row_index = EegFun._tf_build_row_index(df)

        mat1 = EegFun._tf_df_to_matrix(df, :Ch1, freqs, times, row_index)
        mat2 = EegFun._tf_df_to_matrix(df, :Ch2, freqs, times, row_index)

        @test mat1[1, 1] == 1.0
        @test mat2[1, 1] == 10.0
        @test mat1[2, 2] == 4.0
        @test mat2[2, 2] == 40.0
    end
end
