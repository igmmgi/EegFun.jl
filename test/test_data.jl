using Test
using DataFrames
using eegfun

@testset "Data Utilities" begin

    # === TEST DATA SETUP ===
    function create_test_continuous_data()
        # Create test data with metadata, channels, and extra columns
        df = DataFrame(
            time = (0:9) ./ 1000,          # metadata
            triggers = [0, 1, 0, 0, 1, 0, 0, 0, 1, 0],  # metadata
            Fp1 = 1:10,                    # channel
            Fp2 = 11:20,                   # channel
            F3 = 21:30,                    # channel
            F4 = 31:40,                    # channel
            vEOG = 101:110,                # extra
            hEOG = 111:120                 # extra
        )
        layout = eegfun.Layout(
            DataFrame(label = [:Fp1, :Fp2, :F3, :F4], 
                     inc = [30, 30, 45, 45], 
                     azi = [0, 90, 180, 270]), 
            nothing, nothing
        )
        analysis_info = eegfun.AnalysisInfo(reference = :avg, hp_filter = 0.1, lp_filter = 30.0)
        return eegfun.ContinuousData(df, layout, 1000, analysis_info)
    end

    function create_test_epoch_data()
        # Create test epoch data
        epoch1 = DataFrame(
            time = [0.1, 0.2, 0.3],
            triggers = [1, 0, 0],
            Fp1 = [1, 2, 3],
            Fp2 = [4, 5, 6],
            vEOG = [7, 8, 9]
        )
        epoch2 = DataFrame(
            time = [0.1, 0.2, 0.3], 
            triggers = [1, 0, 0],
            Fp1 = [10, 20, 30],
            Fp2 = [40, 50, 60],
            vEOG = [70, 80, 90]
        )
        layout = eegfun.Layout(
            DataFrame(label = [:Fp1, :Fp2], 
                     inc = [30, 30], 
                     azi = [0, 90]), 
            nothing, nothing
        )
        analysis_info = eegfun.AnalysisInfo()
        return eegfun.EpochData([epoch1, epoch2], layout, 1000, analysis_info)
    end

    function create_test_erp_data()
        # Create test ERP data
        df = DataFrame(
            time = [-0.1, 0.0, 0.1, 0.2],
            Fp1 = [1.0, 2.0, 3.0, 4.0],
            Fp2 = [5.0, 6.0, 7.0, 8.0]
        )
        layout = eegfun.Layout(
            DataFrame(label = [:Fp1, :Fp2], 
                     inc = [30, 30], 
                     azi = [0, 90]), 
            nothing, nothing
        )
        analysis_info = eegfun.AnalysisInfo()
        return eegfun.ErpData(df, layout, 1000, analysis_info, 50)  # 50 epochs averaged
    end

    # === COLUMN IDENTIFICATION TESTS ===
    @testset "get_cols_by_group" begin
        dat = create_test_continuous_data()

        # Test channels group
        channels = eegfun.get_cols_by_group(dat, :channels)
        @test channels == [:Fp1, :Fp2, :F3, :F4]

        # Test metadata group  
        metadata = eegfun.get_cols_by_group(dat, :metadata)
        @test metadata == [:time, :triggers]

        # Test extra group
        extra = eegfun.get_cols_by_group(dat, :extra)
        @test extra == [:vEOG, :hEOG]

        # Test with EpochData
        epoch_dat = create_test_epoch_data()
        channels = eegfun.get_cols_by_group(epoch_dat, :channels)
        @test channels == [:Fp1, :Fp2]
        metadata = eegfun.get_cols_by_group(epoch_dat, :metadata)
        @test metadata == [:time, :triggers]
        extra = eegfun.get_cols_by_group(epoch_dat, :extra)
        @test extra == [:vEOG]

        # Test with ErpData 
        erp_dat = create_test_erp_data()
        channels = eegfun.get_cols_by_group(erp_dat, :channels)
        @test channels == [:Fp1, :Fp2]
        metadata = eegfun.get_cols_by_group(erp_dat, :metadata)
        @test metadata == [:time]
        extra = eegfun.get_cols_by_group(erp_dat, :extra)
        @test extra == Symbol[]

        # Test edge cases
        empty_layout = eegfun.Layout(DataFrame(label = Symbol[]), nothing, nothing)
        empty_dat = eegfun.ContinuousData(DataFrame(time = Float64[]), empty_layout, 1000, eegfun.AnalysisInfo())
        @test eegfun.get_cols_by_group(empty_dat, :channels) == Symbol[]
        @test eegfun.get_cols_by_group(empty_dat, :metadata) == Symbol[]
        @test eegfun.get_cols_by_group(empty_dat, :extra) == Symbol[]
        
        # Test invalid group (should return nothing and log error)
        result = eegfun.get_cols_by_group(dat, :invalid_group)
        @test result === nothing
    end

    # === DATA ACCESS FUNCTION TESTS ===
    @testset "Data Access Functions" begin
        
        @testset "all_data and all_labels" begin
            # Test ContinuousData
            cont_dat = create_test_continuous_data()
            all_data_result = eegfun.all_data(cont_dat)
            @test all_data_result isa DataFrame
            @test size(all_data_result) == (10, 8)
            @test propertynames(all_data_result) == [:time, :triggers, :Fp1, :Fp2, :F3, :F4, :vEOG, :hEOG]

            all_labels_result = eegfun.all_labels(cont_dat)
            @test all_labels_result == [:time, :triggers, :Fp1, :Fp2, :F3, :F4, :vEOG, :hEOG]

            # Test EpochData
            epoch_dat = create_test_epoch_data()
            all_data_result = eegfun.all_data(epoch_dat)
            @test all_data_result isa DataFrame
            @test size(all_data_result) == (6, 5)  # 2 epochs * 3 rows each

            all_labels_result = eegfun.all_labels(epoch_dat)
            @test all_labels_result == [:time, :triggers, :Fp1, :Fp2, :vEOG]

            # Test all_labels with specific epoch
            epoch_labels = eegfun.all_labels(epoch_dat, 1)
            @test epoch_labels == [:time, :triggers, :Fp1, :Fp2, :vEOG]

            # Test ErpData
            erp_dat = create_test_erp_data()
            all_data_result = eegfun.all_data(erp_dat)
            @test all_data_result isa DataFrame
            @test size(all_data_result) == (4, 3)
            
            # Test all_labels with DataFrame
            df = DataFrame(a = [1, 2], b = [3, 4])
            @test eegfun.all_labels(df) == [:a, :b]
        end

        @testset "meta_labels and meta_data" begin
            cont_dat = create_test_continuous_data()
            
            # Test meta_labels
            meta_labels = eegfun.meta_labels(cont_dat)
            @test meta_labels == [:time, :triggers]

            # Test meta_data for ContinuousData
            meta_data = eegfun.meta_data(cont_dat)
            @test meta_data isa DataFrame
            @test propertynames(meta_data) == [:time, :triggers]
            @test size(meta_data) == (10, 2)

            # Test meta_data for EpochData
            epoch_dat = create_test_epoch_data()
            meta_data = eegfun.meta_data(epoch_dat)
            @test meta_data isa DataFrame
            @test propertynames(meta_data) == [:time, :triggers]
            @test size(meta_data) == (6, 2)  # 2 epochs * 3 rows each

            # Test meta_data for specific epoch
            meta_data_epoch = eegfun.meta_data(epoch_dat, 1)
            @test meta_data_epoch isa DataFrame
            @test size(meta_data_epoch) == (3, 2)
        end

        @testset "channel_labels and channel_data" begin
            cont_dat = create_test_continuous_data()
            
            # Test channel_labels
            chan_labels = eegfun.channel_labels(cont_dat)
            @test chan_labels == [:Fp1, :Fp2, :F3, :F4]

            # Test channel_labels with indices
            chan_labels_subset = eegfun.channel_labels(cont_dat, [1, 3])
            @test chan_labels_subset == [:Fp1, :F3]

            chan_labels_single = eegfun.channel_labels(cont_dat, 2)
            @test chan_labels_single == [:Fp2]

            chan_labels_range = eegfun.channel_labels(cont_dat, [1:2])
            @test chan_labels_range == [:Fp1, :Fp2]
            
            # Test with invalid indices (should error)
            @test_throws BoundsError eegfun.channel_labels(cont_dat, [10])  # channel 10 doesn't exist
            @test_throws BoundsError eegfun.channel_labels(cont_dat, [1:10])  # range includes non-existent channels

            # Test channel_data for ContinuousData
            chan_data = eegfun.channel_data(cont_dat)
            @test chan_data isa DataFrame
            @test propertynames(chan_data) == [:Fp1, :Fp2, :F3, :F4]
            @test size(chan_data) == (10, 4)

            # Test channel_data for EpochData
            epoch_dat = create_test_epoch_data()
            chan_data = eegfun.channel_data(epoch_dat)
            @test chan_data isa DataFrame
            @test propertynames(chan_data) == [:Fp1, :Fp2]
            @test size(chan_data) == (6, 2)

            # Test channel_data for specific epoch
            chan_data_epoch = eegfun.channel_data(epoch_dat, 1)
            @test chan_data_epoch isa DataFrame
            @test size(chan_data_epoch) == (3, 2)
        end

        @testset "extra_labels and extra_data" begin
            cont_dat = create_test_continuous_data()
            
            # Test extra_labels
            extra_labels = eegfun.extra_labels(cont_dat)
            @test extra_labels == [:vEOG, :hEOG]

            # Test extra_data for ContinuousData
            extra_data = eegfun.extra_data(cont_dat)
            @test extra_data isa DataFrame
            @test propertynames(extra_data) == [:vEOG, :hEOG]
            @test size(extra_data) == (10, 2)

            # Test extra_data for EpochData
            epoch_dat = create_test_epoch_data()
            extra_data = eegfun.extra_data(epoch_dat)
            @test extra_data isa DataFrame
            @test propertynames(extra_data) == [:vEOG]
            @test size(extra_data) == (6, 1)

            # Test extra_data for specific epoch
            extra_data_epoch = eegfun.extra_data(epoch_dat, 1)
            @test extra_data_epoch isa DataFrame
            @test size(extra_data_epoch) == (3, 1)

            # Test ErpData with no extra columns
            erp_dat = create_test_erp_data()
            extra_labels = eegfun.extra_labels(erp_dat)
            @test extra_labels == Symbol[]
            extra_data = eegfun.extra_data(erp_dat)
            @test size(extra_data) == (0, 0)
        end
    end

    # === CONVENIENCE FUNCTION TESTS ===
    @testset "Convenience Functions" begin
        
        @testset "Basic Information Functions" begin
            cont_dat = create_test_continuous_data()
            epoch_dat = create_test_epoch_data()
            erp_dat = create_test_erp_data()

            # Test sample_rate
            @test eegfun.sample_rate(cont_dat) == 1000
            @test eegfun.sample_rate(epoch_dat) == 1000
            @test eegfun.sample_rate(erp_dat) == 1000

            # Test sample_rate from DataFrame
            df_with_time = DataFrame(time = [0.0, 0.001, 0.002, 0.003])
            @test eegfun.sample_rate(df_with_time) == 1000

            # Test reference
            @test eegfun.reference(cont_dat) == :avg
            @test eegfun.reference(cont_dat.analysis_info) == :avg
            @test eegfun.reference(epoch_dat) == :none  # default
            @test eegfun.reference(erp_dat) == :none

            # Test filter_info
            filter_info = eegfun.filter_info(cont_dat.analysis_info)
            @test filter_info == [0.1, 30.0]
            filter_info_default = eegfun.filter_info(epoch_dat.analysis_info)
            @test filter_info_default == [0.0, 0.0]  # default values
        end

        @testset "Size Functions" begin
            cont_dat = create_test_continuous_data()
            epoch_dat = create_test_epoch_data()
            erp_dat = create_test_erp_data()

            # Test n_samples
            @test eegfun.n_samples(cont_dat) == 10
            @test eegfun.n_samples(epoch_dat) == 3  # samples per epoch
            @test eegfun.n_samples(epoch_dat, 1) == 3
            @test eegfun.n_samples(epoch_dat, 2) == 3
            @test eegfun.n_samples(erp_dat) == 4

            # Test n_samples from DataFrame
            @test eegfun.n_samples(cont_dat.data) == 10

            # Test n_channels
            @test eegfun.n_channels(cont_dat) == 4
            @test eegfun.n_channels(epoch_dat) == 2
            @test eegfun.n_channels(erp_dat) == 2

            # Test n_epochs
            @test eegfun.n_epochs(cont_dat) == 1
            @test eegfun.n_epochs(epoch_dat) == 2
            @test eegfun.n_epochs(erp_dat) == 50  # as set in create_test_erp_data

            # Test n_layout
            @test eegfun.n_layout(cont_dat.layout) == 4
            @test eegfun.n_layout(epoch_dat.layout) == 2
        end

        @testset "Duration Functions" begin
            cont_dat = create_test_continuous_data()
            epoch_dat = create_test_epoch_data()

            # Test duration for ContinuousData
            duration_cont = eegfun.duration(cont_dat)
            @test duration_cont ≈ 0.009  # (9 - 0) / 1000

            # Test duration for EpochData
            duration_epoch = eegfun.duration(epoch_dat)
            @test duration_epoch ≈ 0.2  # 0.3 - 0.1

            # Test duration for specific epoch
            duration_epoch1 = eegfun.duration(epoch_dat, 1)
            @test duration_epoch1 ≈ 0.2
            duration_epoch2 = eegfun.duration(epoch_dat, 2)
            @test duration_epoch2 ≈ 0.2

            # Test empty data duration
            empty_df = DataFrame(time = Float64[])
            empty_layout = eegfun.Layout(DataFrame(label = Symbol[]), nothing, nothing)
            empty_dat = eegfun.ContinuousData(empty_df, empty_layout, 1000, eegfun.AnalysisInfo())
            @test eegfun.duration(empty_dat) == 0.0
            
            # Test data without time column
            no_time_df = DataFrame(Fp1 = [1, 2, 3], Fp2 = [4, 5, 6])
            no_time_layout = eegfun.Layout(DataFrame(label = [:Fp1, :Fp2]), nothing, nothing)
            no_time_dat = eegfun.ContinuousData(no_time_df, no_time_layout, 1000, eegfun.AnalysisInfo())
            @test eegfun.duration(no_time_dat) == 0.0
        end
    end

    # === UTILITY FUNCTION TESTS ===
    @testset "Utility Functions" begin
        
        @testset "Channel Utilities" begin
            cont_dat = create_test_continuous_data()
            epoch_dat = create_test_epoch_data()

            # Test has_channels
            @test eegfun.has_channels(cont_dat, [:Fp1, :Fp2])
            @test eegfun.has_channels(cont_dat, [:Fp1])
            @test !eegfun.has_channels(cont_dat, [:C3, :C4])  # not present
            @test !eegfun.has_channels(cont_dat, [:Fp1, :C3])  # partial match

            @test eegfun.has_channels(epoch_dat, [:Fp1, :Fp2])
            @test !eegfun.has_channels(epoch_dat, [:F3, :F4])

            # Test common_channels
            common = eegfun.common_channels(cont_dat, epoch_dat)
            @test common == [:Fp1, :Fp2]

            # Test with same data
            common_same = eegfun.common_channels(cont_dat, cont_dat)
            @test common_same == [:Fp1, :Fp2, :F3, :F4]

            # Test with no common channels
            erp_dat = create_test_erp_data()
            common_subset = eegfun.common_channels(cont_dat, erp_dat)
            @test common_subset == [:Fp1, :Fp2]
        end

        @testset "EpochData Convenience Functions" begin
            # Create epoch data with condition info
            epoch1 = DataFrame(
                time = [0.1, 0.2, 0.3],
                triggers = [1, 0, 0],
                condition = [1, 1, 1],
                condition_name = ["test_condition", "test_condition", "test_condition"],
                file = ["test_file.bdf", "test_file.bdf", "test_file.bdf"],
                Fp1 = [1, 2, 3],
                Fp2 = [4, 5, 6]
            )
            epoch2 = DataFrame(
                time = [0.1, 0.2, 0.3], 
                triggers = [1, 0, 0],
                condition = [1, 1, 1],
                condition_name = ["test_condition", "test_condition", "test_condition"],
                file = ["test_file.bdf", "test_file.bdf", "test_file.bdf"],
                Fp1 = [10, 20, 30],
                Fp2 = [40, 50, 60]
            )
            layout = eegfun.Layout(
                DataFrame(label = [:Fp1, :Fp2], inc = [30, 30], azi = [0, 90]), 
                nothing, nothing
            )
            analysis_info = eegfun.AnalysisInfo()
            epoch_dat = eegfun.EpochData([epoch1, epoch2], layout, 1000, analysis_info)

            # Test condition_number
            @test eegfun.condition_number(epoch_dat) == 1

            # Test condition_name
            @test eegfun.condition_name(epoch_dat) == "test_condition"

            # Test file_name
            @test eegfun.file_name(epoch_dat) == "test_file.bdf"
        end

        @testset "Selection Function Helpers" begin
            # Test component selection functions
            @test eegfun.components()(1:5) == [true, true, true, true, true]
            @test eegfun.components([1, 3])(1:5) == [true, false, true, false, false]
            @test eegfun.components(2)(1:5) == [false, true, false, false, false]
            @test eegfun.components_not([1, 3])(1:5) == [false, true, false, true, true]
            @test eegfun.components_not(2)(1:5) == [true, false, true, true, true]

            # Test get_selected_components (requires mock ICA result)
            # We'll skip this as it requires ICA infrastructure
        end

        @testset "Convert Function" begin
            epoch_dat = create_test_epoch_data()
            
            # Test valid conversion
            converted = eegfun.convert(epoch_dat, 1)
            @test converted isa eegfun.ContinuousData
            @test eegfun.n_samples(converted) == 3
            @test eegfun.n_channels(converted) == 2
            @test converted.data.Fp1 == [1, 2, 3]
            
            # Test second epoch
            converted2 = eegfun.convert(epoch_dat, 2)
            @test converted2.data.Fp1 == [10, 20, 30]
            
            # Test invalid epoch index
            @test_throws Exception eegfun.convert(epoch_dat, 0)
            @test_throws Exception eegfun.convert(epoch_dat, 3)  # only 2 epochs
        end
    end

    # Test viewer, head and tail
    @testset "viewer, head and tail" begin
        eeg = create_test_continuous_data()
        
        # Test viewer (basic functionality - should not error)
        @test nothing === eegfun.viewer(eeg)  # viewer returns nothing
        
        # Test viewer with DataFrame
        df = DataFrame(a = [1, 2], b = [3, 4])
        @test nothing === eegfun.viewer(df)
        eeg = create_test_continuous_data()

        # Test head
        head_result = eegfun.head(eeg)  # Default n=5
        @test head_result isa DataFrame
        @test size(head_result) == (5, 8)
        @test propertynames(head_result) == [:time, :triggers, :Fp1, :Fp2, :F3, :F4, :vEOG, :hEOG]

        head_result_3 = eegfun.head(eeg, n = 3)
        @test head_result_3 isa DataFrame
        @test size(head_result_3) == (3, 8)

        # Test tail
        tail_result = eegfun.tail(eeg)  # Default n=5
        @test tail_result isa DataFrame
        @test size(tail_result) == (5, 8)

        tail_result_3 = eegfun.tail(eeg, n = 3)
        @test tail_result_3 isa DataFrame
        @test size(tail_result_3) == (3, 8)
        
        # Test with very small data
        small_df = DataFrame(time = [1.0], Fp1 = [10.0])
        small_layout = eegfun.Layout(DataFrame(label = [:Fp1]), nothing, nothing)
        small_dat = eegfun.ContinuousData(small_df, small_layout, 1000, eegfun.AnalysisInfo())
        
        # Request more rows than available
        head_small = eegfun.head(small_dat, n = 10)
        @test size(head_small) == (1, 2)  # Should return all available rows
        
        tail_small = eegfun.tail(small_dat, n = 10)
        @test size(tail_small) == (1, 2)  # Should return all available rows
    end

    # Test datarange
    @testset "datarange" begin
        @test eegfun.datarange([1.0, 2.0, 3.0]) == 2.0
        @test eegfun.datarange([-1.0, 0.0, 1.0]) == 2.0
        @test eegfun.datarange([1.0]) == 0.0
    end

    # Test colmeans
    @testset "colmeans" begin
        # Test DataFrame
        df = DataFrame(a = [1.0, 2.0, 3.0], b = [4.0, 5.0, 6.0])
        @test eegfun.colmeans(df, [:a]) ≈ [1.0, 2.0, 3.0]  # Each row's mean of column a
        @test eegfun.colmeans(df, [:b]) ≈ [4.0, 5.0, 6.0]  # Each row's mean of column b
        @test eegfun.colmeans(df, [:a, :b]) ≈ [2.5, 3.5, 4.5]  # Each row's mean across both columns

        # Test Matrix
        mat = [1.0 2.0; 3.0 4.0; 5.0 6.0]  # 3×2 matrix
        @test eegfun.colmeans(mat) ≈ [1.5, 3.5, 5.5]  # Each row's mean across both columns
        @test eegfun.colmeans(mat, [1]) ≈ [1.0, 3.0, 5.0]  # Each row's mean of first column
        @test eegfun.colmeans(mat, [2]) ≈ [2.0, 4.0, 6.0]  # Each row's mean of second column
        
        # Test with empty data
        empty_df = DataFrame(a = Float64[], b = Float64[])
        empty_result = eegfun.colmeans(empty_df, [:a])
        @test isempty(empty_result) || all(isnan, empty_result)  # Should return empty or NaN for empty DataFrame
        
        empty_mat = Matrix{Float64}(undef, 0, 2)
        empty_mat_result = eegfun.colmeans(empty_mat)
        @test isempty(empty_mat_result) || all(isnan, empty_mat_result)  # Should return empty or NaN for empty Matrix
    end

    # Test data_limits
    @testset "data_limits" begin
        df = DataFrame(time = [1.0, 2.0, 3.0], value = [4.0, 5.0, 6.0])

        # Test data_limits_x
        @test eegfun.data_limits_x(df) == (1.0, 3.0)
        @test eegfun.data_limits_x(df, col = :value) == (4.0, 6.0)
        
        # Test with non-existent column (should error)
        @test_throws ArgumentError eegfun.data_limits_x(df, col = :nonexistent)

        # Test data_limits_y
        @test eegfun.data_limits_y(df, :value) == [4.0, 6.0]
        @test eegfun.data_limits_y(df, :time) == [1.0, 3.0]  # Single column
        
        # Test data_limits_y with multiple columns
        @test eegfun.data_limits_y(df, [:time, :value]) == [1.0, 6.0]  # min across all columns, max across all columns
        
        # Test with non-existent column (should error)
        @test_throws ArgumentError eegfun.data_limits_y(df, :nonexistent)
        @test_throws ArgumentError eegfun.data_limits_y(df, [:time, :nonexistent])

        # Test empty data
        empty_df = DataFrame(time = Float64[], value = Float64[])
        @test eegfun.data_limits_x(empty_df) === nothing
        @test eegfun.data_limits_y(empty_df, :value) === nothing
        @test eegfun.data_limits_y(empty_df, [:time, :value]) === nothing
    end

    # Test to_data_frame
    @testset "to_data_frame" begin
        # Use our improved test data
        epoch_data = create_test_epoch_data()

        # Test single EpochData
        result = eegfun.to_data_frame(epoch_data)
        @test size(result) == (6, 5)  # 2 epochs * 3 rows each
        @test result.time == [0.1, 0.2, 0.3, 0.1, 0.2, 0.3]
        @test result.Fp1 == [1, 2, 3, 10, 20, 30]

        # Test Vector of EpochData
        epoch_data_vec = [epoch_data, epoch_data]
        result = eegfun.to_data_frame(epoch_data_vec)
        @test size(result) == (12, 5)  # 2 * (2 epochs * 3 rows each)

        # Test empty EpochData
        empty_layout = eegfun.Layout(DataFrame(label = Symbol[]), nothing, nothing)
        empty_epoch = eegfun.EpochData(DataFrame[], empty_layout, 1000, eegfun.AnalysisInfo())
        @test size(eegfun.to_data_frame(empty_epoch)) == (0, 0)
        @test size(eegfun.to_data_frame([empty_epoch])) == (0, 0)
    end

    # Test ylimits functions
    @testset "ylimits Functions" begin
        
        @testset "ylimits for ErpData" begin
            erp_dat = create_test_erp_data()
            
            # Test basic ylimits
            limits = eegfun.ylimits(erp_dat)
            @test limits isa Tuple{Float64, Float64}
            @test limits[1] < limits[2]  # min < max
            @test limits[1] == 1.0  # minimum value in data
            @test limits[2] == 8.0  # maximum value in data
            
            # Test with channel selection
            limits_fp1 = eegfun.ylimits(erp_dat, channel_selection = eegfun.channels([:Fp1]))
            @test limits_fp1 == (1.0, 4.0)  # Fp1 data range
            
            # Test with sample selection
            limits_positive = eegfun.ylimits(erp_dat, sample_selection = x -> x.time .>= 0.0)
            @test limits_positive[1] >= 2.0  # Should exclude negative time samples
        end
        
        @testset "ylimits for EpochData" begin
            epoch_dat = create_test_epoch_data()
            
            # Test basic ylimits
            limits = eegfun.ylimits(epoch_dat)
            @test limits isa Tuple{Float64, Float64}
            @test limits[1] < limits[2]
            @test limits[1] == 1.0  # minimum value across all epochs
            @test limits[2] == 60.0  # maximum value across all epochs (Fp2 in epoch 2)
            
            # Test with channel selection
            limits_fp2 = eegfun.ylimits(epoch_dat, channel_selection = eegfun.channels([:Fp2]))
            @test limits_fp2[1] >= 4.0  # Fp2 has higher values
            
            # Test with sample selection
            limits_triggers = eegfun.ylimits(epoch_dat, sample_selection = x -> x.triggers .== 1)
            @test limits_triggers isa Tuple{Float64, Float64}
        end
    end

    # Test log_pretty_table (only function remaining in data.jl)
    @testset "log_pretty_table" begin
        df = DataFrame(a = [1, 2], b = [3, 4])
        @test eegfun.log_pretty_table(df) === nothing
    end

    # === SUBSET FUNCTION TESTS ===
    @testset "Subset Functions" begin
        
        @testset "subset_dataframe and subset_dataframes" begin
            # Test basic DataFrame subsetting
            df = DataFrame(
                time = [1, 2, 3, 4],
                Fp1 = [10, 20, 30, 40],
                Fp2 = [100, 200, 300, 400]
            )
            
            result = eegfun.subset_dataframe(df, [:time, :Fp1], [1, 3])
            @test size(result) == (2, 2)
            @test result.time == [1, 3]
            @test result.Fp1 == [10, 30]
            
            # Test with empty DataFrame
            empty_df = DataFrame(time = Float64[], Fp1 = Float64[])
            empty_result = eegfun.subset_dataframe(empty_df, [:time, :Fp1], Int[])
            @test size(empty_result) == (0, 2)
            
            # Test subset_dataframes with Vector{DataFrame}
            df1 = DataFrame(time = [1, 2], Fp1 = [10, 20])
            df2 = DataFrame(time = [3, 4], Fp1 = [30, 40])
            df3 = DataFrame(time = [5, 6], Fp1 = [50, 60])
            
            dataframes = [df1, df2, df3]
            result = eegfun.subset_dataframes(dataframes, [1, 3], [:time], [1])
            @test length(result) == 2
            @test result[1].time == [1]
            @test result[2].time == [5]
        end

        @testset "subset for ContinuousData" begin
            dat = create_test_continuous_data()
            
            # Test basic subset - all data
            subset_all = eegfun.subset(dat)
            @test eegfun.n_samples(subset_all) == 10
            @test eegfun.n_channels(subset_all) == 4
            @test typeof(subset_all) == typeof(dat)
            
            # Test channel selection
            subset_chans = eegfun.subset(dat, channel_selection = eegfun.channels([:Fp1, :Fp2]))
            @test eegfun.n_channels(subset_chans) == 2
            @test eegfun.channel_labels(subset_chans) == [:Fp1, :Fp2]
            
            # Test sample selection
            subset_samples = eegfun.subset(dat, sample_selection = x -> x.triggers .== 1)
            @test eegfun.n_samples(subset_samples) == 3  # 3 trigger events
            
            # Test combined selection
            subset_combined = eegfun.subset(dat, 
                channel_selection = eegfun.channels([:Fp1]), 
                sample_selection = x -> x.triggers .== 1)
            @test eegfun.n_channels(subset_combined) == 1
            @test eegfun.n_samples(subset_combined) == 3
            
            # Test include_extra
            subset_extra = eegfun.subset(dat, 
                channel_selection = eegfun.channels([:vEOG]), 
                include_extra = true)
            @test :vEOG in eegfun.all_labels(subset_extra)
        end

        @testset "subset for EpochData" begin
            dat = create_test_epoch_data()
            
            # Test basic subset - all data
            subset_all = eegfun.subset(dat)
            @test eegfun.n_epochs(subset_all) == 2
            @test eegfun.n_channels(subset_all) == 2
            @test typeof(subset_all) == typeof(dat)
            
            # Test epoch selection
            subset_epoch = eegfun.subset(dat, epoch_selection = eegfun.epochs([1]))
            @test eegfun.n_epochs(subset_epoch) == 1
            @test length(subset_epoch.data) == 1
            
            # Test channel selection
            subset_chans = eegfun.subset(dat, channel_selection = eegfun.channels([:Fp1]))
            @test eegfun.n_channels(subset_chans) == 1
            @test eegfun.channel_labels(subset_chans) == [:Fp1]
            
            # Test sample selection
            subset_samples = eegfun.subset(dat, sample_selection = x -> x.triggers .== 1)
            @test eegfun.n_samples(subset_samples, 1) == 1  # Only first sample has trigger
            
            # Test combined selection
            subset_combined = eegfun.subset(dat, 
                epoch_selection = eegfun.epochs([2]),
                channel_selection = eegfun.channels([:Fp2]),
                sample_selection = x -> x.triggers .== 0)
            @test eegfun.n_epochs(subset_combined) == 1
            @test eegfun.n_channels(subset_combined) == 1
            @test eegfun.n_samples(subset_combined, 1) == 2  # 2 non-trigger samples
            
            # Test include_extra
            subset_extra = eegfun.subset(dat, 
                channel_selection = eegfun.channels([:vEOG]), 
                include_extra = true)
            @test :vEOG in eegfun.all_labels(subset_extra)
        end

        @testset "subset for ErpData" begin
            dat = create_test_erp_data()
            
            # Test basic subset - all data
            subset_all = eegfun.subset(dat)
            @test eegfun.n_samples(subset_all) == 4
            @test eegfun.n_channels(subset_all) == 2
            @test eegfun.n_epochs(subset_all) == 50  # Preserved n_epochs
            @test typeof(subset_all) == typeof(dat)
            
            # Test channel selection
            subset_chans = eegfun.subset(dat, channel_selection = eegfun.channels([:Fp1]))
            @test eegfun.n_channels(subset_chans) == 1
            @test eegfun.channel_labels(subset_chans) == [:Fp1]
            @test eegfun.n_epochs(subset_chans) == 50
            
            # Test sample selection - select positive time points
            subset_samples = eegfun.subset(dat, sample_selection = x -> x.time .>= 0.0)
            @test eegfun.n_samples(subset_samples) == 3  # times: 0.0, 0.1, 0.2
            @test eegfun.n_epochs(subset_samples) == 50
            
            # Test combined selection
            subset_combined = eegfun.subset(dat, 
                channel_selection = eegfun.channels([:Fp2]),
                sample_selection = x -> x.time .>= 0.0)
            @test eegfun.n_channels(subset_combined) == 1
            @test eegfun.n_samples(subset_combined) == 3
            @test eegfun.n_epochs(subset_combined) == 50
        end

        @testset "Helper Functions" begin
            # Note: _subset_common and _create_subset are internal functions
            # but we can test them indirectly through the main subset functions
            
            # Test that subsetting preserves data structure integrity
            cont_dat = create_test_continuous_data()
            subset_result = eegfun.subset(cont_dat, channel_selection = eegfun.channels([:Fp1, :Fp2]))
            
            # Verify metadata is preserved
            @test subset_result.sample_rate == cont_dat.sample_rate
            @test subset_result.analysis_info.reference == cont_dat.analysis_info.reference
            @test subset_result.analysis_info.hp_filter == cont_dat.analysis_info.hp_filter
            @test subset_result.analysis_info.lp_filter == cont_dat.analysis_info.lp_filter
            
            # Verify layout is properly subset
            @test length(eegfun.channel_labels(subset_result.layout)) == 2
            @test eegfun.channel_labels(subset_result.layout) == [:Fp1, :Fp2]
            
            # Test with EpochData
            epoch_dat = create_test_epoch_data()
            subset_result = eegfun.subset(epoch_dat, epoch_selection = eegfun.epochs([1]))
            
            # Verify epochs are properly selected
            @test length(subset_result.data) == 1
            @test subset_result.data[1].time == epoch_dat.data[1].time
            @test subset_result.data[1].Fp1 == epoch_dat.data[1].Fp1
        end

        @testset "Edge Cases" begin
            # Test empty selection
            dat = create_test_continuous_data()
            
            # Empty channel selection should result in no channel data, but metadata preserved
            subset_empty_channels = eegfun.subset(dat, channel_selection = eegfun.channels(Symbol[]))
            @test eegfun.n_channels(subset_empty_channels) == 0
            @test :time in eegfun.all_labels(subset_empty_channels)  # metadata preserved
            
            # Empty sample selection
            subset_empty_samples = eegfun.subset(dat, sample_selection = x -> falses(nrow(x)))
            @test eegfun.n_samples(subset_empty_samples) == 0
            
            # Test with EpochData empty selection
            epoch_dat = create_test_epoch_data()
            subset_empty_epochs = eegfun.subset(epoch_dat, epoch_selection = eegfun.epochs(Int[]))
            @test eegfun.n_epochs(subset_empty_epochs) == 0
            @test length(subset_empty_epochs.data) == 0
            
            # Test with invalid epoch indices (should return empty EpochData)
            invalid_epoch_result = eegfun.subset(epoch_dat, epoch_selection = eegfun.epochs([10]))  # epoch 10 doesn't exist
            @test eegfun.n_epochs(invalid_epoch_result) == 0
            @test length(invalid_epoch_result.data) == 0
        end

    # Test rename_channel functionality with EEG data types
    @testset "rename_channel with EEG data" begin
        # Test basic renaming with ContinuousData
        dat_cont = create_test_continuous_data()
        dat_rename = copy(dat_cont)
        eegfun.rename_channel!(dat_rename, Dict(:Fp1 => :Fp1_new))
        @test :Fp1_new ∈ propertynames(dat_rename.data)
        @test :Fp1 ∉ propertynames(dat_rename.data)
        @test :Fp1_new in dat_rename.layout.data.label
        
        # Test swap behavior with ContinuousData
        dat_swap = copy(dat_cont)
        eegfun.rename_channel!(dat_swap, Dict(:Fp1 => :Fp2, :Fp2 => :Fp1))
        @test :Fp1 in dat_swap.layout.data.label
        @test :Fp2 in dat_swap.layout.data.label
        @test :Fp1 in propertynames(dat_swap.data)
        @test :Fp2 in propertynames(dat_swap.data)
        
        # Test that duplicate names are prevented
        @test_throws Any eegfun.rename_channel!(copy(dat_cont), Dict(:Fp1 => :X, :Fp2 => :X))
        
        # Test non-mutating version with ContinuousData
        dat_copy = copy(dat_cont)
        new_dat = eegfun.rename_channel(dat_copy, Dict(:Fp1 => :Fp1_new))
        @test :Fp1_new ∈ propertynames(new_dat.data)
        @test :Fp1 ∉ propertynames(new_dat.data)
        @test :Fp1 ∈ propertynames(dat_copy.data)  # Original unchanged
        @test :Fp1_new ∉ propertynames(dat_copy.data)  # Original unchanged
        
        # Test with EpochData
        dat_epoch = create_test_epoch_data()
        eegfun.rename_channel!(dat_epoch, Dict(:Fp1 => :Fp1_new))
        @test :Fp1_new ∈ propertynames(dat_epoch.data[1])
        @test :Fp1_new ∈ propertynames(dat_epoch.data[2])
        @test :Fp1 ∉ propertynames(dat_epoch.data[1])
        @test :Fp1 ∉ propertynames(dat_epoch.data[2])
        
        # Test with ErpData
        dat_erp = create_test_erp_data()
        eegfun.rename_channel!(dat_erp, Dict(:Fp1 => :Fp1_new))
        @test :Fp1_new ∈ propertynames(dat_erp.data)
        @test :Fp1 ∉ propertynames(dat_erp.data)
    end


    end
end
