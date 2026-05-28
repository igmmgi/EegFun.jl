using Test
using EegFun

@testset "Read Raw Data" begin
    # Base path for resources
    data_dir = joinpath(dirname(dirname(@__DIR__)), "resources", "data")

    @testset "BDF Import" begin
        test_file = joinpath(data_dir, "bdf", "example1.bdf")
        if isfile(test_file)
            result = EegFun.read_raw_data(test_file)
            @test result isa EegFun.BiosemiDataFormat.BiosemiData
            @test size(result.data, 1) > 0
            @test length(result.time) > 0
            @test result.header.sample_rate[1] > 0
            
            # Test ContinuousData construction
            dat = EegFun.create_eegfun_data(result)
            @test dat isa EegFun.ContinuousData
            @test EegFun.sample_rate(dat) > 0
            @test size(EegFun.all_data(dat), 1) > 0
        else
            EegFun.@minimal_warning "Test file not found: $test_file"
        end
    end

    @testset "BrainVision Import" begin
        test_file = joinpath(data_dir, "brainvision", "example1.vhdr")
        if isfile(test_file)
            result = EegFun.read_raw_data(test_file)
            @test result isa EegFun.BrainVisionDataFormat.BrainVisionData
            @test size(result.data, 1) > 0
            @test result.header.Fs > 0
            
            # Test ContinuousData construction
            dat = EegFun.create_eegfun_data(result)
            @test dat isa EegFun.ContinuousData
            @test EegFun.sample_rate(dat) > 0
            @test size(EegFun.all_data(dat), 1) > 0
        else
            EegFun.@minimal_warning "Test file not found: $test_file"
        end
    end

    @testset "EDF Import" begin
        test_file = joinpath(data_dir, "edf", "test.edf")
        if isfile(test_file)
            result = EegFun.read_raw_data(test_file)
            @test result isa EegFun.EuropeanDataFormat.EdfData
            @test size(result.data, 1) > 0
            @test length(result.time) > 0
            @test result.header.sample_rate[1] > 0
            
            # Test ContinuousData construction
            dat = EegFun.create_eegfun_data(result)
            @test dat isa EegFun.ContinuousData
            @test EegFun.sample_rate(dat) > 0
            @test size(EegFun.all_data(dat), 1) > 0
        else
            EegFun.@minimal_warning "Test file not found: $test_file"
        end

        test_file5 = joinpath(data_dir, "edf", "test5.edf")
        if isfile(test_file5)
            result5 = EegFun.read_raw_data(test_file5)
            @test result5 isa EegFun.EuropeanDataFormat.EdfData
            @test size(result5.data, 1) > 0
            @test result5.header.sample_rate[1] == 200
            
            dat5 = EegFun.create_eegfun_data(result5)
            @test dat5 isa EegFun.ContinuousData
            @test size(EegFun.all_data(dat5), 1) > 0
        else
            EegFun.@minimal_warning "Test file not found: $test_file5"
        end
    end

    @testset "FIF Import" begin
        test_file = joinpath(data_dir, "fif", "test_raw.fif")
        if isfile(test_file)
            result = EegFun.read_raw_data(test_file)
            @test result isa EegFun.FunctionalImageFormat.FifData
            @test size(result.data, 1) > 0
            @test result.header.sfreq > 0
            
            # Test ContinuousData construction
            dat = EegFun.create_eegfun_data(result)
            @test dat isa EegFun.ContinuousData
            @test EegFun.sample_rate(dat) > 0
            @test size(EegFun.all_data(dat), 1) > 0
        else
            EegFun.@minimal_warning "Test file not found: $test_file"
        end
        
        test_epo_file = joinpath(data_dir, "fif", "test_epochs.fif")
        if isfile(test_epo_file)
            result_epo = EegFun.read_raw_data(test_epo_file)
            @test result_epo isa EegFun.FunctionalImageFormat.FifEpochs
            @test size(result_epo.data, 3) > 0  # number of epochs
            
            # Test EpochData construction
            epo_dat = EegFun.create_eegfun_data(result_epo)
            @test epo_dat isa Union{EegFun.EpochData, EegFun.ErpData}
            @test EegFun.sample_rate(epo_dat) > 0
            
            # Verify trigger mapping was successfully extracted
            # Check that the trigger column exists and sum > 0 if triggers are present
            first_epoch_df = epo_dat.data[1]
            @test "trigger" in names(first_epoch_df)
        else
            EegFun.@minimal_warning "Test epoch file not found: $test_epo_file"
        end
    end
    @testset "XDF Import" begin
        test_file = joinpath(data_dir, "xdf", "test1.xdf")
        if isfile(test_file)
            result = EegFun.read_raw_data(test_file)
            @test result isa EegFun.ExtensibleDataFormat.XdfData
            
            eeg_streams = [s for s in values(result.streams) if s.header.type == "EEG"]
            @test !isempty(eeg_streams)
            @test size(eeg_streams[1].time_series, 1) > 0
            
            # Test ContinuousData construction
            dat = EegFun.create_eegfun_data(result)
            @test dat isa EegFun.ContinuousData
            @test EegFun.sample_rate(dat) > 0
            @test size(EegFun.all_data(dat), 1) > 0
            @test "trigger" in names(EegFun.all_data(dat))
            
            # Test trigger_count
            info = EegFun.trigger_count(dat)
            @test info isa EegFun.TriggerInfo
        else
            EegFun.@minimal_warning "Test file not found: $test_file"
        end
    end
end
