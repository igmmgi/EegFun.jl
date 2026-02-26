using Test
using Dates

@testset "Logging Utilities" begin
    # Create temporary directory for log files
    test_dir = mktempdir()

    @testset "Duration formatting" begin
        # Test _format_duration with different durations
        @test EegFun._format_duration(Millisecond(5000)) == "5 seconds"
        @test EegFun._format_duration(Millisecond(30000)) == "30 seconds"
        @test EegFun._format_duration(Millisecond(90000)) == "1 minute, 30 seconds"
        @test EegFun._format_duration(Millisecond(3661000)) == "1 hour, 1 minute, 1 second"
        @test EegFun._format_duration(Millisecond(7200000)) == "2 hours"
    end

    @testset "Global logging setup and teardown" begin
        log_file = joinpath(test_dir, "global_test.log")

        # Test setup_global_logging
        global_log_handle = EegFun.setup_global_logging(log_file)
        @test global_log_handle isa IO
        @test isfile(log_file)

        # Test that logging works (basic functionality)
        @info "Test global logging message"

        # Test close_global_logging
        EegFun.close_global_logging()

        # Verify log file exists and has some content
        @test isfile(log_file)
        log_content = read(log_file, String)
        @test !isempty(log_content)
    end

    @testset "File logging setup and teardown" begin
        log_file = joinpath(test_dir, "file_test.log")

        # Test setup_logging
        EegFun.setup_logging(log_file)
        @test isfile(log_file)

        # Test that logging works (basic functionality)
        @info "Test file logging message"

        # Test close_logging
        EegFun.close_logging()

        # Verify log file exists and has some content
        @test isfile(log_file)
        log_content = read(log_file, String)
        @test !isempty(log_content)
    end

    @testset "Logging state management" begin
        # Test that we can switch between global and file logging
        global_log_file = joinpath(test_dir, "global_state.log")
        file_log_file = joinpath(test_dir, "file_state.log")

        # Setup global logging
        EegFun.setup_global_logging(global_log_file)
        @info "Global message 1"

        # Setup file logging (should work alongside global)
        EegFun.setup_logging(file_log_file)
        @info "File message 1"

        # Close file logging
        EegFun.close_logging()
        @info "Global message 2"

        # Close global logging
        EegFun.close_global_logging()

        # Verify both log files exist and have content
        @test isfile(global_log_file)
        @test isfile(file_log_file)

        global_content = read(global_log_file, String)
        file_content = read(file_log_file, String)

        @test contains(global_content, "Global message 1")
        @test contains(global_content, "Global message 2")
        @test contains(file_content, "File message 1")
    end

    @testset "Multiple file logging sessions" begin
        # Test that we can have multiple file logging sessions
        log_file1 = joinpath(test_dir, "session1.log")
        log_file2 = joinpath(test_dir, "session2.log")

        # First session
        EegFun.setup_logging(log_file1)
        @info "Session 1 message"
        EegFun.close_logging()

        # Second session
        EegFun.setup_logging(log_file2)
        @info "Session 2 message"
        EegFun.close_logging()

        # Verify both files exist and have correct content
        @test isfile(log_file1)
        @test isfile(log_file2)

        content1 = read(log_file1, String)
        content2 = read(log_file2, String)

        @test contains(content1, "Session 1 message")
        @test contains(content2, "Session 2 message")
        @test !contains(content1, "Session 2 message")
        @test !contains(content2, "Session 1 message")
    end

    @testset "Error handling" begin
        # Test that closing non-existent logging doesn't cause errors
        EegFun.close_logging()  # Should not throw error
        EegFun.close_global_logging()  # Should not throw error

        # Test that we can call setup multiple times safely
        log_file = joinpath(test_dir, "multiple_setup.log")
        EegFun.setup_logging(log_file)
        EegFun.setup_logging(log_file)  # Should not throw error
        EegFun.close_logging()
    end




    @testset "Duration formatting edge cases" begin
        # Test very short durations
        @test EegFun._format_duration(Millisecond(0)) == "0 seconds"
        @test EegFun._format_duration(Millisecond(1)) == "0 seconds"  # 1ms rounds to 0 seconds

        # Test exactly 1 minute
        @test EegFun._format_duration(Millisecond(60000)) == "1 minute"

        # Test exactly 1 hour
        @test EegFun._format_duration(Millisecond(3600000)) == "1 hour"

        # Test very long duration
        @test EegFun._format_duration(Millisecond(3661000)) == "1 hour, 1 minute, 1 second"
    end

    # Clean up
    rm(test_dir, recursive = true, force = true)
end
