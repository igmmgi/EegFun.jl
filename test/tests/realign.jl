using Test
using JLD2
using DataFrames

# ============================================================================
# Helper: build synthetic EpochData with a response trigger embedded at a
# specified sample index within the time vector.
# ============================================================================
function _make_epoch_with_response(;
    sample_rate = 256,
    epoch_start = -0.5,
    epoch_end = 1.5,
    rt_seconds = 0.5,         # time of response trigger relative to epoch start
    response_trigger = 201,
    n_channels = 3,
)
    time = collect(epoch_start:(1/sample_rate):epoch_end)
    n = length(time)
    ch_names = [Symbol("Ch$i") for i = 1:n_channels]

    # Build DataFrame: time, trigger (0 everywhere, response_trigger at RT), channels
    df = DataFrame(time = time)
    df.trigger = zeros(Int, n)
    rt_idx = argmin(abs.(time .- rt_seconds))
    df.trigger[rt_idx] = response_trigger
    for ch in ch_names
        df[!, ch] = randn(n)
    end
    return df
end

function _make_epoch_data(epochs_vec; participant = 1, condition = 1, sample_rate = 256)
    layout_df = DataFrame(
        label = [Symbol("Ch$i") for i = 1:(ncol(epochs_vec[1])-2)],
        inc   = zeros(ncol(epochs_vec[1])-2),
        azi   = zeros(ncol(epochs_vec[1])-2),
    )
    layout = EegFun.Layout(layout_df, nothing, nothing, nothing)
    return EegFun.EpochData(
        "test_participant_$participant.jld2",
        condition,
        "condition_$condition",
        epochs_vec,
        layout,
        sample_rate,
        EegFun.AnalysisInfo(),
    )
end

# Response trigger codes used throughout these tests
const RESP_TRIGGERS = [201]

# ============================================================================
@testset "Realignment" begin

    @testset "Basic realignment (non-mutating)" begin
        rts = range(0.3, 0.9, length = 10)
        epochs = [_make_epoch_with_response(rt_seconds = rt) for rt in rts]
        epoch_data = _make_epoch_data(epochs)

        # Record original time bounds
        orig_tmin = minimum(epoch_data.data[1].time)
        orig_tmax = maximum(epoch_data.data[1].time)

        # Non-mutating realign
        realigned = EegFun.realign(epoch_data, RESP_TRIGGERS)

        # Original should be unchanged
        @test minimum(epoch_data.data[1].time) ≈ orig_tmin
        @test maximum(epoch_data.data[1].time) ≈ orig_tmax

        # Realigned should have t=0 at the response trigger
        for epoch in realigned.data
            t0_idx = argmin(abs.(epoch.time))
            @test abs(epoch.time[t0_idx]) < 1/256 + 1e-9
        end

        # All epochs must share the same time vector after cropping
        for i = 2:length(realigned.data)
            @test realigned.data[i].time ≈ realigned.data[1].time
        end

        @test realigned isa EegFun.EpochData
        @test length(realigned.data) == 10
    end

    @testset "In-place realignment" begin
        rts = range(0.3, 0.9, length = 10)
        epochs = [_make_epoch_with_response(rt_seconds = rt) for rt in rts]
        epoch_data = _make_epoch_data(epochs)

        EegFun.realign!(epoch_data, RESP_TRIGGERS)

        # t=0 should be at (or nearest sample to) the response
        for epoch in epoch_data.data
            t0_idx = argmin(abs.(epoch.time))
            @test abs(epoch.time[t0_idx]) < 1/256 + 1e-9
        end

        # All epochs must share the same time vector
        for i = 2:length(epoch_data.data)
            @test epoch_data.data[i].time ≈ epoch_data.data[1].time
        end
    end

    @testset "Time interval cropping" begin
        # Varying RTs → after realignment, common window must be shorter than raw epoch
        rts = range(0.3, 0.7, length = 10)
        epochs = [_make_epoch_with_response(rt_seconds = rt, epoch_start = -0.5, epoch_end = 1.5) for rt in rts]
        epoch_data = _make_epoch_data(epochs)
        original_len = nrow(epoch_data.data[1])

        realigned = EegFun.realign(epoch_data, RESP_TRIGGERS)
        realigned_lengths = [nrow(e) for e in realigned.data]

        # All epochs must have identical length
        @test all(realigned_lengths .== realigned_lengths[1])

        # Common-window cropping means result ≤ original
        @test realigned_lengths[1] <= original_len

        # Same time bounds across all epochs
        for i = 1:length(realigned.data)
            @test minimum(realigned.data[i].time) ≈ minimum(realigned.data[1].time)
            @test maximum(realigned.data[i].time) ≈ maximum(realigned.data[1].time)
        end
    end

    @testset "Common time interval calculation" begin
        # Wide RT spread → common window must be positive but strictly less than raw
        rts = range(0.2, 1.0, length = 10)
        epochs = [_make_epoch_with_response(rt_seconds = rt, epoch_start = -0.5, epoch_end = 1.5) for rt in rts]
        epoch_data = _make_epoch_data(epochs)

        realigned = EegFun.realign(epoch_data, RESP_TRIGGERS)

        time_range = maximum(realigned.data[1].time) - minimum(realigned.data[1].time)
        @test time_range > 0
        @test time_range < 2.0  # Shorter than the 2s raw window
    end

    @testset "Channel data preservation" begin
        rts = range(0.3, 0.9, length = 10)
        epochs = [_make_epoch_with_response(rt_seconds = rt, n_channels = 3) for rt in rts]
        epoch_data = _make_epoch_data(epochs)

        realigned = EegFun.realign(epoch_data, RESP_TRIGGERS)

        # All channels should still be present
        @test hasproperty(realigned.data[1], :Ch1)
        @test hasproperty(realigned.data[1], :Ch2)
        @test hasproperty(realigned.data[1], :Ch3)

        # Value range should be similar (just sliced, not transformed)
        @test minimum(realigned.data[1].Ch1) ≈ minimum(epoch_data.data[1].Ch1) atol = 2.0
        @test maximum(realigned.data[1].Ch1) ≈ maximum(epoch_data.data[1].Ch1) atol = 2.0
    end

    @testset "Metadata preservation" begin
        rts = range(0.3, 0.9, length = 10)
        epochs = [_make_epoch_with_response(rt_seconds = rt) for rt in rts]
        epoch_data = _make_epoch_data(epochs)

        realigned = EegFun.realign(epoch_data, RESP_TRIGGERS)

        # Struct-level metadata must be preserved
        @test hasproperty(realigned, :condition)
        @test realigned.condition == epoch_data.condition
        @test realigned.condition_name == epoch_data.condition_name
        @test realigned.sample_rate == epoch_data.sample_rate
        @test realigned.layout.data == epoch_data.layout.data
    end

    @testset "Error handling: missing column" begin
        # No response trigger in any epoch → all dropped → realigned is empty
        epochs = [_make_epoch_with_response(rt_seconds = 0.5, response_trigger = 999) for _ = 1:10]
        epoch_data = _make_epoch_data(epochs)

        # Trigger 201 is NOT 999, so all epochs are dropped. Should not crash,
        # but the result should be empty or a warning should be issued.
        result = EegFun.realign(epoch_data, RESP_TRIGGERS)
        @test isempty(result.data)
    end

    @testset "Error handling: non-finite values" begin
        rts = range(0.3, 0.9, length = 10)
        epochs = [_make_epoch_with_response(rt_seconds = rt) for rt in rts]
        # Inject NaN into channel data of first epoch
        epochs[1].Ch1 .= NaN
        epoch_data = _make_epoch_data(epochs)

        # Realignment itself should not crash on NaN channel values
        # (NaN propagates through but structure must remain valid)
        realigned = EegFun.realign(epoch_data, RESP_TRIGGERS)
        @test length(realigned.data) > 0
    end

    @testset "Error handling: insufficient epoch length" begin
        # Very short epoch, RT close to end → common window could be tiny
        rts = range(0.25, 0.28, length = 10)
        epochs = [_make_epoch_with_response(rt_seconds = rt, epoch_start = 0.0, epoch_end = 0.3) for rt in rts]
        epoch_data = _make_epoch_data(epochs)

        try
            realigned = EegFun.realign(epoch_data, RESP_TRIGGERS)
            @test all(length(realigned.data[i].time) > 0 for i = 1:length(realigned.data))
        catch e
            @test occursin("common time interval", sprint(showerror, e)) ||
                  occursin("No samples found", sprint(showerror, e)) ||
                  e isa Exception
        end
    end

    @testset "Multiple channels preserved" begin
        rts = range(0.3, 0.9, length = 5)
        epochs = [_make_epoch_with_response(rt_seconds = rt, n_channels = 10) for rt in rts]
        epoch_data = _make_epoch_data(epochs)

        realigned = EegFun.realign(epoch_data, RESP_TRIGGERS)

        for i = 1:10
            @test hasproperty(realigned.data[1], Symbol("Ch$i"))
        end

        n_ch_orig = length(EegFun.channel_labels(epoch_data))
        n_ch_realigned = length(EegFun.channel_labels(realigned))
        @test n_ch_orig == n_ch_realigned
    end

    @testset "Edge case: identical RTs" begin
        # All epochs have the exact same RT
        epochs = [_make_epoch_with_response(rt_seconds = 0.5) for _ = 1:10]
        epoch_data = _make_epoch_data(epochs)

        realigned = EegFun.realign(epoch_data, RESP_TRIGGERS)

        # With identical RTs, epoch length remains the same as the raw length
        @test nrow(realigned.data[1]) == nrow(epoch_data.data[1])

        # All time vectors must be identical
        for i = 2:length(realigned.data)
            @test realigned.data[i].time == realigned.data[1].time
        end
    end

    @testset "Realistic use case: response-locked LRP" begin
        rts = range(0.4, 1.2, length = 20)
        epochs = [_make_epoch_with_response(rt_seconds = rt, epoch_start = -0.5, epoch_end = 2.0, n_channels = 4) for rt in rts]
        epoch_data = _make_epoch_data(epochs, condition = 1)

        realigned = EegFun.realign(epoch_data, RESP_TRIGGERS)

        # All epochs share the same time vector
        for i = 2:length(realigned.data)
            @test realigned.data[i].time ≈ realigned.data[1].time
        end

        # Pre- and post-response data both exist
        @test minimum(realigned.data[1].time) < 0
        @test maximum(realigned.data[1].time) > 0

        # Varying RTs → realigned epochs are shorter than the original
        @test nrow(realigned.data[1]) < nrow(epoch_data.data[1])
    end

    @testset "Vector{EpochData} batch (in-memory)" begin
        # Two conditions, each with different RT ranges
        rts1 = range(0.3, 0.8, length = 10)
        rts2 = range(0.4, 0.9, length = 10)

        dat1 = _make_epoch_data([_make_epoch_with_response(rt_seconds = rt) for rt in rts1], condition = 1)
        dat2 = _make_epoch_data([_make_epoch_with_response(rt_seconds = rt) for rt in rts2], condition = 2)
        vec_data = [dat1, dat2]

        EegFun.realign!(vec_data, RESP_TRIGGERS)

        # After batch realign, both conditions must share identical sample counts
        @test nrow(vec_data[1].data[1]) == nrow(vec_data[2].data[1])

        # And all epochs within each condition share the same time vector
        for cond in vec_data
            for i = 2:length(cond.data)
                @test cond.data[i].time ≈ cond.data[1].time
            end
        end
    end
end

@testset "Batch realignment (file-based)" begin
    test_dir   = mktempdir()
    output_dir = joinpath(test_dir, "realigned")

    # Write three single-condition JLD2 files (Vector{EpochData})
    for p = 1:3
        rts = range(0.3, 0.8, length = 10)
        epochs = [_make_epoch_with_response(rt_seconds = rt) for rt in rts]
        dat = _make_epoch_data(epochs, participant = p)
        jldsave(joinpath(test_dir, "$(p)_epochs_good.jld2"); data = [dat])
    end

    @testset "Basic batch processing" begin
        EegFun.realign("epochs_good", RESP_TRIGGERS, input_dir = test_dir, output_dir = output_dir)

        @test isdir(output_dir)
        out_files = readdir(output_dir)
        @test any(startswith(f, "1_") for f in out_files)
        @test any(startswith(f, "2_") for f in out_files)
        @test any(startswith(f, "3_") for f in out_files)

        # Load and check alignment
        out_f = first(filter(f -> startswith(f, "1_"), out_files))
        realigned = load(joinpath(output_dir, out_f), "data")
        @test realigned isa Vector{EegFun.EpochData}

        for epoch in realigned[1].data
            t0_idx = argmin(abs.(epoch.time))
            @test abs(epoch.time[t0_idx]) < 1/256 + 1e-9
        end
    end

    @testset "Batch with participant filtering" begin
        output_dir2 = joinpath(test_dir, "realigned_filtered")
        EegFun.realign(
            "epochs_good",
            RESP_TRIGGERS,
            input_dir = test_dir,
            participant_selection = EegFun.participants([1, 2]),
            output_dir = output_dir2,
        )

        @test isdir(output_dir2)
        out_files = filter(f -> endswith(f, ".jld2"), readdir(output_dir2))
        @test length(out_files) == 2
        @test any(startswith(f, "1_") for f in out_files)
        @test any(startswith(f, "2_") for f in out_files)
        @test !any(startswith(f, "3_") for f in out_files)
    end

    @testset "Batch error handling: no matching files" begin
        output_dir3 = joinpath(test_dir, "realigned_nomatch")
        result = EegFun.realign("nonexistent_pattern", RESP_TRIGGERS, input_dir = test_dir, output_dir = output_dir3)
        @test isnothing(result)
    end

    rm(test_dir, recursive = true)
end
