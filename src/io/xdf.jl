"""
    read_xdf(filepath::String; select_streams = nothing, sync = true, dejitter_timestamps = true, handle_clock_resets = true) -> ExtensibleDataFormat.XdfData

Read an Extensible Data Format (`.xdf`, `.xdfz`, or `.xdf.gz`) Lab Streaming Layer container file.
Pass the returned `XdfData` to [`create_eegfun_data`](@ref) to convert it into `ContinuousData`.

# Keyword Arguments
- `select_streams::Union{Nothing,Vector{Int}} = nothing`: Stream IDs to load (`nothing` loads all streams)
- `sync::Bool = true`: Enable cross-stream clock synchronization
- `dejitter_timestamps::Bool = true`: Remove timestamp jitter from regularly sampled streams
- `handle_clock_resets::Bool = true`: Detect and correct clock resets within the recording
"""
function read_xdf(
    filepath::String;
    select_streams::Union{Nothing,Vector{Int}} = nothing,
    sync::Bool = true,
    dejitter_timestamps::Bool = true,
    handle_clock_resets::Bool = true,
)
    @info "Reading XDF data: $(basename(filepath))"
    return ExtensibleDataFormat.read_xdf(
        filepath;
        select_streams = select_streams,
        sync = sync,
        dejitter_timestamps = dejitter_timestamps,
        handle_clock_resets = handle_clock_resets,
    )
end
