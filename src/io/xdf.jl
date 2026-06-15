# =============================================================================
# XDF IMPORT
# =============================================================================

"""
    read_xdf(filepath::String; select_streams, sync, dejitter_timestamps, handle_clock_resets) -> ExtensibleDataFormat.XdfData

Read an Extensible Data Format (`.xdf` / `.xdfz` / `.xdf.gz`) file.

XDF is a universal container format used by Lab Streaming Layer (LSL). A single file
may contain multiple streams (e.g., EEG, Markers, physiological signals). EegFun
automatically identifies the `\"EEG\"` stream and any `\"Markers\"` stream when
converting to `ContinuousData` via [`create_eegfun_data`](@ref).

# Arguments
- `filepath::String`: Path to the `.xdf`, `.xdfz`, or `.xdf.gz` file.

# Keyword Arguments
- `select_streams::Union{Nothing,Vector{Int}}=nothing`: Only load the given stream IDs.
  Pass `nothing` (default) to load all streams.
- `sync::Bool=true`: Enable clock synchronisation across streams.
- `dejitter_timestamps::Bool=true`: Remove timestamp jitter from regularly-sampled streams.
- `handle_clock_resets::Bool=true`: Detect and correct computer clock resets within the
  recording.

# Returns
- `ExtensibleDataFormat.XdfData`: Raw XDF data containing all parsed streams. Pass the
  result to [`create_eegfun_data`](@ref) to convert to a `ContinuousData` object.

# Examples
```julia
using EegFun

# Read all streams
raw = EegFun.read_xdf("recording.xdf")
dat = EegFun.create_eegfun_data(raw)

# Or use read_raw_data for automatic format detection
dat = EegFun.create_eegfun_data(EegFun.read_raw_data("recording.xdf"))

# Read only specific stream IDs (e.g., stream 1)
raw = EegFun.read_xdf("recording.xdf"; select_streams = [1])
```

# See Also
- [`read_raw_data`](@ref): Format-agnostic reader that dispatches to this function.
- [`create_eegfun_data`](@ref): Converts `XdfData` to `ContinuousData`.
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
