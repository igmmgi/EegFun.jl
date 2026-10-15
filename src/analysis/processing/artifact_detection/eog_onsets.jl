"""
    detect_eog_onsets!(dat::ContinuousData, criterion::Float64, channel_in::Symbol, channel_out::Symbol)

Detects EOG (electrooculogram) onsets in the EEG data based on a specified criterion.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `criterion::Real`: The threshold for detecting EOG onsets.
- `channel_in::Symbol`: The channel from which to detect EOG onsets.
- `channel_out::Symbol`: The channel where the detected EOG onsets will be recorded as boolean values, indicating the presence of an EOG event.

# Returns
Nothing. The function modifies the input data in place.

# Examples
```julia
detect_eog_onsets!(dat, 50.0, :vEOG, :is_vEOG) # Detect vertical EOG onsets
detect_eog_onsets!(dat, 30.0, :hEOG, :is_hEOG) # Detect horizontal EOG onsets
```
"""

function detect_eog_onsets!(dat::ContinuousData, criterion::Real, channel_in::Symbol, channel_out::Symbol; step_size::Int = 20)
    @info "Detecting EOG onsets in channel $(channel_in) with stepsize criterion $(criterion) μV"
    channel_in ∉ propertynames(dat.data) && @minimal_error("channel $(channel_in) not found in data")

    step_samples = div(dat.sample_rate, step_size)
    eog_diff = diff(dat.data[1:step_samples:end, channel_in])
    eog_idx = findall(x -> abs(x) >= criterion, eog_diff)
    eog_idx = [idx for (i, idx) in enumerate(eog_idx) if i == 1 || (idx - eog_idx[i-1] > 2)] .* step_samples
    filter!(idx -> idx <= n_samples(dat), eog_idx)
    dat.data[!, channel_out] = falses(n_samples(dat))
    dat.data[eog_idx, channel_out] .= true
    return nothing
end


"""
    detect_eog_onsets!(dat::EegData, eog_cfg::Dict)

Detect EOG onsets for both vertical and horizontal EOG channels based on configuration.

# Arguments
- `dat::EegData`: The EEG data object
- `eog_cfg::Dict`: EOG configuration dictionary containing vEOG and hEOG settings

# Example
```julia
eog_cfg = Dict(
    "vEOG_criterion" => 50.0,
    "hEOG_criterion" => 50.0,
    "vEOG_channels" => [["Fp1"], ["IO1"], ["vEOG"]],
    "hEOG_channels" => [["F9"], ["F10"], ["hEOG"]]
)
detect_eog_onsets!(dat, eog_cfg)
```
"""
function detect_eog_onsets!(dat::EegData, eog_cfg::Dict)
    vEOG_cfg = eog_cfg["vEOG_channels"]
    detect_eog_onsets!(dat, eog_cfg["vEOG_criterion"], Symbol(vEOG_cfg[3][1]), Symbol("is_" * vEOG_cfg[3][1]))
    hEOG_cfg = eog_cfg["hEOG_channels"]
    detect_eog_onsets!(dat, eog_cfg["hEOG_criterion"], Symbol(hEOG_cfg[3][1]), Symbol("is_" * hEOG_cfg[3][1]))
end


