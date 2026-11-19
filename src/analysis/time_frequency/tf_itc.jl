# Inter-trial phase coherence (ITC) helpers.
# ITC from the decomposition kernels is requested with `tf_morlet(...; output = :itc)` / `tf_stft(...; output = :itc)`.

# ITC results are ordinary TimeFreqData tagged via the `method` field (e.g. :wavelet_itc)
_is_itc(tf::Union{TimeFreqData,TimeFreqEpochData}) = endswith(String(tf.method), "_itc")

# Label for the quantity stored in `data_power` (used by colorbars and plot titles)
_tf_quantity_label(tf::Union{TimeFreqData,TimeFreqEpochData}) = _is_itc(tf) ? "ITC" : "Power"

"""
    tf_itc(tf::TimeFreqEpochData) -> TimeFreqData
    tf_itc(tfs::AbstractVector{<:TimeFreqEpochData}) -> Vector{TimeFreqData}

Compute inter-trial phase coherence (ITC) from single-trial phases, i.e. ``|\\frac{1}{N}\\sum_k e^{i\\varphi_k}|``
at every time-frequency point. Requires data decomposed with `return_trials = true, return_phase = true`.

Returns a `TimeFreqData` whose `data_power` holds ITC values (0–1), `data_phase` the preferred (mean) phase,
and whose `method` is tagged with an `_itc` suffix. This gives the same result as
`tf_morlet(...; output = :itc)` / `tf_stft(...; output = :itc)`, which avoid storing every trial and are
preferred when single-trial data are not otherwise needed.

ITC is biased upwards by small trial counts (chance level ≈ `sqrt(π / (4N))`), so compare conditions with
equal trial numbers.

# Examples
```julia
tf_trials = tf_morlet(epochs; frequencies = logrange(2, 40, length = 20), cycles = (3, 10), return_trials = true, return_phase = true)
itc = tf_itc(tf_trials)
```
"""
function tf_itc(tf::TimeFreqEpochData)::TimeFreqData
    n_trials = length(tf.data_phase)
    n_trials < 2 && @minimal_error("ITC requires at least 2 trials, got $n_trials")
    _is_itc(tf) && @minimal_error("Data are already ITC values")
    if n_trials < 10
        @minimal_warning "Small trial count ($n_trials < 10): ITC is positively biased at low trial counts (chance ≈ $(round(sqrt(pi / (4 * n_trials)); digits = 3)))."
    end

    ref = tf.data_phase[1]
    n_rows = nrow(ref)
    all(df -> nrow(df) == n_rows, tf.data_phase) || @minimal_error("All trials must share the same time-frequency grid")
    has_power = length(tf.data_power) == n_trials && all(df -> nrow(df) == n_rows, tf.data_power)

    itc_df = DataFrame(time = copy(ref.time), freq = copy(ref.freq))
    phase_df = DataFrame(time = copy(ref.time), freq = copy(ref.freq))
    phasor_sum = Vector{ComplexF64}(undef, n_rows)
    inv_n = 1.0 / n_trials
    raw_unbaselined = isnothing(tf.baseline)

    for ch in channel_labels(tf)
        all(isnan, ref[!, ch]) &&
            @minimal_error("tf_itc requires single-trial phase: rerun the decomposition with `return_trials = true, return_phase = true`")
        fill!(phasor_sum, 0.0im)
        for tr = 1:n_trials
            pcol = tf.data_phase[tr][!, ch]
            if has_power && hasproperty(tf.data_power[tr], ch)
                pwcol = tf.data_power[tr][!, ch]
                @inbounds for i in eachindex(phasor_sum, pcol, pwcol)
                    ph = pcol[i]
                    pw = pwcol[i]
                    if isnan(ph) || isnan(pw)
                        phasor_sum[i] = NaN + NaN * im
                    elseif isfinite(ph) && isfinite(pw) && (!raw_unbaselined || pw > 0.0)
                        phasor_sum[i] += cis(ph)
                    end
                end
            else
                @inbounds for i in eachindex(phasor_sum, pcol)
                    ph = pcol[i]
                    phasor_sum[i] += isnan(ph) ? (NaN + NaN * im) : cis(ph)
                end
            end
        end
        itc_df[!, ch] = abs.(phasor_sum) .* inv_n
        phase_df[!, ch] = angle.(phasor_sum)
    end

    return TimeFreqData(
        tf.file,
        tf.condition,
        tf.condition_name,
        itc_df,
        phase_df,
        copy(tf.layout),
        tf.sample_rate,
        Symbol(tf.method, :_itc),
        nothing,
        copy(tf.analysis_info),
    )
end

tf_itc(tfs::AbstractVector{<:TimeFreqEpochData}) = [tf_itc(tf) for tf in tfs]
