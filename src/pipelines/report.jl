# Self-contained HTML quality-control report generator for preprocessing pipelines.

######################################
######### Types & Constants ##########
######################################

const _BASE64_CHARS = UInt8.(collect("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"))
const _FRONTAL_ROI_CANDIDATES = (:Fp1, :Fp2, :Fpz, :AF7, :AF8, :AF3, :AF4, :F1, :F2, :Fz)
const _POSTERIOR_ROI_CANDIDATES = (:O1, :O2, :Oz, :Iz, :PO7, :PO8, :PO3, :PO4, :POz, :P1, :P2, :Pz)

struct _PsdSummary
    freqs::Vector{Float64}
    mean_db::Vector{Float64}
    posterior_db::Union{Vector{Float64},Nothing}
    frontal_db::Union{Vector{Float64},Nothing}
    is_epoch_level::Bool
end

struct _ParticipantReportData
    stem::String
    artifact_info::Union{ArtifactInfo,Nothing}
    ica::Union{InfoIca,Nothing}
    erps::Union{Vector{ErpData},Nothing}
    erps_raw::Union{Vector{ErpData},Nothing}
    psd::Union{_PsdSummary,Nothing}
    repaired_channels::Vector{Symbol}
    skipped_channels::Vector{Symbol}
    repair_method::String
    epoch_channel_repairs::Dict{Symbol,Int}
    ic_counts::NamedTuple{(:total, :vEOG, :hEOG, :ECG, :line_noise, :channel_noise, :manual),NTuple{7,Int}}
    condition_df::DataFrame
    total_original::Int
    total_kept::Int
    overall_kept_pct::Float64
    is_flagged::Bool
    flag_reasons::Vector{String}
end

const _REPORT_CSS = """
:root {
    --bg: #f8fafc;
    --card: #ffffff;
    --border: #e2e8f0;
    --text: #0f172a;
    --muted: #475569;
    --accent: #2563eb;
    --ok-bg: #dcfce7;
    --ok-fg: #166534;
    --warn-bg: #fef3c7;
    --warn-fg: #92400e;
    --bad-bg: #fee2e2;
    --bad-fg: #991b1b;
}
* { box-sizing: border-box; }
body {
    margin: 0;
    padding: 24px;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    background: var(--bg);
    color: var(--text);
    line-height: 1.5;
}
.container { max-width: 1320px; margin: 0 auto; }
header.hero {
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: 12px;
    padding: 20px 24px;
    margin-bottom: 20px;
}
header.hero h1 { margin: 0 0 6px 0; font-size: 1.6rem; }
header.hero .meta { color: var(--muted); font-size: 0.9rem; }
.kpi-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
    gap: 14px;
    margin-bottom: 20px;
}
.kpi {
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 16px 18px;
}
.kpi .label { font-size: 0.82rem; text-transform: uppercase; letter-spacing: 0.04em; color: var(--muted); }
.kpi .value { font-size: 1.65rem; font-weight: 700; margin-top: 4px; }
.kpi .sub { font-size: 0.82rem; color: var(--muted); margin-top: 2px; }
section.panel {
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: 12px;
    padding: 20px 24px;
    margin-bottom: 20px;
}
section.panel h2 { margin: 0 0 14px 0; font-size: 1.2rem; border-bottom: 1px solid var(--border); padding-bottom: 8px; }
section.panel h3 { margin: 16px 0 10px 0; font-size: 1.02rem; }
.section-head {
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    justify-content: space-between;
    gap: 12px;
    border-bottom: 1px solid var(--border);
    padding-bottom: 8px;
    margin-bottom: 14px;
}
.section-head h2 { margin: 0; border-bottom: none; padding-bottom: 0; }
.btn-group { display: flex; flex-wrap: wrap; gap: 8px; align-items: center; }
button.btn-ctrl, a.btn-ctrl {
    appearance: none;
    display: inline-flex;
    align-items: center;
    background: #f8fafc;
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 5px 12px;
    font-size: 0.82rem;
    font-weight: 600;
    color: var(--text);
    text-decoration: none;
    cursor: pointer;
}
button.btn-ctrl:hover, a.btn-ctrl:hover { background: #e2e8f0; }
button.btn-ctrl.active { background: var(--warn-bg); border-color: #f59e0b; color: var(--warn-fg); }
.methods-box {
    margin-top: 14px;
    padding: 14px 16px;
    background: #f8fafc;
    border: 1px solid var(--border);
    border-left: 4px solid var(--accent);
    border-radius: 8px;
    font-size: 0.9rem;
}
.methods-head {
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 10px;
    margin-bottom: 6px;
    font-weight: 600;
    font-size: 0.84rem;
    color: var(--muted);
    text-transform: uppercase;
    letter-spacing: 0.03em;
}
.methods-box p { margin: 0; color: var(--text); line-height: 1.6; }
.cfg-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
    gap: 8px 20px;
    font-size: 0.9rem;
}
.cfg-item span.k { color: var(--muted); font-weight: 600; }
table {
    width: 100%;
    border-collapse: collapse;
    font-size: 0.88rem;
}
th, td {
    padding: 8px 10px;
    border-bottom: 1px solid var(--border);
    text-align: left;
    vertical-align: middle;
}
th { background: #f1f5f9; font-weight: 600; color: #334155; }
td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
tr.flagged-row { background: #fffbeb; }
body.flagged-only tr.ok-row { display: none; }
body.flagged-only details.subj-card:not(.flagged) { display: none; }
a.subj-link { color: var(--accent); text-decoration: none; font-weight: 600; }
a.subj-link:hover { text-decoration: underline; }
a.top-link {
    color: var(--muted);
    text-decoration: none;
    font-size: 0.76rem;
    font-weight: 600;
    padding: 1px 7px;
    border: 1px solid var(--border);
    border-radius: 4px;
    background: #ffffff;
}
a.top-link:hover { color: var(--accent); border-color: var(--accent); }
.badge {
    display: inline-block;
    padding: 2px 8px;
    border-radius: 999px;
    font-size: 0.76rem;
    font-weight: 600;
}
.badge-ok { background: var(--ok-bg); color: var(--ok-fg); }
.badge-warn { background: var(--warn-bg); color: var(--warn-fg); }
.badge-bad { background: var(--bad-bg); color: var(--bad-fg); }
.fig-stack {
    display: flex;
    flex-direction: column;
    gap: 14px;
    margin-top: 14px;
}
.fig-row-top {
    display: grid;
    grid-template-columns: minmax(260px, 1fr) minmax(420px, 2fr);
    gap: 16px;
    align-items: stretch;
}
.fig-row-single {
    display: grid;
    grid-template-columns: 1fr;
    gap: 16px;
}
@media (max-width: 860px) {
    .fig-row-top { grid-template-columns: 1fr; }
}
.fig-box {
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 10px;
    background: #ffffff;
    text-align: center;
    display: flex;
    flex-direction: column;
    justify-content: space-between;
}
.fig-box img { max-width: 100%; height: auto; display: block; margin: auto; }
.fig-box .cap { font-size: 0.82rem; color: var(--muted); margin-top: 6px; }
details.subj-card {
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: 10px;
    margin-bottom: 12px;
    overflow: hidden;
}
details.subj-card > summary {
    padding: 12px 18px;
    cursor: pointer;
    font-weight: 600;
    display: flex;
    flex-wrap: wrap;
    gap: 12px;
    align-items: center;
    justify-content: space-between;
    background: #f8fafc;
}
details.subj-card.flagged > summary { background: #fffbeb; }
details.subj-card .card-body { padding: 16px 20px; border-top: 1px solid var(--border); }
.summary-pills { display: flex; flex-wrap: wrap; gap: 8px 14px; font-weight: 400; font-size: 0.84rem; color: var(--muted); align-items: center; }
@media print {
    body { padding: 0; background: #fff; }
    details.subj-card { break-inside: avoid; }
    .btn-group { display: none; }
}
"""

######################################
########### User Functions ###########
######################################

"""
    pipeline_report(output_dir = "./preprocessed_files"; output_file = "pipeline_report.html", save_csv = true, title = "EegFun Preprocessing QC Report", low_epoch_threshold = 70.0, max_ica_pct = 15.0, include_plots = true) -> String

Generate a self-contained HTML quality-control report (and optional `pipeline_qc_summary.csv`) from `EegFun` pipeline outputs in `output_dir`.

Reads available `config.toml`, `neighbors.toml`, `epoch_summary.jld2`, `*_artifact_info.jld2`, `*_ica.jld2`,
`*_erps.jld2`, `*_erps_raw.jld2`, `*_epochs.jld2`, and `*_epochs_raw.jld2` files to summarize configuration,
methods/attrition text, participant flags, sensor repairs, grand-average ERPs/GFP and topographies, Welch PSDs,
and pre- vs. post-cleaning overlays.

# Arguments
- `output_dir::AbstractString = "./preprocessed_files"`: Directory containing saved preprocessing pipeline outputs.

# Keyword Arguments
- `output_file::AbstractString = "pipeline_report.html"`: Output HTML filename (resolved inside `output_dir` when relative, or used directly when absolute).
- `save_csv::Bool = true`: Also write `pipeline_qc_summary.csv` alongside `output_file` for direct tabular analysis in Julia, R, Python, or spreadsheets.
- `title::AbstractString = "EegFun Preprocessing QC Report"`: Report heading and HTML page `<title>`.
- `low_epoch_threshold::Real = 70.0`: Retention percentage threshold below which participants/conditions are flagged for review.
- `max_ica_pct::Real = 15.0`: Maximum percentage of ICA components removed (for participants with ≥ 3 removed components) before flagging for review.
- `include_plots::Bool = true`: Embed Base64 PNG figures rendered via Makie (requires a loaded backend such as `CairoMakie`; falls back to tabular-only output with a warning if no backend is active).

# Returns
- `String`: Absolute path to the generated `.html` report file.

# Examples
```julia
using EegFun, CairoMakie

EegFun.preprocess("pipeline.toml")
report_path = EegFun.pipeline_report("./preprocessed_files")
```

# See Also
- [`preprocess`](@ref), [`summarize_electrode_repairs`](@ref), [`summarize_ica_components`](@ref), [`compare_rejections`](@ref)
"""
function pipeline_report(
    output_dir::AbstractString = "./preprocessed_files";
    output_file::AbstractString = "pipeline_report.html",
    save_csv::Bool = true,
    title::AbstractString = "EegFun Preprocessing QC Report",
    low_epoch_threshold::Real = 70.0,
    max_ica_pct::Real = 15.0,
    include_plots::Bool = true,
)::String
    dir_str = String(output_dir)
    isdir(dir_str) || @minimal_error("Output directory does not exist: $dir_str")
    0.0 <= low_epoch_threshold <= 100.0 || @minimal_error("low_epoch_threshold must be between 0 and 100, got $low_epoch_threshold")
    0.0 <= max_ica_pct <= 100.0 || @minimal_error("max_ica_pct must be between 0 and 100, got $max_ica_pct")

    stems = _discover_pipeline_stems(dir_str)
    isempty(stems) && @minimal_error("No preprocessed EegFun `.jld2` files found in: $dir_str")

    can_plot = include_plots && _has_renderable_makie_backend()
    if include_plots && !can_plot
        @minimal_warning "No active Makie backend detected (run `using CairoMakie` to embed plots); generating tabular HTML report."
    end

    cfg_path = joinpath(dir_str, "config.toml")
    pipeline_cfg = isfile(cfg_path) ? _safe_load_pipeline_cfg(cfg_path) : nothing
    epoch_summary_df = _safe_load_jld2_df(joinpath(dir_str, "epoch_summary.jld2"))

    records = _load_participant_records(dir_str, stems, epoch_summary_df, Float64(low_epoch_threshold), Float64(max_ica_pct), can_plot)
    _, electrode_freq_df = summarize_electrode_repairs("_artifact_info"; input_dir = dir_str)

    report_path = isabspath(output_file) ? abspath(String(output_file)) : abspath(joinpath(dir_str, String(output_file)))
    mkpath(dirname(report_path))

    csv_text = _build_qc_summary_csv(records)
    if save_csv
        csv_path = joinpath(dirname(report_path), "pipeline_qc_summary.csv")
        write(csv_path, csv_text)
    end

    html = _build_pipeline_report_html(
        dir_str,
        String(title),
        Float64(low_epoch_threshold),
        Float64(max_ica_pct),
        pipeline_cfg,
        records,
        electrode_freq_df,
        can_plot,
        csv_text,
    )
    write(report_path, html)
    @info "Generated preprocessing QC report: $report_path"
    return report_path
end

######################################
######### Internal Functions #########
######################################

function _base64_encode(bytes::AbstractVector{UInt8})::String
    n = length(bytes)
    n == 0 && return ""
    out_len = 4 * cld(n, 3)
    out = Vector{UInt8}(undef, out_len)
    i = 1
    o = 1
    @inbounds while i + 2 <= n
        b1, b2, b3 = bytes[i], bytes[i+1], bytes[i+2]
        out[o] = _BASE64_CHARS[(b1>>2)+1]
        out[o+1] = _BASE64_CHARS[(((b1&0x03)<<4)|(b2>>4))+1]
        out[o+2] = _BASE64_CHARS[(((b2&0x0f)<<2)|(b3>>6))+1]
        out[o+3] = _BASE64_CHARS[(b3&0x3f)+1]
        i += 3
        o += 4
    end
    rem_bytes = n - i + 1
    @inbounds if rem_bytes == 1
        b1 = bytes[i]
        out[o] = _BASE64_CHARS[(b1>>2)+1]
        out[o+1] = _BASE64_CHARS[((b1&0x03)<<4)+1]
        out[o+2] = UInt8('=')
        out[o+3] = UInt8('=')
    elseif rem_bytes == 2
        b1, b2 = bytes[i], bytes[i+1]
        out[o] = _BASE64_CHARS[(b1>>2)+1]
        out[o+1] = _BASE64_CHARS[(((b1&0x03)<<4)|(b2>>4))+1]
        out[o+2] = _BASE64_CHARS[((b2&0x0f)<<2)+1]
        out[o+3] = UInt8('=')
    end
    return String(out)
end

function _html_escape(s::AbstractString)::String
    io = IOBuffer(; sizehint = sizeof(s) + 16)
    for c in s
        if c == '&'
            write(io, "&amp;")
        elseif c == '<'
            write(io, "&lt;")
        elseif c == '>'
            write(io, "&gt;")
        elseif c == '"'
            write(io, "&quot;")
        elseif c == '\''
            write(io, "&#39;")
        else
            write(io, c)
        end
    end
    return String(take!(io))
end

function _csv_escape(s::AbstractString)::String
    if occursin(r"[\",\r\n]", s)
        return "\"" * replace(s, "\"" => "\"\"") * "\""
    end
    return String(s)
end

_slug_id(s::AbstractString) = "subj-" * replace(s, r"[^A-Za-z0-9_-]" => "_")

function _unique_slugs(records::Vector{_ParticipantReportData})::Vector{String}
    seen = Dict{String,Int}()
    out = Vector{String}(undef, length(records))
    for (i, r) in enumerate(records)
        base = _slug_id(r.stem)
        cnt = get(seen, base, 0) + 1
        seen[base] = cnt
        out[i] = cnt == 1 ? base : "$base-$cnt"
    end
    return out
end

function _has_renderable_makie_backend()::Bool
    backend = try
        Makie.current_backend()
    catch
        missing
    end
    return !ismissing(backend) && !isnothing(backend)
end

function _fig_to_data_uri(fig::Figure; px_per_unit::Float64 = 1.0)::Union{String,Nothing}
    try
        saved_bytes = mktempdir() do tmp
            png_path = joinpath(tmp, "fig.png")
            Makie.save(png_path, fig; px_per_unit = px_per_unit)
            read(png_path)
        end
        !isempty(saved_bytes) && return "data:image/png;base64," * _base64_encode(saved_bytes)
    catch
    end
    try
        buf = IOBuffer()
        if showable(MIME"image/png"(), fig)
            show(buf, MIME"image/png"(), fig)
            bytes = take!(buf)
            !isempty(bytes) && return "data:image/png;base64," * _base64_encode(bytes)
        end
        return nothing
    catch e
        @debug "Failed to encode figure for pipeline report" exception = e
        return nothing
    end
end

function _discover_pipeline_stems(output_dir::String)::Vector{String}
    # Order longest suffixes first so `_erps_raw.jld2` never creates a `foo_erps_raw` stem
    suffixes = (
        "_artifact_info.jld2",
        "_continuous_corrected.jld2",
        "_continuous_raw.jld2",
        "_epochs_corrected.jld2",
        "_epochs_raw.jld2",
        "_erps_corrected.jld2",
        "_erps_raw.jld2",
        "_epochs.jld2",
        "_erps.jld2",
        "_ica.jld2",
    )
    raw_only_suffixes = ("_continuous_corrected.jld2", "_continuous_raw.jld2", "_epochs_raw.jld2", "_erps_raw.jld2")
    stems = Set{String}()
    for f in readdir(output_dir)
        f in ("epoch_summary.jld2", "file_summary.jld2") && continue
        any(s -> endswith(f, s), raw_only_suffixes) && continue
        for suf in suffixes
            if endswith(f, suf)
                stem = f[1:(end-ncodeunits(suf))]
                !isempty(stem) && push!(stems, stem)
                break
            end
        end
    end
    return sort!(collect(stems); by = _natural_sort_key)
end

function _safe_load_pipeline_cfg(cfg_path::String)::Union{PipelineConfig,Nothing}
    try
        return PipelineConfig(TOML.parsefile(cfg_path))
    catch e
        @minimal_warning "Could not parse $cfg_path for report: $(sprint(showerror, e))"
        return nothing
    end
end

function _safe_load_jld2_df(path::String)::Union{DataFrame,Nothing}
    isfile(path) || return nothing
    try
        val = JLD2.load(path, "data")
        return val isa DataFrame ? val : nothing
    catch
        return nothing
    end
end

function _safe_read_jld2_typed(path::String, ::Type{T}) where {T}
    isfile(path) || return nothing
    try
        val = read_data(path)
        return val isa T ? val : nothing
    catch
        return nothing
    end
end

function _coerce_erp_vec(val)::Union{Vector{ErpData},Nothing}
    if val isa ErpData
        return ErpData[val]
    elseif val isa AbstractVector && !isempty(val) && all(x -> x isa ErpData, val)
        return ErpData[x for x in val]
    end
    return nothing
end

function _coerce_epoch_vec(val)::Union{Vector{EpochData},Nothing}
    if val isa EpochData
        return EpochData[val]
    elseif val isa AbstractVector && !isempty(val) && all(x -> x isa EpochData, val)
        return EpochData[x for x in val]
    end
    return nothing
end

function _compute_psd_from_dfs(dfs::Vector{DataFrame}, fs::Float64, chs::Vector{Symbol}, is_epoch_level::Bool)::Union{_PsdSummary,Nothing}
    (isempty(dfs) || isempty(chs) || fs <= 0.0) && return nothing
    n_s = nrow(dfs[1])
    n_s < 16 && return nothing

    freqs_full = collect(FFTW.rfftfreq(n_s, fs))
    max_f = min(60.0, fs / 2.0 - 0.5)
    f_idx = findall(f -> 1.0 <= f <= max_f, freqs_full)
    length(f_idx) < 3 && return nothing

    step = max(1, fld(length(dfs), 80))
    chosen_indices = 1:step:length(dfs)

    w = DSP.hanning(n_s)
    win_norm = fs * sum(abs2, w)
    n_chs = length(chs)
    buf = Matrix{Float64}(undef, n_s, n_chs)
    power_acc = zeros(Float64, length(f_idx), n_chs)
    n_used = 0

    for idx in chosen_indices
        df = dfs[idx]
        nrow(df) == n_s || continue
        valid_df = true
        @inbounds for (c_i, ch) in enumerate(chs)
            if !hasproperty(df, ch)
                valid_df = false
                break
            end
            col = df[!, ch]
            if !(eltype(col) <: Real)
                valid_df = false
                break
            end
            sum_fin = 0.0
            cnt_fin = 0
            for s = 1:n_s
                val = Float64(col[s])
                if isfinite(val)
                    sum_fin += val
                    cnt_fin += 1
                end
            end
            m = cnt_fin > 0 ? (sum_fin / cnt_fin) : 0.0
            for s = 1:n_s
                val = Float64(col[s])
                buf[s, c_i] = isfinite(val) ? w[s] * (val - m) : 0.0
            end
        end
        valid_df || continue
        spec = FFTW.rfft(buf, 1)
        @inbounds for c_i = 1:n_chs
            for (k_out, k_in) in enumerate(f_idx)
                power_acc[k_out, c_i] += abs2(spec[k_in, c_i])
            end
        end
        n_used += 1
    end
    n_used == 0 && return nothing

    inv_scale = 1.0 / (n_used * win_norm)
    nyq = fs / 2.0
    @inbounds for c_i = 1:n_chs
        for (k_out, k_in) in enumerate(f_idx)
            fk = freqs_full[k_in]
            onesided = (fk > 0.0 && fk < nyq) ? 2.0 : 1.0
            power_acc[k_out, c_i] *= onesided * inv_scale
        end
    end

    mean_p = vec(mean(power_acc; dims = 2))
    mean_db = [10.0 * log10(max(v, 1e-12)) for v in mean_p]

    post_idx = findall(in(_POSTERIOR_ROI_CANDIDATES), chs)
    front_idx = findall(in(_FRONTAL_ROI_CANDIDATES), chs)
    post_db = isempty(post_idx) ? nothing : [10.0 * log10(max(v, 1e-12)) for v in vec(mean(@view(power_acc[:, post_idx]); dims = 2))]
    front_db = isempty(front_idx) ? nothing : [10.0 * log10(max(v, 1e-12)) for v in vec(mean(@view(power_acc[:, front_idx]); dims = 2))]

    return _PsdSummary(freqs_full[f_idx], mean_db, post_db, front_db, is_epoch_level)
end

function _compute_psd_summary(epochs::Vector{EpochData})::Union{_PsdSummary,Nothing}
    isempty(epochs) && return nothing
    chs = get_selected_channels(epochs[1], channels(); include_meta = false, include_extra = false)
    all_dfs = DataFrame[]
    for ep in epochs
        append!(all_dfs, ep.data)
    end
    return _compute_psd_from_dfs(all_dfs, Float64(epochs[1].sample_rate), chs, true)
end

function _compute_psd_summary(erps::Vector{ErpData})::Union{_PsdSummary,Nothing}
    isempty(erps) && return nothing
    chs = get_selected_channels(erps[1], channels(); include_meta = false, include_extra = false)
    dfs = [e.data for e in erps]
    return _compute_psd_from_dfs(dfs, Float64(erps[1].sample_rate), chs, false)
end

function _extract_ic_counts(ca::Union{ArtifactComponents,Nothing})
    isnothing(ca) && return (total = 0, vEOG = 0, hEOG = 0, ECG = 0, line_noise = 0, channel_noise = 0, manual = 0)
    veog = length(get(ca.eog, :vEOG, Int[]))
    heog = length(get(ca.eog, :hEOG, Int[]))
    ecg = length(ca.ecg)
    ln = length(ca.line_noise)
    cn = length(ca.channel_noise)
    man = length(ca.manual)
    tot = length(unique(get_all_ica_components(ca)))
    return (total = tot, vEOG = veog, hEOG = heog, ECG = ecg, line_noise = ln, channel_noise = cn, manual = man)
end

function _build_condition_rejection_df(
    stem::String,
    artifact_info::Union{ArtifactInfo,Nothing},
    erps::Union{Vector{ErpData},Nothing},
    epochs::Union{Vector{EpochData},Nothing},
    epoch_summary_df::Union{DataFrame,Nothing},
)::DataFrame
    if !isnothing(artifact_info) && !isempty(artifact_info.epoch_rejections)
        rejs = artifact_info.epoch_rejections
        step1 = filter(r -> r.name == "rejection_step1", rejs)
        step2 = filter(r -> r.name == "rejection_step2", rejs)
        if !isempty(step1) && length(step1) == length(step2)
            cmp = compare_rejections(step1, step2)
            n_kept = cmp.n_original .- cmp.rejected_after_repair
            return DataFrame(
                condition = cmp.condition,
                condition_name = cmp.condition_name,
                n_original = cmp.n_original,
                n_repaired = cmp.successfully_repaired,
                n_rejected = cmp.rejected_after_repair,
                n_kept = n_kept,
                pct_kept = cmp.pct_final_kept,
            )
        else
            conds = [r.info.number for r in rejs]
            cnames = [r.info.name for r in rejs]
            n_orig = [r.info.n for r in rejs]
            n_rej = [length(unique_epochs(r)) for r in rejs]
            n_rep = [isnothing(r.repaired) ? 0 : count(!isempty, values(r.repaired)) for r in rejs]
            n_kept = n_orig .- n_rej
            pct = [n_orig[i] > 0 ? round(100.0 * n_kept[i] / n_orig[i]; digits = 1) : 0.0 for i in eachindex(n_orig)]
            return DataFrame(
                condition = conds,
                condition_name = cnames,
                n_original = n_orig,
                n_repaired = n_rep,
                n_rejected = n_rej,
                n_kept = n_kept,
                pct_kept = pct,
            )
        end
    end

    if !isnothing(epoch_summary_df) && hasproperty(epoch_summary_df, :file)
        sub = epoch_summary_df[epoch_summary_df.file .== stem, :]
        if !isempty(sub)
            n_orig = hasproperty(sub, :original) ? sub.original : sub.final
            n_kept = hasproperty(sub, :final) ? sub.final : n_orig
            pct =
                hasproperty(sub, :percentage) ? Float64.(sub.percentage) :
                [o > 0 ? round(100.0 * k / o; digits = 1) : 0.0 for (o, k) in zip(n_orig, n_kept)]
            return DataFrame(
                condition = hasproperty(sub, :condition_num) ? sub.condition_num : collect(1:nrow(sub)),
                condition_name = hasproperty(sub, :condition_name) ? string.(sub.condition_name) : fill("Condition", nrow(sub)),
                n_original = n_orig,
                n_repaired = zeros(Int, nrow(sub)),
                n_rejected = n_orig .- n_kept,
                n_kept = n_kept,
                pct_kept = pct,
            )
        end
    end

    if !isnothing(erps) && !isempty(erps)
        return DataFrame(
            condition = [e.condition for e in erps],
            condition_name = [e.condition_name for e in erps],
            n_original = [e.n_epochs for e in erps],
            n_repaired = zeros(Int, length(erps)),
            n_rejected = zeros(Int, length(erps)),
            n_kept = [e.n_epochs for e in erps],
            pct_kept = fill(100.0, length(erps)),
        )
    elseif !isnothing(epochs) && !isempty(epochs)
        cnts = [length(ep.data) for ep in epochs]
        return DataFrame(
            condition = [ep.condition for ep in epochs],
            condition_name = [ep.condition_name for ep in epochs],
            n_original = cnts,
            n_repaired = zeros(Int, length(epochs)),
            n_rejected = zeros(Int, length(epochs)),
            n_kept = cnts,
            pct_kept = fill(100.0, length(epochs)),
        )
    end

    return DataFrame(
        condition = Int[],
        condition_name = String[],
        n_original = Int[],
        n_repaired = Int[],
        n_rejected = Int[],
        n_kept = Int[],
        pct_kept = Float64[],
    )
end

function _extract_epoch_channel_repairs(art::Union{ArtifactInfo,Nothing})::Dict{Symbol,Int}
    counts = Dict{Symbol,Int}()
    (isnothing(art) || isempty(art.epoch_rejections)) && return counts
    rejs = art.epoch_rejections
    step1 = filter(r -> r.name == "rejection_step1", rejs)
    step2 = filter(r -> r.name == "rejection_step2", rejs)
    target_rejs = !isempty(step1) ? step1 : rejs
    for (idx, r) in enumerate(target_rejs)
        isnothing(r.repaired) && continue
        rej_epochs = (!isempty(step1) && idx <= length(step2)) ? Set(unique_epochs(step2[idx])) : Set{Int}()
        for (ep_i, chs) in r.repaired
            ep_i in rej_epochs && continue
            for ch in chs
                counts[ch] = get(counts, ch, 0) + 1
            end
        end
    end
    return counts
end

function _load_participant_records(
    output_dir::String,
    stems::Vector{String},
    epoch_summary_df::Union{DataFrame,Nothing},
    low_epoch_threshold::Float64,
    max_ica_pct::Float64,
    can_plot::Bool,
)::Vector{_ParticipantReportData}
    records = Vector{_ParticipantReportData}(undef, length(stems))
    for (idx, stem) in enumerate(stems)
        art = _safe_read_jld2_typed(joinpath(output_dir, "$(stem)_artifact_info.jld2"), ArtifactInfo)
        ica = _safe_read_jld2_typed(joinpath(output_dir, "$(stem)_ica.jld2"), InfoIca)
        erps = _coerce_erp_vec(_safe_read_jld2_typed(joinpath(output_dir, "$(stem)_erps.jld2"), Any))

        # Stream epochs one participant at a time for PSD (and ERP fallback) without retaining epoch arrays in RAM
        ep_path = joinpath(output_dir, "$(stem)_epochs.jld2")
        epochs_tmp = (isnothing(erps) || can_plot) ? _coerce_epoch_vec(_safe_read_jld2_typed(ep_path, Any)) : nothing
        isnothing(erps) && !isnothing(epochs_tmp) && (erps = average_epochs(epochs_tmp))
        psd = if !can_plot
            nothing
        elseif !isnothing(epochs_tmp)
            _compute_psd_summary(epochs_tmp)
        elseif !isnothing(erps)
            _compute_psd_summary(erps)
        else
            nothing
        end

        erps_raw = if can_plot
            raw_e = _coerce_erp_vec(_safe_read_jld2_typed(joinpath(output_dir, "$(stem)_erps_raw.jld2"), Any))
            if isnothing(raw_e)
                raw_ep = _coerce_epoch_vec(_safe_read_jld2_typed(joinpath(output_dir, "$(stem)_epochs_raw.jld2"), Any))
                !isnothing(raw_ep) ? average_epochs(raw_ep) : nothing
            else
                raw_e
            end
        else
            nothing
        end

        rep_chs = Symbol[]
        skip_chs = Symbol[]
        rep_method = "—"
        if !isnothing(art) && !isempty(art.continuous_repairs)
            for cr in art.continuous_repairs
                append!(rep_chs, cr.repaired)
                append!(skip_chs, cr.skipped)
                rep_method = string(cr.method)
            end
            unique!(rep_chs)
            unique!(skip_chs)
        end

        ep_ch_repairs = _extract_epoch_channel_repairs(art)
        ic_counts = _extract_ic_counts(isnothing(art) ? nothing : art.ica_components)
        cond_df = _build_condition_rejection_df(stem, art, erps, epochs_tmp, epoch_summary_df)
        epochs_tmp = nothing
        tot_orig = isempty(cond_df) ? 0 : sum(cond_df.n_original)
        tot_kept = isempty(cond_df) ? 0 : sum(cond_df.n_kept)
        overall_pct = isempty(cond_df) ? NaN : round(mean(cond_df.pct_kept); digits = 1)

        reasons = String[]
        !isempty(skip_chs) && push!(reasons, "unrepairable bad channel(s): $(join(skip_chs, ", "))")
        if !isempty(cond_df)
            any(==(0), cond_df.n_kept) && push!(reasons, "0 clean epochs in ≥1 condition")
            low_conds = cond_df.condition_name[cond_df.pct_kept .< low_epoch_threshold]
            !isempty(low_conds) && push!(reasons, "retention < $(round(low_epoch_threshold; digits = 1))% ($(join(low_conds, ", ")))")
        end

        n_ica_ch = !isnothing(ica) ? size(ica.unmixing, 1) : (!isnothing(erps) && !isempty(erps) ? nrow(erps[1].layout.data) : 0)
        if n_ica_ch > 0 && ic_counts.total >= 3
            ica_pct = round(100.0 * ic_counts.total / n_ica_ch; digits = 1)
            if ica_pct > max_ica_pct
                push!(
                    reasons,
                    "high ICA removal: $(ic_counts.total)/$n_ica_ch components ($(ica_pct)% > $(round(max_ica_pct; digits = 1))%)",
                )
            end
        end

        records[idx] = _ParticipantReportData(
            stem,
            art,
            ica,
            erps,
            erps_raw,
            psd,
            rep_chs,
            skip_chs,
            rep_method,
            ep_ch_repairs,
            ic_counts,
            cond_df,
            tot_orig,
            tot_kept,
            overall_pct,
            !isempty(reasons),
            reasons,
        )
    end

    # Cohort-level check for truncated recordings (<80% of cohort median original trials)
    pos_totals = [r.total_original for r in records if r.total_original > 0]
    if length(pos_totals) >= 3
        med_orig = median(pos_totals)
        for i in eachindex(records)
            r = records[i]
            if 0 < r.total_original < 0.8 * med_orig
                new_reasons = vcat(
                    r.flag_reasons,
                    ["low trial count: $(r.total_original) original epochs (< 80% of cohort median $(round(Int, med_orig)))"],
                )
                records[i] = _ParticipantReportData(
                    r.stem,
                    r.artifact_info,
                    r.ica,
                    r.erps,
                    r.erps_raw,
                    r.psd,
                    r.repaired_channels,
                    r.skipped_channels,
                    r.repair_method,
                    r.epoch_channel_repairs,
                    r.ic_counts,
                    r.condition_df,
                    r.total_original,
                    r.total_kept,
                    r.overall_kept_pct,
                    true,
                    new_reasons,
                )
            end
        end
    end

    return records
end

function _participant_layout(rec::_ParticipantReportData)::Union{Layout,Nothing}
    !isnothing(rec.erps) && !isempty(rec.erps) && return rec.erps[1].layout
    !isnothing(rec.ica) && return rec.ica.layout
    return nothing
end

function _cohort_layout(records::Vector{_ParticipantReportData})::Union{Layout,Nothing}
    for rec in records
        lay = _participant_layout(rec)
        _has_spatial_layout(lay) && return lay
    end
    return nothing
end

function _has_spatial_layout(lay::Union{Layout,Nothing})::Bool
    (isnothing(lay) || isempty(lay.data)) && return false
    return has_valid_coordinates(lay) || (has_2d_coords(lay) && (_any_nonzero(lay.data.x2) || _any_nonzero(lay.data.y2)))
end

function _baseline_for_display(erps::Vector{ErpData})::Vector{ErpData}
    isempty(erps) && return erps
    out = Vector{ErpData}(undef, length(erps))
    for (i, erp) in enumerate(erps)
        tvec = time_vector(erp)
        needs_bl = !isempty(tvec) && first(tvec) < 0.0 && last(tvec) > 0.0 && any(tvec .<= 0.0)
        needs_rename = startswith(erp.condition_name, "grand_avg_")
        if needs_bl || needs_rename
            c = copy(erp)
            needs_rename && (c.condition_name = replace(c.condition_name, r"^grand_avg_" => ""))
            if needs_bl
                bl_mask = tvec .<= 0.0
                for lbl in c.layout.data.label
                    hasproperty(c.data, lbl) || continue
                    col = c.data[!, lbl]
                    eltype(col) <: Real || continue
                    m = mean(@view col[bl_mask])
                    isfinite(m) && (c.data[!, lbl] = col .- m)
                end
            end
            out[i] = c
        else
            out[i] = erp
        end
    end
    return out
end

function _quiet_gfp_vector(erp::ErpData)::Vector{Float64}
    chs = get_selected_channels(erp, channels(); include_meta = false, include_extra = false)
    n_t = nrow(erp.data)
    (isempty(chs) || n_t == 0) && return zeros(Float64, n_t)
    mat = Matrix{Float64}(undef, n_t, length(chs))
    @inbounds for (j, ch) in enumerate(chs)
        col = erp.data[!, ch]
        for i = 1:n_t
            mat[i, j] = Float64(col[i])
        end
    end
    out = Vector{Float64}(undef, n_t)
    @inbounds for i = 1:n_t
        out[i] = std(@view(mat[i, :]); corrected = false)
    end
    return out
end

function _frontal_roi_waveform(erps::Vector{ErpData})::Tuple{Vector{Float64},Vector{Float64},String}
    isempty(erps) && return Float64[], Float64[], ""
    disp = _baseline_for_display(erps)
    chs = get_selected_channels(disp[1], channels(); include_meta = false, include_extra = false)
    isempty(chs) && return Float64[], Float64[], ""
    roi = [ch for ch in _FRONTAL_ROI_CANDIDATES if ch in chs]
    isempty(roi) && (roi = first(chs, min(2, length(chs))))
    tvec = collect(Float64, time_vector(disp[1]))
    n_t = length(tvec)
    acc = zeros(Float64, n_t)
    n_used = 0
    for erp in disp
        nrow(erp.data) == n_t || continue
        for ch in roi
            hasproperty(erp.data, ch) || continue
            col = erp.data[!, ch]
            @inbounds for i = 1:n_t
                acc[i] += Float64(col[i])
            end
            n_used += 1
        end
    end
    n_used > 0 && (acc ./= n_used)
    lbl = length(roi) <= 3 ? join(roi, "/") : "$(first(roi))..$(last(roi)) ($(length(roi)) ch)"
    return tvec, acc, lbl
end

function _mean_gfp_across_conditions(erps::Vector{ErpData})::Tuple{Vector{Float64},Vector{Float64}}
    isempty(erps) && return Float64[], Float64[]
    disp = _baseline_for_display(erps)
    tvec = collect(Float64, time_vector(disp[1]))
    n_t = length(tvec)
    acc = zeros(Float64, n_t)
    n_used = 0
    for erp in disp
        nrow(erp.data) == n_t || continue
        acc .+= _quiet_gfp_vector(erp)
        n_used += 1
    end
    n_used > 0 && (acc ./= n_used)
    return tvec, acc
end

function _sensor_status_caption(rec::_ParticipantReportData)::String
    n_ep_ch = length(rec.epoch_channel_repairs)
    if !isempty(rec.repaired_channels) || !isempty(rec.skipped_channels)
        ep_part = n_ep_ch > 0 ? "; $n_ep_ch channel(s) interpolated across epochs" : ""
        method_part = rec.repair_method == "—" ? "" : " ($(rec.repair_method))"
        return "Continuous channel status$method_part$ep_part."
    elseif n_ep_ch > 0
        sorted_ep = sort!(collect(rec.epoch_channel_repairs); by = x -> (-x[2], string(x[1])))
        top_str = join(["$(ch):$(cnt)" for (ch, cnt) in first(sorted_ep, min(3, length(sorted_ep)))], ", ")
        return "No continuous bad channels; $n_ep_ch channel(s) interpolated across epochs (top: $top_str)."
    else
        return "No continuous or epoch-level channel interpolations."
    end
end

function _render_cohort_sensor_heatmap_uri(records::Vector{_ParticipantReportData})::Union{String,Nothing}
    lay = _cohort_layout(records)
    _has_spatial_layout(lay) || return nothing
    return try
        lay_copy = copy(lay)
        _ensure_coordinates_2d!(lay_copy)
        labels = lay_copy.data.label
        n_subj = length(records)

        cont_counts = zeros(Int, length(labels))
        ep_counts = zeros(Int, length(labels))
        lbl_to_idx = Dict{Symbol,Int}(lbl => i for (i, lbl) in enumerate(labels))

        for rec in records
            for ch in union(rec.repaired_channels, rec.skipped_channels)
                idx = get(lbl_to_idx, ch, 0)
                idx > 0 && (cont_counts[idx] += 1)
            end
            for (ch, cnt) in rec.epoch_channel_repairs
                idx = get(lbl_to_idx, ch, 0)
                idx > 0 && (ep_counts[idx] += cnt)
            end
        end

        (all(==(0), cont_counts) && all(==(0), ep_counts)) && return nothing

        fig = _eegfun_figure(size = (900, 360))
        ax_c = Axis(
            fig[1, 1],
            aspect = DataAspect(),
            title = "Continuous Repairs / Skips (Participants / $n_subj)",
            titlesize = 13,
            xautolimitmargin = (0.14, 0.16),
            yautolimitmargin = (0.14, 0.16),
        )
        plot_layout_2d!(fig, ax_c, lay_copy; point_markersize = 6, point_color = :lightgray, label_plot = false)
        hidedecorations!(ax_c)
        hidespines!(ax_c)

        nz_c = findall(>(0), cont_counts)
        max_c = isempty(nz_c) ? 1 : max(1, maximum(cont_counts))
        cb_max_c = max(2, max_c)
        if !isempty(nz_c)
            sc_c = scatter!(
                ax_c,
                lay_copy.data.x2[nz_c],
                lay_copy.data.y2[nz_c];
                color = Float64.(cont_counts[nz_c]),
                colormap = :YlOrRd,
                colorrange = (0.0, Float64(cb_max_c)),
                markersize = [12.0 + 5.0 * (cont_counts[i] / cb_max_c) for i in nz_c],
                strokecolor = :black,
                strokewidth = 1.2,
            )
            Colorbar(fig[1, 2], sc_c; width = 12, label = "Subjects", ticks = collect(0:max(1, cld(cb_max_c, 5)):cb_max_c))
            top_c = sort(nz_c; by = i -> -cont_counts[i])[1:min(5, length(nz_c))]
            text!(
                ax_c,
                lay_copy.data.x2[top_c],
                lay_copy.data.y2[top_c];
                text = ["$(labels[i]) ($(cont_counts[i]))" for i in top_c],
                fontsize = 9,
                font = :bold,
                color = :darkred,
                align = (:center, :bottom),
                offset = (0, 8),
            )
        end

        ax_e = Axis(
            fig[1, 3],
            aspect = DataAspect(),
            title = "Epoch-Level Interpolations (Total Epochs)",
            titlesize = 13,
            xautolimitmargin = (0.14, 0.16),
            yautolimitmargin = (0.14, 0.16),
        )
        plot_layout_2d!(fig, ax_e, lay_copy; point_markersize = 6, point_color = :lightgray, label_plot = false)
        hidedecorations!(ax_e)
        hidespines!(ax_e)

        nz_e = findall(>(0), ep_counts)
        max_e = isempty(nz_e) ? 1 : max(1, maximum(ep_counts))
        if !isempty(nz_e)
            sc_e = scatter!(
                ax_e,
                lay_copy.data.x2[nz_e],
                lay_copy.data.y2[nz_e];
                color = Float64.(ep_counts[nz_e]),
                colormap = :YlOrRd,
                colorrange = (0.0, Float64(max_e)),
                markersize = [9.0 + 7.0 * (ep_counts[i] / max_e) for i in nz_e],
                strokecolor = :saddlebrown,
                strokewidth = 1.1,
            )
            Colorbar(fig[1, 4], sc_e; width = 12, label = "Epochs")
            top_e = sort(nz_e; by = i -> -ep_counts[i])[1:min(5, length(nz_e))]
            text!(
                ax_e,
                lay_copy.data.x2[top_e],
                lay_copy.data.y2[top_e];
                text = ["$(labels[i]) ($(ep_counts[i]))" for i in top_e],
                fontsize = 9,
                font = :bold,
                color = :saddlebrown,
                align = (:center, :bottom),
                offset = (0, 8),
            )
        end

        _fig_to_data_uri(fig)
    catch e
        @debug "Skipping cohort sensor heatmap in pipeline report" exception = e
        nothing
    end
end

function _render_sensor_status_uri(rec::_ParticipantReportData)::Union{String,Nothing}
    lay = _participant_layout(rec)
    _has_spatial_layout(lay) || return nothing
    return try
        lay_copy = copy(lay)
        _ensure_coordinates_2d!(lay_copy)
        show_all_labels = nrow(lay_copy.data) <= 64
        cont_set = union(Set(rec.repaired_channels), Set(rec.skipped_channels))
        ep_only = [ch => cnt for (ch, cnt) in rec.epoch_channel_repairs if ch ∉ cont_set && cnt > 0]
        sort!(ep_only; by = x -> (-x[2], string(x[1])))

        title_str = if isempty(rec.repaired_channels) && isempty(rec.skipped_channels) && !isempty(ep_only)
            "Channel Status (0 cont., $(length(ep_only)) epoch interp.)"
        else
            "Channel Status ($(length(rec.repaired_channels)) repaired, $(length(rec.skipped_channels)) skipped)"
        end

        fig = _eegfun_figure(size = (400, 390))
        ax = Axis(
            fig[1, 1],
            aspect = DataAspect(),
            title = title_str,
            titlesize = 13,
            xautolimitmargin = (0.14, 0.16),
            yautolimitmargin = (0.14, 0.14),
        )
        plot_layout_2d!(
            fig,
            ax,
            lay_copy;
            point_markersize = show_all_labels ? 8 : 6,
            point_color = :slategray,
            label_plot = show_all_labels,
            label_fontsize = 9,
        )
        hidedecorations!(ax)
        hidespines!(ax)

        has_legend = false
        if !isempty(ep_only)
            ep_dict = Dict(ep_only)
            idx_ep = [i for (i, lbl) in enumerate(lay_copy.data.label) if haskey(ep_dict, lbl)]
            if !isempty(idx_ep)
                max_cnt = maximum(last, ep_only)
                sizes = [10.0 + 5.0 * (ep_dict[lay_copy.data.label[i]] / max(1, max_cnt)) for i in idx_ep]
                scatter!(
                    ax,
                    lay_copy.data.x2[idx_ep],
                    lay_copy.data.y2[idx_ep];
                    color = (:darkorange, 0.85),
                    markersize = sizes,
                    strokecolor = :saddlebrown,
                    strokewidth = 1.2,
                    label = "Epoch interp.",
                )
                has_legend = true
                if !show_all_labels
                    top_syms = Set(first.(first(ep_only, min(3, length(ep_only)))))
                    idx_top = [i for i in idx_ep if lay_copy.data.label[i] in top_syms]
                    if !isempty(idx_top)
                        text!(
                            ax,
                            lay_copy.data.x2[idx_top],
                            lay_copy.data.y2[idx_top];
                            text = ["$(lay_copy.data.label[i]) ($(ep_dict[lay_copy.data.label[i]]))" for i in idx_top],
                            color = :saddlebrown,
                            fontsize = 9,
                            font = :bold,
                            align = (:center, :bottom),
                            offset = (0, 7),
                        )
                    end
                end
            end
        end

        if !isempty(rec.repaired_channels)
            idx_rep = findall(in(Set(rec.repaired_channels)), lay_copy.data.label)
            if !isempty(idx_rep)
                scatter!(
                    ax,
                    lay_copy.data.x2[idx_rep],
                    lay_copy.data.y2[idx_rep];
                    color = :dodgerblue,
                    markersize = 15,
                    strokecolor = :navy,
                    strokewidth = 1.5,
                    label = "Cont. repaired",
                )
                has_legend = true
                if !show_all_labels
                    text!(
                        ax,
                        lay_copy.data.x2[idx_rep],
                        lay_copy.data.y2[idx_rep];
                        text = string.(lay_copy.data.label[idx_rep]),
                        color = :navy,
                        fontsize = 10,
                        font = :bold,
                        align = (:center, :bottom),
                        offset = (0, 8),
                    )
                end
            end
        end

        if !isempty(rec.skipped_channels)
            idx_skip = findall(in(Set(rec.skipped_channels)), lay_copy.data.label)
            if !isempty(idx_skip)
                scatter!(
                    ax,
                    lay_copy.data.x2[idx_skip],
                    lay_copy.data.y2[idx_skip];
                    color = :crimson,
                    marker = :xcross,
                    markersize = 17,
                    strokewidth = 2,
                    label = "Skipped",
                )
                has_legend = true
                if !show_all_labels
                    text!(
                        ax,
                        lay_copy.data.x2[idx_skip],
                        lay_copy.data.y2[idx_skip];
                        text = string.(lay_copy.data.label[idx_skip]),
                        color = :crimson,
                        fontsize = 10,
                        font = :bold,
                        align = (:center, :bottom),
                        offset = (0, 9),
                    )
                end
            end
        end

        has_legend && axislegend(ax; position = :rb, labelsize = 10)
        _fig_to_data_uri(fig)
    catch e
        @debug "Skipping sensor status figure in pipeline report" exception = e
        nothing
    end
end

function _render_ica_artifacts_uri(rec::_ParticipantReportData)::Union{String,Nothing}
    (isnothing(rec.ica) || isnothing(rec.artifact_info) || isnothing(rec.artifact_info.ica_components)) && return nothing
    (rec.ic_counts.total == 0 || !_has_spatial_layout(rec.ica.layout)) && return nothing
    return try
        res = plot_artifact_components(
            rec.ica,
            rec.artifact_info.ica_components;
            label_plot = false,
            point_markersize = 5,
            display_plot = false,
            interactive = false,
        )
        isempty(res.axes) ? nothing : _fig_to_data_uri(res.fig)
    catch e
        @debug "Skipping ICA artifacts figure in pipeline report" exception = e
        nothing
    end
end

function _render_erp_gfp_uri(
    erps::Vector{ErpData};
    plot_title::String = "Condition ERPs & GFP",
    fig_size::Tuple{Int,Int} = (880, 410),
)::Union{String,Nothing}
    isempty(erps) && return nothing
    return try
        disp_erps = _baseline_for_display(erps)
        res = plot_gfp(disp_erps; show_erp_traces = true, normalize = false, display_plot = false, plot_title = plot_title)
        resize!(res.fig, fig_size...)
        _fig_to_data_uri(res.fig)
    catch e
        @debug "Skipping ERP/GFP figure in pipeline report" exception = e
        nothing
    end
end

function _render_grand_topo_uri(grand_erps::Vector{ErpData})::Union{String,Nothing}
    (isempty(grand_erps) || !_has_spatial_layout(grand_erps[1].layout)) && return nothing
    disp_erps = _baseline_for_display(grand_erps)
    tvec = time_vector(disp_erps[1])
    isempty(tvec) && return nothing

    t_min, t_max = first(tvec), last(tvec)
    t_start = t_max > 0.0 ? max(0.0, t_min) : t_min
    t_end = (t_max > 0.0 && t_max >= 0.5 && t_start == 0.0) ? 0.5 : t_max
    t_end <= t_start && ((t_start, t_end) = (t_min, t_max))

    n_in_win = count(t -> t_start <= t <= t_end, tvec)
    n_bins = n_in_win >= 10 ? 5 : (n_in_win >= 4 ? 2 : 1)
    bin_edges = range(t_start, t_end, length = n_bins + 1)

    max_abs = 0.0
    for erp in disp_erps
        t_erp = time_vector(erp)
        for b = 1:n_bins
            mask = (t_erp .>= bin_edges[b]) .& (t_erp .<= bin_edges[b+1])
            any(mask) || continue
            for lbl in erp.layout.data.label
                hasproperty(erp.data, lbl) || continue
                col = erp.data[!, lbl]
                eltype(col) <: Real || continue
                v = abs(mean(@view col[mask]))
                isfinite(v) && v > max_abs && (max_abs = v)
            end
        end
    end
    lim = max_abs > 0.0 ? ceil(max_abs * 1.05; digits = 2) : 1.0
    cr = (-lim, lim)

    return try
        res = if n_bins > 1
            plot_topography(
                disp_erps;
                interval_selection = (t_start, t_end),
                n_topo = n_bins,
                time_unit = :ms,
                label_plot = false,
                point_markersize = 4,
                head_linewidth = 1.5,
                plot_title_fontsize = 13,
                colorrange = cr,
                colorbar_plot_numbers = collect(1:length(disp_erps)),
                interactive = false,
                display_plot = false,
            )
        else
            plot_topography(
                disp_erps;
                interval_selection = (t_start, t_end),
                label_plot = false,
                point_markersize = 5,
                head_linewidth = 1.5,
                plot_title_fontsize = 14,
                colorrange = cr,
                interactive = false,
                display_plot = false,
            )
        end
        isnothing(res) && return nothing
        if n_bins > 1 && length(res.axes) == length(disp_erps) * n_bins
            for b = 1:(n_bins-1)
                2 * b <= size(res.fig.layout)[2] && colsize!(res.fig.layout, 2 * b, Fixed(6))
            end
            for r in eachindex(disp_erps)
                for b = 1:n_bins
                    ax = res.axes[(r-1)*n_bins+b]
                    ax.title = r == 1 ? @sprintf("%.0f – %.0f ms", bin_edges[b] * 1000, bin_edges[b+1] * 1000) : ""
                    ax.titlesize = 13
                end
                raw_cond = strip(disp_erps[r].condition_name)
                cond_lbl = isempty(raw_cond) ? "Condition $(disp_erps[r].condition)" : replace(raw_cond, r"^grand_avg_" => "")
                Label(res.fig[r, 0], cond_lbl; rotation = pi / 2, fontsize = 13, font = :bold, tellheight = false)
            end
        end
        w = n_bins > 1 ? max(820, 215 * n_bins + 80) : max(440, 360 * length(disp_erps))
        h = n_bins > 1 ? max(260, 235 * length(disp_erps)) : 320
        resize!(res.fig, w, h)
        _fig_to_data_uri(res.fig)
    catch e
        @debug "Skipping grand-average topography figure in pipeline report" exception = e
        nothing
    end
end

function _cohort_mean_psd(records::Vector{_ParticipantReportData})::Union{Tuple{Vector{Float64},Vector{Float64}},Nothing}
    psds = [r.psd for r in records if !isnothing(r.psd)]
    isempty(psds) && return nothing
    f_ref = psds[1].freqs
    matching = [p.mean_db for p in psds if p.freqs == f_ref]
    isempty(matching) && return nothing
    return f_ref, vec(mean(reduce(hcat, matching); dims = 2))
end

function _render_cohort_psd_uri(records::Vector{_ParticipantReportData})::Union{String,Nothing}
    psd_recs = [r for r in records if !isnothing(r.psd)]
    isempty(psd_recs) && return nothing
    return try
        f_ref = psd_recs[1].psd.freqs
        same_grid = [r for r in psd_recs if r.psd.freqs == f_ref]
        isempty(same_grid) && return nothing
        kind_lbl = any(r -> r.psd.is_epoch_level, same_grid) ? "Cleaned Single-Trial Welch PSD" : "Evoked ERP PSD"

        fig = _eegfun_figure(size = (1040, 330))
        ax1 =
            Axis(fig[1, 1], title = "$kind_lbl — Cohort Spaghetti", titlesize = 13, xlabel = "Frequency (Hz)", ylabel = "Power (dB μV²/Hz)")
        for r in same_grid
            lines!(
                ax1,
                f_ref,
                r.psd.mean_db;
                color = r.is_flagged ? (:crimson, 0.45) : (:steelblue, 0.25),
                linewidth = r.is_flagged ? 1.2 : 0.9,
                linestyle = r.is_flagged ? :dash : :solid,
            )
        end
        mean_mat = reduce(hcat, [r.psd.mean_db for r in same_grid])
        cohort_m = vec(mean(mean_mat; dims = 2))
        lines!(ax1, f_ref, cohort_m; color = :navy, linewidth = 2.4, label = "Cohort Mean")
        10.0 <= last(f_ref) && vlines!(ax1, [10.0]; color = (:seagreen, 0.6), linestyle = :dot, linewidth = 1.3, label = "10 Hz (Alpha)")
        50.0 <= last(f_ref) && vlines!(ax1, [50.0]; color = (:darkorange, 0.6), linestyle = :dot, linewidth = 1.3, label = "50 Hz (Line)")
        axislegend(ax1; position = :rt, labelsize = 10)

        post_vecs = [r.psd.posterior_db for r in same_grid if !isnothing(r.psd.posterior_db)]
        front_vecs = [r.psd.frontal_db for r in same_grid if !isnothing(r.psd.frontal_db)]
        if !isempty(post_vecs) || !isempty(front_vecs)
            ax2 = Axis(
                fig[1, 2],
                title = "Regional PSD (Mean ± 1 SEM)",
                titlesize = 13,
                xlabel = "Frequency (Hz)",
                ylabel = "Power (dB μV²/Hz)",
            )
            lines!(ax2, f_ref, cohort_m; color = (:slategray, 0.7), linewidth = 1.6, linestyle = :dash, label = "Whole Scalp")
            if !isempty(post_vecs)
                pmat = reduce(hcat, post_vecs)
                pm = vec(mean(pmat; dims = 2))
                if length(post_vecs) >= 2
                    psem = vec(std(pmat; dims = 2)) ./ sqrt(length(post_vecs))
                    band!(ax2, f_ref, pm .- psem, pm .+ psem; color = (:seagreen, 0.2))
                end
                lines!(ax2, f_ref, pm; color = :seagreen, linewidth = 2.2, label = "Posterior (O/PO/P)")
            end
            if !isempty(front_vecs)
                fmat = reduce(hcat, front_vecs)
                fm = vec(mean(fmat; dims = 2))
                if length(front_vecs) >= 2
                    fsem = vec(std(fmat; dims = 2)) ./ sqrt(length(front_vecs))
                    band!(ax2, f_ref, fm .- fsem, fm .+ fsem; color = (:royalblue, 0.2))
                end
                lines!(ax2, f_ref, fm; color = :royalblue, linewidth = 2.2, label = "Frontal (Fp/AF/F)")
            end
            10.0 <= last(f_ref) && vlines!(ax2, [10.0]; color = (:gray, 0.5), linestyle = :dot, linewidth = 1.2)
            50.0 <= last(f_ref) && vlines!(ax2, [50.0]; color = (:gray, 0.5), linestyle = :dot, linewidth = 1.2)
            axislegend(ax2; position = :rt, labelsize = 10)
        end

        _fig_to_data_uri(fig)
    catch e
        @debug "Skipping cohort PSD figure in pipeline report" exception = e
        nothing
    end
end

function _render_participant_psd_uri(
    rec::_ParticipantReportData,
    cohort_psd::Union{Tuple{Vector{Float64},Vector{Float64}},Nothing},
)::Union{String,Nothing}
    isnothing(rec.psd) && return nothing
    return try
        p = rec.psd
        fig = _eegfun_figure(size = (420, 320))
        ax = Axis(
            fig[1, 1],
            title = p.is_epoch_level ? "Cleaned Welch PSD" : "Evoked ERP PSD",
            titlesize = 13,
            xlabel = "Frequency (Hz)",
            ylabel = "Power (dB μV²/Hz)",
        )
        if !isnothing(cohort_psd) && cohort_psd[1] == p.freqs
            lines!(ax, p.freqs, cohort_psd[2]; color = (:slategray, 0.7), linestyle = :dash, linewidth = 1.5, label = "Cohort Mean")
        end
        !isnothing(p.posterior_db) && lines!(ax, p.freqs, p.posterior_db; color = :seagreen, linewidth = 1.8, label = "Posterior")
        !isnothing(p.frontal_db) && lines!(ax, p.freqs, p.frontal_db; color = :darkorange, linewidth = 1.8, label = "Frontal")
        lines!(ax, p.freqs, p.mean_db; color = :navy, linewidth = 2.1, label = "All Channels")
        10.0 <= last(p.freqs) && vlines!(ax, [10.0]; color = (:gray, 0.45), linestyle = :dot, linewidth = 1.1)
        50.0 <= last(p.freqs) && vlines!(ax, [50.0]; color = (:gray, 0.45), linestyle = :dot, linewidth = 1.1)
        axislegend(ax; position = :rt, labelsize = 9)
        _fig_to_data_uri(fig)
    catch e
        @debug "Skipping participant PSD figure in pipeline report" exception = e
        nothing
    end
end

function _matched_raw_clean_erps(rec::_ParticipantReportData)::Tuple{Vector{ErpData},Vector{ErpData}}
    (isnothing(rec.erps_raw) || isnothing(rec.erps) || isempty(rec.erps_raw) || isempty(rec.erps)) && return ErpData[], ErpData[]
    common = Set(intersect([e.condition for e in rec.erps_raw], [e.condition for e in rec.erps]))
    isempty(common) && return rec.erps_raw, rec.erps
    return [e for e in rec.erps_raw if e.condition in common], [e for e in rec.erps if e.condition in common]
end

function _render_cohort_raw_vs_clean_uri(records::Vector{_ParticipantReportData})::Union{String,Nothing}
    paired = [r for r in records if !isnothing(r.erps_raw) && !isnothing(r.erps) && !isempty(r.erps_raw) && !isempty(r.erps)]
    isempty(paired) && return nothing
    return try
        r1_raw, r1_cln = _matched_raw_clean_erps(paired[1])
        t_ref, raw1_f, roi_lbl = _frontal_roi_waveform(r1_raw)
        _, cln1_f, _ = _frontal_roi_waveform(r1_cln)
        (isempty(t_ref) || length(raw1_f) != length(cln1_f)) && return nothing
        n_t = length(t_ref)

        raw_front_list = Vector{Float64}[]
        cln_front_list = Vector{Float64}[]
        raw_gfp_list = Vector{Float64}[]
        cln_gfp_list = Vector{Float64}[]

        for r in paired
            r_raw, r_cln = _matched_raw_clean_erps(r)
            tr, rf, _ = _frontal_roi_waveform(r_raw)
            tc, cf, _ = _frontal_roi_waveform(r_cln)
            _, rg = _mean_gfp_across_conditions(r_raw)
            _, cg = _mean_gfp_across_conditions(r_cln)
            if length(tr) == n_t && length(tc) == n_t && length(rg) == n_t && length(cg) == n_t
                push!(raw_front_list, rf)
                push!(cln_front_list, cf)
                push!(raw_gfp_list, rg)
                push!(cln_gfp_list, cg)
            end
        end
        isempty(raw_front_list) && return nothing

        fig = _eegfun_figure(size = (1040, 330))
        ax1 = Axis(
            fig[1, 1],
            title = "Frontal Ocular ROI ($roi_lbl) — Raw vs. Cleaned",
            titlesize = 13,
            xlabel = "Time (s)",
            ylabel = "Amplitude (μV)",
        )
        rf_mat = reduce(hcat, raw_front_list)
        cf_mat = reduce(hcat, cln_front_list)
        rf_m = vec(mean(rf_mat; dims = 2))
        cf_m = vec(mean(cf_mat; dims = 2))
        if length(raw_front_list) >= 2
            rf_s = vec(std(rf_mat; dims = 2)) ./ sqrt(length(raw_front_list))
            cf_s = vec(std(cf_mat; dims = 2)) ./ sqrt(length(cln_front_list))
            band!(ax1, t_ref, rf_m .- rf_s, rf_m .+ rf_s; color = (:crimson, 0.16))
            band!(ax1, t_ref, cf_m .- cf_s, cf_m .+ cf_s; color = (:dodgerblue, 0.2))
        end
        vlines!(ax1, [0.0]; color = :black, linestyle = :dash, linewidth = 1)
        lines!(ax1, t_ref, rf_m; color = :crimson, linestyle = :dash, linewidth = 1.9, label = "Raw (Pre-ICA)")
        lines!(ax1, t_ref, cf_m; color = :dodgerblue, linewidth = 2.2, label = "Cleaned (Post-ICA)")
        axislegend(ax1; position = :rt, labelsize = 10)

        ax2 = Axis(fig[1, 2], title = "Global Field Power — Raw vs. Cleaned", titlesize = 13, xlabel = "Time (s)", ylabel = "GFP (μV)")
        rg_mat = reduce(hcat, raw_gfp_list)
        cg_mat = reduce(hcat, cln_gfp_list)
        rg_m = vec(mean(rg_mat; dims = 2))
        cg_m = vec(mean(cg_mat; dims = 2))
        if length(raw_gfp_list) >= 2
            rg_s = vec(std(rg_mat; dims = 2)) ./ sqrt(length(raw_gfp_list))
            cg_s = vec(std(cg_mat; dims = 2)) ./ sqrt(length(cln_gfp_list))
            band!(ax2, t_ref, rg_m .- rg_s, rg_m .+ rg_s; color = (:crimson, 0.16))
            band!(ax2, t_ref, cg_m .- cg_s, cg_m .+ cg_s; color = (:dodgerblue, 0.2))
        end
        vlines!(ax2, [0.0]; color = :black, linestyle = :dash, linewidth = 1)
        lines!(ax2, t_ref, rg_m; color = :crimson, linestyle = :dash, linewidth = 1.9, label = "Raw (Pre-ICA)")
        lines!(ax2, t_ref, cg_m; color = :dodgerblue, linewidth = 2.2, label = "Cleaned (Post-ICA)")
        axislegend(ax2; position = :rt, labelsize = 10)

        _fig_to_data_uri(fig)
    catch e
        @debug "Skipping cohort raw-vs-clean figure in pipeline report" exception = e
        nothing
    end
end

function _render_participant_raw_vs_clean_uri(rec::_ParticipantReportData)::Union{String,Nothing}
    (isnothing(rec.erps_raw) || isnothing(rec.erps) || isempty(rec.erps_raw) || isempty(rec.erps)) && return nothing
    return try
        r_raw, r_cln = _matched_raw_clean_erps(rec)
        tr, rf, roi_lbl = _frontal_roi_waveform(r_raw)
        tc, cf, _ = _frontal_roi_waveform(r_cln)
        _, rg = _mean_gfp_across_conditions(r_raw)
        _, cg = _mean_gfp_across_conditions(r_cln)
        (isempty(tr) || length(tr) != length(tc) || length(rg) != length(tr) || length(cg) != length(tr)) && return nothing

        fig = _eegfun_figure(size = (880, 320))
        ax1 = Axis(
            fig[1, 1],
            title = "Frontal ROI ($roi_lbl): Raw vs. Cleaned",
            titlesize = 13,
            xlabel = "Time (s)",
            ylabel = "Amplitude (μV)",
        )
        vlines!(ax1, [0.0]; color = :black, linestyle = :dash, linewidth = 1)
        lines!(ax1, tr, rf; color = :crimson, linestyle = :dash, linewidth = 1.7, label = "Raw (Pre-ICA)")
        lines!(ax1, tc, cf; color = :dodgerblue, linewidth = 2.0, label = "Cleaned")
        axislegend(ax1; position = :rt, labelsize = 9)

        ax2 = Axis(fig[1, 2], title = "Mean GFP: Raw vs. Cleaned", titlesize = 13, xlabel = "Time (s)", ylabel = "GFP (μV)")
        vlines!(ax2, [0.0]; color = :black, linestyle = :dash, linewidth = 1)
        lines!(ax2, tr, rg; color = :crimson, linestyle = :dash, linewidth = 1.7, label = "Raw (Pre-ICA)")
        lines!(ax2, tc, cg; color = :dodgerblue, linewidth = 2.0, label = "Cleaned")
        axislegend(ax2; position = :rt, labelsize = 9)

        _fig_to_data_uri(fig)
    catch e
        @debug "Skipping participant raw-vs-clean figure in pipeline report" exception = e
        nothing
    end
end

function _format_ic_breakdown(ic)::String
    ic.total == 0 && return "0"
    parts = String[]
    ic.vEOG > 0 && push!(parts, "vEOG:$(ic.vEOG)")
    ic.hEOG > 0 && push!(parts, "hEOG:$(ic.hEOG)")
    ic.ECG > 0 && push!(parts, "ECG:$(ic.ECG)")
    ic.line_noise > 0 && push!(parts, "Line:$(ic.line_noise)")
    ic.channel_noise > 0 && push!(parts, "Chan:$(ic.channel_noise)")
    ic.manual > 0 && push!(parts, "Manual:$(ic.manual)")
    cat_sum = ic.vEOG + ic.hEOG + ic.ECG + ic.line_noise + ic.channel_noise + ic.manual
    cat_sum > ic.total && push!(parts, "$(cat_sum - ic.total) shared")
    return isempty(parts) ? string(ic.total) : "$(ic.total) ($(join(parts, ", ")))"
end

function _resolve_avg_neighbors(
    output_dir::String,
    cfg::Union{PipelineConfig,Nothing},
    records::Vector{_ParticipantReportData},
)::Union{Float64,Nothing}
    nb_path = joinpath(output_dir, "neighbors.toml")
    if isfile(nb_path)
        try
            nb_toml = TOML.parsefile(nb_path)
            meta = get(nb_toml, "metadata", nothing)
            if meta isa AbstractDict && haskey(meta, "average_neighbors_per_electrode")
                val = Float64(meta["average_neighbors_per_electrode"])
                isfinite(val) && return round(val; digits = 2)
            end
        catch
        end
    end

    lay = _cohort_layout(records)
    if !isnothing(lay)
        if has_neighbors(lay) && !isempty(lay.neighbors)
            return round(mean(length(nb.channels) for nb in values(lay.neighbors)); digits = 2)
        elseif !isnothing(cfg) && _has_spatial_layout(lay)
            try
                lay_tmp = copy(lay)
                _ensure_coordinates_2d!(lay_tmp)
                get_neighbors_xy!(lay_tmp, cfg.preprocess.neighbor_criterion)
                if has_neighbors(lay_tmp) && !isempty(lay_tmp.neighbors)
                    return round(mean(length(nb.channels) for nb in values(lay_tmp.neighbors)); digits = 2)
                end
            catch
            end
        end
    end
    return nothing
end

function _build_methods_boilerplate(
    cfg::Union{PipelineConfig,Nothing},
    records::Vector{_ParticipantReportData},
    ver_str::String,
    low_epoch_threshold::Float64,
    max_ica_pct::Float64,
    avg_neighbors::Union{Float64,Nothing} = nothing,
)::String
    n_subj = length(records)
    rep_counts = [length(r.repaired_channels) for r in records]
    skip_counts = [length(r.skipped_channels) for r in records]
    ic_totals = [r.ic_counts.total for r in records]
    valid_pcts = [r.overall_kept_pct for r in records if !isnan(r.overall_kept_pct)]

    rep_m = round(mean(rep_counts); digits = 2)
    rep_sd = n_subj > 1 ? round(std(rep_counts); digits = 2) : 0.0
    ic_m = round(mean(ic_totals); digits = 2)
    ic_sd = n_subj > 1 ? round(std(ic_totals); digits = 2) : 0.0

    parts = String[]
    if !isnothing(cfg)
        p = cfg.preprocess
        filt_bits = String[]
        p.filter.highpass.apply && push!(
            filt_bits,
            "high-pass filtered at $(p.filter.highpass.freq) Hz ($(p.filter.highpass.method), order $(p.filter.highpass.order))",
        )
        p.filter.lowpass.apply && push!(
            filt_bits,
            "low-pass filtered at $(p.filter.lowpass.freq) Hz ($(p.filter.lowpass.method), order $(p.filter.lowpass.order))",
        )
        p.resample.apply && push!(filt_bits, "resampled to $(p.resample.target_rate) Hz")
        p.cleanline.apply && push!(filt_bits, "cleaned of $(join(p.cleanline.line_frequencies, "/")) Hz line noise via CleanLine")
        filt_clause = isempty(filt_bits) ? "imported without continuous frequency filtering" : join(filt_bits, ", ")
        ref_str = p.reference_channel isa AbstractVector ? join(p.reference_channel, ", ") : string(p.reference_channel)
        nb_clause =
            isnothing(avg_neighbors) ? "neighbor criterion = $(p.neighbor_criterion)" :
            "neighbor criterion = $(p.neighbor_criterion), avg neighbors = $avg_neighbors"

        push!(
            parts,
            "EEG recordings (N = $n_subj; layout: $(cfg.files.input.layout_file)) were preprocessed in Julia using EegFun.jl (v$ver_str). Continuous signals were $filt_clause, bad channels were interpolated via $(p.channel_repair_method) ($nb_clause; M = $rep_m ± $rep_sd channels/participant, range $(minimum(rep_counts))–$(maximum(rep_counts)); $(sum(skip_counts)) skipped), and data were re-referenced to $ref_str.",
        )

        if p.ica.apply
            ica_filts = String[]
            p.filter.ica_highpass.apply && push!(ica_filts, "$(p.filter.ica_highpass.freq) Hz high-pass")
            p.filter.ica_lowpass.apply && push!(ica_filts, "$(p.filter.ica_lowpass.freq) Hz low-pass")
            ica_filt_str = isempty(ica_filts) ? "" : " after $(join(ica_filts, " and ")) filtering"
            push!(
                parts,
                "Independent Component Analysis (Infomax ICA$ica_filt_str; $(p.ica.component_method) artifact identification on $(p.ica.percentage_of_data)% of clean samples) removed M = $ic_m ± $ic_sd artifact components per participant (range $(minimum(ic_totals))–$(maximum(ic_totals))).",
            )
        end

        rej_desc = if p.autoreject.apply
            "Autoreject ($(p.autoreject.mode) mode, $(p.autoreject.n_folds) folds)"
        else
            "amplitude threshold screening (±$(p.eeg.artifact_value_abs_criterion) μV; extreme ±$(p.eeg.extreme_value_abs_criterion) μV) with per-epoch channel interpolation"
        end
        push!(
            parts,
            "Data were segmented from $(p.epoch_start) s to $(p.epoch_end) s relative to event onset, baseline-corrected, and cleaned via $rej_desc.",
        )
    else
        push!(
            parts,
            "EEG datasets (N = $n_subj) were summarized using EegFun.jl (v$ver_str), with M = $rep_m ± $rep_sd continuous channels interpolated (range $(minimum(rep_counts))–$(maximum(rep_counts))) and M = $ic_m ± $ic_sd ICA components removed (range $(minimum(ic_totals))–$(maximum(ic_totals))).",
        )
    end

    if !isempty(valid_pcts)
        pct_m = round(mean(valid_pcts); digits = 1)
        pct_sd = length(valid_pcts) > 1 ? round(std(valid_pcts); digits = 1) : 0.0
        cond_names = String[]
        cond_kept = Dict{String,Vector{Int}}()
        cond_pct = Dict{String,Vector{Float64}}()
        for r in records, row in eachrow(r.condition_df)
            cn = row.condition_name
            if !haskey(cond_kept, cn)
                push!(cond_names, cn)
                cond_kept[cn] = Int[]
                cond_pct[cn] = Float64[]
            end
            push!(cond_kept[cn], row.n_kept)
            push!(cond_pct[cn], row.pct_kept)
        end
        cond_clauses = [
            "$(cn): M = $(round(mean(cond_kept[cn]); digits = 1)) ± $(length(cond_kept[cn]) > 1 ? round(std(cond_kept[cn]); digits = 1) : 0.0) clean epochs ($(round(mean(cond_pct[cn]); digits = 1))%)"
            for cn in cond_names
        ]
        cond_tail = isempty(cond_clauses) ? "" : " ($(join(cond_clauses, "; ")))"
        push!(
            parts,
            "Overall epoch retention averaged $(pct_m)% ± $(pct_sd)% across participants (range $(minimum(valid_pcts))%–$(maximum(valid_pcts))%)$cond_tail.",
        )
    end

    n_flagged = count(r -> r.is_flagged, records)
    push!(
        parts,
        "Automated quality-control screening (epoch retention < $(round(low_epoch_threshold; digits = 1))% or ICA removal > $(round(max_ica_pct; digits = 1))%) flagged $n_flagged / $n_subj participant(s) for manual inspection.",
    )
    return join(parts, " ")
end

function _participant_n_channels(r::_ParticipantReportData)::Union{Int,Nothing}
    if !isnothing(r.erps) && !isempty(r.erps)
        n_lay = nrow(r.erps[1].layout.data)
        return n_lay > 0 ? n_lay : length(channel_labels(r.erps[1]))
    elseif !isnothing(r.ica)
        return nrow(r.ica.layout.data)
    end
    return nothing
end

function _participant_ic_available(r::_ParticipantReportData)::Union{Int,Nothing}
    if !isnothing(r.ica)
        return size(r.ica.unmixing, 1)
    elseif !isnothing(r.artifact_info) && !isnothing(r.artifact_info.ica_components)
        return _participant_n_channels(r)
    end
    return nothing
end

function _cohort_condition_columns(records::Vector{_ParticipantReportData})::Vector{Tuple{Int,String}}
    cond_names = OrderedDict{Int,String}()
    for r in records, row in eachrow(r.condition_df)
        cid = Int(row.condition)
        cname = strip(string(row.condition_name))
        if !haskey(cond_names, cid) || (isempty(cond_names[cid]) && !isempty(cname))
            cond_names[cid] = cname
        end
    end
    cids = sort!(collect(keys(cond_names)))
    raw_slugs = Dict{Int,String}()
    slug_counts = Dict{String,Int}()
    for cid in cids
        s = strip(replace(cond_names[cid], r"[^A-Za-z0-9]+" => "_"), '_')
        isempty(s) && (s = "cond_$cid")
        raw_slugs[cid] = s
        slug_counts[s] = get(slug_counts, s, 0) + 1
    end
    return Tuple{Int,String}[(cid, slug_counts[raw_slugs[cid]] > 1 ? "$(raw_slugs[cid])_$cid" : raw_slugs[cid]) for cid in cids]
end

function _build_qc_summary_csv(records::Vector{_ParticipantReportData})::String
    cond_cols = _cohort_condition_columns(records)
    headers = String["participant", "status", "epochs_original", "epochs_repaired", "epochs_rejected", "epochs_kept", "epochs_kept_pct"]
    for (_, slug) in cond_cols
        push!(headers, "$(slug)_original", "$(slug)_kept", "$(slug)_kept_pct")
    end
    append!(
        headers,
        (
            "n_channels",
            "n_repaired_channels",
            "n_skipped_channels",
            "n_epoch_interp_channels",
            "n_epoch_channel_interps",
            "ic_available",
            "ic_removed",
            "ic_removed_pct",
            "ic_vEOG",
            "ic_hEOG",
            "ic_ECG",
            "ic_line_noise",
            "ic_channel_noise",
            "ic_manual",
            "repaired_channels",
            "skipped_channels",
            "top_epoch_interp_channels",
            "flag_reasons",
        ),
    )

    io = IOBuffer()
    println(io, join(headers, ","))
    fields = Vector{String}(undef, length(headers))

    for r in records
        k = 0
        fields[k+=1] = _csv_escape(r.stem)
        fields[k+=1] = r.is_flagged ? "Review" : "OK"

        has_epochs = !isempty(r.condition_df) || r.total_original > 0
        n_rep_ep = isempty(r.condition_df) ? 0 : sum(r.condition_df.n_repaired)
        fields[k+=1] = has_epochs ? string(r.total_original) : ""
        fields[k+=1] = has_epochs ? string(n_rep_ep) : ""
        fields[k+=1] = has_epochs ? string(r.total_original - r.total_kept) : ""
        fields[k+=1] = has_epochs ? string(r.total_kept) : ""
        fields[k+=1] = r.total_original > 0 ? string(round(100.0 * r.total_kept / r.total_original; digits = 1)) : ""

        cond_lookup = Dict{Int,Int}()
        for (row_i, cid) in enumerate(r.condition_df.condition)
            haskey(cond_lookup, cid) || (cond_lookup[cid] = row_i)
        end
        for (cid, _) in cond_cols
            row_i = get(cond_lookup, cid, 0)
            if row_i > 0
                fields[k+=1] = string(r.condition_df.n_original[row_i])
                fields[k+=1] = string(r.condition_df.n_kept[row_i])
                fields[k+=1] = string(r.condition_df.pct_kept[row_i])
            else
                fields[k+=1] = ""
                fields[k+=1] = ""
                fields[k+=1] = ""
            end
        end

        n_ch = _participant_n_channels(r)
        n_ep_ch = length(r.epoch_channel_repairs)
        n_ep_interps = isempty(r.epoch_channel_repairs) ? 0 : sum(values(r.epoch_channel_repairs))
        fields[k+=1] = isnothing(n_ch) ? "" : string(n_ch)
        fields[k+=1] = string(length(r.repaired_channels))
        fields[k+=1] = string(length(r.skipped_channels))
        fields[k+=1] = string(n_ep_ch)
        fields[k+=1] = string(n_ep_interps)

        ic_avail = _participant_ic_available(r)
        has_ica = !isnothing(ic_avail)
        fields[k+=1] = has_ica ? string(ic_avail) : ""
        fields[k+=1] = has_ica ? string(r.ic_counts.total) : (r.ic_counts.total > 0 ? string(r.ic_counts.total) : "0")
        fields[k+=1] = (has_ica && ic_avail > 0) ? string(round(100.0 * r.ic_counts.total / ic_avail; digits = 1)) : ""
        fields[k+=1] = string(r.ic_counts.vEOG)
        fields[k+=1] = string(r.ic_counts.hEOG)
        fields[k+=1] = string(r.ic_counts.ECG)
        fields[k+=1] = string(r.ic_counts.line_noise)
        fields[k+=1] = string(r.ic_counts.channel_noise)
        fields[k+=1] = string(r.ic_counts.manual)

        ep_sorted = sort!(collect(r.epoch_channel_repairs); by = x -> (-x[2], string(x[1])))
        top_n = min(5, length(ep_sorted))
        top_ep_str = join(["$(ep_sorted[i][1]):$(ep_sorted[i][2])" for i = 1:top_n], "; ")
        length(ep_sorted) > 5 && (top_ep_str *= " (+$(length(ep_sorted) - 5) more)")

        fields[k+=1] = _csv_escape(join(r.repaired_channels, "; "))
        fields[k+=1] = _csv_escape(join(r.skipped_channels, "; "))
        fields[k+=1] = _csv_escape(top_ep_str)
        fields[k+=1] = _csv_escape(join(r.flag_reasons, "; "))

        println(io, join(fields, ","))
    end
    return String(take!(io))
end

function _write_config_and_methods_section!(
    io::IOBuffer,
    output_dir::String,
    cfg::Union{PipelineConfig,Nothing},
    records::Vector{_ParticipantReportData},
    ver_str::String,
    low_epoch_threshold::Float64,
    max_ica_pct::Float64,
)
    avg_nb = _resolve_avg_neighbors(output_dir, cfg, records)
    methods_txt = _build_methods_boilerplate(cfg, records, ver_str, low_epoch_threshold, max_ica_pct, avg_nb)
    write(io, "<section class=\"panel\">")
    if !isnothing(cfg)
        p = cfg.preprocess
        hp_str =
            p.filter.highpass.apply ? "$(p.filter.highpass.freq) Hz ($(p.filter.highpass.method), order $(p.filter.highpass.order))" : "off"
        lp_str =
            p.filter.lowpass.apply ? "$(p.filter.lowpass.freq) Hz ($(p.filter.lowpass.method), order $(p.filter.lowpass.order))" : "off"
        ica_str =
            p.ica.apply ?
            "on ($(p.ica.component_method), $(p.ica.percentage_of_data)% data, HP $(p.filter.ica_highpass.apply ? "$(p.filter.ica_highpass.freq) Hz" : "off"))" :
            "off"
        cl_str = p.cleanline.apply ? "on ($(join(p.cleanline.line_frequencies, ", ")) Hz)" : "off"
        rs_str = p.resample.apply ? "on ($(p.resample.target_rate) Hz)" : "off"
        rej_str = if p.autoreject.apply
            "Autoreject (mode=$(p.autoreject.mode), folds=$(p.autoreject.n_folds), κ=$(p.autoreject.consensus_perc), ρ=$(p.autoreject.max_repair_perc))"
        else
            z_part = p.eeg.artifact_value_z_criterion > 0 ? ", z=$(p.eeg.artifact_value_z_criterion)" : ""
            "Threshold (abs=$(p.eeg.artifact_value_abs_criterion) μV$z_part, extreme=$(p.eeg.extreme_value_abs_criterion) μV)"
        end
        ref_str = p.reference_channel isa AbstractVector ? join(p.reference_channel, ", ") : string(p.reference_channel)
        nb_str =
            isnothing(avg_nb) ? "neighbor criterion = $(p.neighbor_criterion)" :
            "neighbor criterion = $(p.neighbor_criterion), avg neighbors = $avg_nb"

        write(io, "<h2>Pipeline Configuration &amp; Methods Summary</h2><div class=\"cfg-grid\">")
        for (k, v) in (
            ("Layout", cfg.files.input.layout_file),
            ("Reference", ref_str),
            ("Epoch Window", "$(p.epoch_start) s to $(p.epoch_end) s"),
            ("High-Pass Filter", hp_str),
            ("Low-Pass Filter", lp_str),
            ("Resampling", rs_str),
            ("CleanLine", cl_str),
            ("ICA", ica_str),
            ("Epoch Rejection", rej_str),
            ("Channel Repair", "$(p.channel_repair_method) ($nb_str)"),
            ("EOG Thresholds", "vEOG = $(p.eog.vEOG_criterion) μV, hEOG = $(p.eog.hEOG_criterion) μV"),
        )
            write(io, "<div class=\"cfg-item\"><span class=\"k\">$(_html_escape(k)):</span> $(_html_escape(v))</div>")
        end
        write(io, "</div>")
    else
        write(io, "<h2>Methods &amp; Attrition Summary</h2>")
    end

    write(
        io,
        """
  <div class="methods-box">
    <div class="methods-head">
      <span>Draft Methods &amp; Attrition Boilerplate</span>
      <button type="button" class="btn-ctrl" onclick="copyMethodsText(this)">Copy Methods Text</button>
    </div>
    <p id="methods-boilerplate">$(_html_escape(methods_txt))</p>
  </div>
</section>
""",
    )
end

function _build_pipeline_report_html(
    output_dir::String,
    title::String,
    low_epoch_threshold::Float64,
    max_ica_pct::Float64,
    pipeline_cfg::Union{PipelineConfig,Nothing},
    records::Vector{_ParticipantReportData},
    electrode_freq_df::DataFrame,
    can_plot::Bool,
    csv_text::String = _build_qc_summary_csv(records),
)::String
    n_subj = length(records)
    valid_pcts = [r.overall_kept_pct for r in records if !isnan(r.overall_kept_pct)]
    mean_kept = isempty(valid_pcts) ? NaN : round(mean(valid_pcts); digits = 1)
    n_flagged = count(r -> r.is_flagged, records)
    mean_rep = round(mean(r -> length(r.repaired_channels), records); digits = 2)
    mean_ics = round(mean(r -> r.ic_counts.total, records); digits = 2)
    ver_str = string(pkgversion(@__MODULE__))
    ts_str = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    csv_uri = "data:text/csv;charset=utf-8;base64," * _base64_encode(codeunits(csv_text))
    slugs = _unique_slugs(records)

    io = IOBuffer(; sizehint = 65536)
    write(
        io,
        """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>$(_html_escape(title))</title>
<style>
$_REPORT_CSS
</style>
</head>
<body>
<div class="container">
<header class="hero">
  <h1>$(_html_escape(title))</h1>
  <div class="meta">Directory: <code>$(_html_escape(abspath(output_dir)))</code> &bull; Generated: $(_html_escape(ts_str)) &bull; EegFun v$(_html_escape(ver_str))</div>
</header>
""",
    )

    # KPI summary cards
    kept_str = isnan(mean_kept) ? "—" : "$(mean_kept)%"
    write(
        io,
        """
<div class="kpi-grid">
  <div class="kpi"><div class="label">Participants</div><div class="value">$n_subj</div><div class="sub">preprocessed datasets</div></div>
  <div class="kpi"><div class="label">Mean Epochs Retained</div><div class="value">$kept_str</div><div class="sub">threshold: $(round(low_epoch_threshold; digits = 1))%</div></div>
  <div class="kpi"><div class="label">Flagged for Review</div><div class="value">$n_flagged / $n_subj</div><div class="sub">$(n_flagged == 0 ? "all within thresholds" : "expanded below")</div></div>
  <div class="kpi"><div class="label">Continuous Repairs</div><div class="value">$mean_rep</div><div class="sub">avg channels / participant ($(nrow(electrode_freq_df)) unique)</div></div>
  <div class="kpi"><div class="label">ICA Removed</div><div class="value">$mean_ics</div><div class="sub">avg components / participant</div></div>
</div>
""",
    )

    _write_config_and_methods_section!(io, output_dir, pipeline_cfg, records, ver_str, low_epoch_threshold, max_ica_pct)

    # Cohort traffic-light overview table + CSV export + Flagged filter
    flagged_btn_html =
        n_flagged > 0 ?
        "<button type=\"button\" class=\"btn-ctrl flag-filter-btn\" onclick=\"toggleFlaggedOnly()\">Flagged Only ($n_flagged)</button>" : ""

    write(
        io,
        """
<section id="cohort-summary" class="panel">
  <div class="section-head">
    <h2>Cohort Quality-Control Summary</h2>
    <div class="btn-group">
      $flagged_btn_html
      <a class="btn-ctrl" href="$csv_uri" download="pipeline_qc_summary.csv">Download CSV</a>
    </div>
  </div>
  <table>
    <thead>
      <tr>
        <th>Participant</th>
        <th>Status</th>
        <th>Repaired Channels</th>
        <th>Skipped Channels</th>
        <th>Removed ICs</th>
        <th class="num">Total Epochs</th>
        <th class="num">Mean Kept (%)</th>
        <th>Condition Breakdown</th>
      </tr>
    </thead>
    <tbody>
""",
    )

    for (r_idx, rec) in enumerate(records)
        slug = slugs[r_idx]
        row_cls = rec.is_flagged ? " class=\"flagged-row\"" : " class=\"ok-row\""
        badge_html = if rec.is_flagged
            has_zero = !isempty(rec.condition_df) && any(==(0), rec.condition_df.n_kept)
            bcls = (!isempty(rec.skipped_channels) || has_zero) ? "badge-bad" : "badge-warn"
            "<span class=\"badge $bcls\" title=\"$(_html_escape(join(rec.flag_reasons, "; ")))\">Review</span>"
        else
            "<span class=\"badge badge-ok\">OK</span>"
        end
        rep_str = isempty(rec.repaired_channels) ? "—" : join(rec.repaired_channels, ", ")
        skip_html =
            isempty(rec.skipped_channels) ? "—" : "<span class=\"badge badge-bad\">$(_html_escape(join(rec.skipped_channels, ", ")))</span>"
        ic_str = _format_ic_breakdown(rec.ic_counts)
        trials_cell = rec.total_original > 0 ? "$(rec.total_kept) / $(rec.total_original)" : "—"
        pct_cell = isnan(rec.overall_kept_pct) ? "—" : "$(rec.overall_kept_pct)%"
        cond_parts =
            ["$(_html_escape(row.condition_name)): $(row.n_kept)/$(row.n_original) ($(row.pct_kept)%)" for row in eachrow(rec.condition_df)]
        cond_cell = isempty(cond_parts) ? "—" : join(cond_parts, " &bull; ")

        write(
            io,
            """      <tr$row_cls>
        <td><a class="subj-link" href="#$slug" onclick="const d=document.getElementById('$slug'); if(d) d.open=true;">$(_html_escape(rec.stem))</a></td>
        <td>$badge_html</td>
        <td>$(_html_escape(rep_str))</td>
        <td>$skip_html</td>
        <td>$(_html_escape(ic_str))</td>
        <td class="num">$trials_cell</td>
        <td class="num">$pct_cell</td>
        <td>$cond_cell</td>
      </tr>
""",
        )
    end
    write(io, "    </tbody>\n  </table>\n")

    cohort_heatmap_uri = can_plot ? _render_cohort_sensor_heatmap_uri(records) : nothing
    if !isnothing(cohort_heatmap_uri)
        write(
            io,
            """  <div class="fig-stack">
    <div class="fig-box"><img src="$cohort_heatmap_uri" alt="Cohort sensor repair and interpolation heatmap"><div class="cap">Cohort-wide electrode repair &amp; epoch-interpolation hotspot maps (check recurring sensors for cap/pin wear or gel bridges).</div></div>
  </div>
""",
        )
    end

    if !isempty(electrode_freq_df)
        write(io, "  <h3>Most Frequently Repaired Electrodes (Continuous Level)</h3>\n  <table>\n")
        write(io, "    <thead><tr><th>Electrode</th><th class=\"num\">Participants Repaired</th></tr></thead>\n    <tbody>\n")
        for row in eachrow(electrode_freq_df)
            write(
                io,
                "      <tr><td>$(_html_escape(string(row.electrode)))</td><td class=\"num\">$(row.n_participants) / $n_subj</td></tr>\n",
            )
        end
        write(io, "    </tbody>\n  </table>\n")
    end
    write(io, "</section>\n")

    # Grand-average & spectral cohort plots
    cohort_psd_ref = can_plot ? _cohort_mean_psd(records) : nothing
    if can_plot
        all_erps = ErpData[]
        for rec in records
            !isnothing(rec.erps) && append!(all_erps, rec.erps)
        end
        ga_erps = if !isempty(all_erps)
            try
                length(records) > 1 ? grand_average(all_erps) : all_erps
            catch
                ErpData[]
            end
        else
            ErpData[]
        end

        ga_erp_uri =
            !isempty(ga_erps) ? _render_erp_gfp_uri(ga_erps; plot_title = "Grand-Average ERPs & GFP", fig_size = (1040, 420)) : nothing
        ga_topo_uri = !isempty(ga_erps) ? _render_grand_topo_uri(ga_erps) : nothing
        ga_raw_vs_clean_uri = _render_cohort_raw_vs_clean_uri(records)
        ga_psd_uri = _render_cohort_psd_uri(records)

        if !isnothing(ga_erp_uri) || !isnothing(ga_topo_uri) || !isnothing(ga_raw_vs_clean_uri) || !isnothing(ga_psd_uri)
            write(io, "<section class=\"panel\"><h2>Cohort Grand-Average &amp; Spectral Preview</h2><div class=\"fig-stack\">")
            if !isnothing(ga_erp_uri)
                write(
                    io,
                    "<div class=\"fig-box\"><img src=\"$ga_erp_uri\" alt=\"Grand-average ERPs and GFP\"><div class=\"cap\">Baseline-corrected grand-average butterfly ERPs and Global Field Power across conditions.</div></div>",
                )
            end
            if !isnothing(ga_topo_uri)
                write(
                    io,
                    "<div class=\"fig-box\"><img src=\"$ga_topo_uri\" alt=\"Grand-average scalp topographies\"><div class=\"cap\">Time-resolved grand-average post-stimulus scalp topographies (shared symmetric μV scale).</div></div>",
                )
            end
            if !isnothing(ga_raw_vs_clean_uri)
                write(
                    io,
                    "<div class=\"fig-box\"><img src=\"$ga_raw_vs_clean_uri\" alt=\"Cohort pre- vs post-cleaning ERPs and GFP\"><div class=\"cap\">Pre- vs. post-preprocessing cohort comparison (±1 SEM) showing frontal ocular artifact suppression and cleaned GFP.</div></div>",
                )
            end
            if !isnothing(ga_psd_uri)
                write(
                    io,
                    "<div class=\"fig-box\"><img src=\"$ga_psd_uri\" alt=\"Cohort cleaned power spectral density\"><div class=\"cap\">Frequency-domain QC: 1/f spectral slope, posterior ~10 Hz alpha peak, and 50 Hz line-noise check across participants.</div></div>",
                )
            end
            write(io, "</div></section>\n")
        end
    end

    # Per-participant collapsible QC cards
    write(
        io,
        """<section class="panel">
  <div class="section-head">
    <h2>Per-Participant Quality-Control Details</h2>
    <div class="btn-group">
      $flagged_btn_html
      <button type="button" class="btn-ctrl" onclick="document.querySelectorAll('details.subj-card').forEach(d => d.open = true)">Expand All</button>
      <button type="button" class="btn-ctrl" onclick="document.querySelectorAll('details.subj-card').forEach(d => d.open = false)">Collapse All</button>
    </div>
  </div>
""",
    )
    for (r_idx, rec) in enumerate(records)
        slug = slugs[r_idx]
        open_attr = rec.is_flagged ? " open" : ""
        card_cls = rec.is_flagged ? "subj-card flagged" : "subj-card"
        badge_html = rec.is_flagged ? "<span class=\"badge badge-warn\">Review</span>" : "<span class=\"badge badge-ok\">OK</span>"
        trials_pill = rec.total_original > 0 ? "<span>Epochs: <strong>$(rec.total_kept)/$(rec.total_original)</strong></span>" : ""
        pct_pill = isnan(rec.overall_kept_pct) ? "" : "<span>Kept: <strong>$(rec.overall_kept_pct)%</strong></span>"

        write(
            io,
            """<details id="$slug" class="$card_cls"$open_attr>
  <summary>
    <span>$(_html_escape(rec.stem)) $badge_html</span>
    <span class="summary-pills">
      <span>Repaired ch: <strong>$(length(rec.repaired_channels))</strong></span>
      <span>Removed ICs: <strong>$(rec.ic_counts.total)</strong></span>
      $trials_pill
      $pct_pill
      <a class="top-link" href="#cohort-summary" onclick="event.stopPropagation();">↑ Summary</a>
    </span>
  </summary>
  <div class="card-body">
""",
        )

        if rec.is_flagged
            write(io, "    <p><strong>Review notes:</strong> $(_html_escape(join(rec.flag_reasons, "; ")))</p>\n")
        end

        if !isempty(rec.condition_df)
            write(
                io,
                """    <table>
      <thead>
        <tr>
          <th class="num">Condition</th>
          <th>Name</th>
          <th class="num">Original Epochs</th>
          <th class="num">Repaired Epochs</th>
          <th class="num">Rejected Epochs</th>
          <th class="num">Final Kept</th>
          <th class="num">Kept (%)</th>
        </tr>
      </thead>
      <tbody>
""",
            )
            for row in eachrow(rec.condition_df)
                write(
                    io,
                    """        <tr>
          <td class="num">$(row.condition)</td>
          <td>$(_html_escape(row.condition_name))</td>
          <td class="num">$(row.n_original)</td>
          <td class="num">$(row.n_repaired)</td>
          <td class="num">$(row.n_rejected)</td>
          <td class="num">$(row.n_kept)</td>
          <td class="num">$(row.pct_kept)%</td>
        </tr>
""",
                )
            end
            write(io, "      </tbody>\n    </table>\n")
        end

        if can_plot
            sensor_uri = _render_sensor_status_uri(rec)
            erp_uri = isnothing(rec.erps) ? nothing : _render_erp_gfp_uri(rec.erps; plot_title = "$(rec.stem) — ERPs & GFP")
            psd_uri = _render_participant_psd_uri(rec, cohort_psd_ref)
            raw_clean_uri = _render_participant_raw_vs_clean_uri(rec)
            ica_uri = _render_ica_artifacts_uri(rec)

            if !isnothing(sensor_uri) || !isnothing(erp_uri) || !isnothing(psd_uri) || !isnothing(raw_clean_uri) || !isnothing(ica_uri)
                write(io, "    <div class=\"fig-stack\">\n")
                if !isnothing(sensor_uri) || !isnothing(erp_uri)
                    top_cls = (!isnothing(sensor_uri) && !isnothing(erp_uri)) ? "fig-row-top" : "fig-row-single"
                    write(io, "      <div class=\"$top_cls\">\n")
                    if !isnothing(sensor_uri)
                        write(
                            io,
                            "        <div class=\"fig-box\"><img src=\"$sensor_uri\" alt=\"Sensor status map for $(_html_escape(rec.stem))\"><div class=\"cap\">$(_html_escape(_sensor_status_caption(rec)))</div></div>\n",
                        )
                    end
                    if !isnothing(erp_uri)
                        write(
                            io,
                            "        <div class=\"fig-box\"><img src=\"$erp_uri\" alt=\"Condition ERPs and GFP for $(_html_escape(rec.stem))\"><div class=\"cap\">Baseline-corrected condition ERP traces and Global Field Power.</div></div>\n",
                        )
                    end
                    write(io, "      </div>\n")
                end
                if !isnothing(psd_uri) || !isnothing(raw_clean_uri)
                    mid_cls = (!isnothing(psd_uri) && !isnothing(raw_clean_uri)) ? "fig-row-top" : "fig-row-single"
                    write(io, "      <div class=\"$mid_cls\">\n")
                    if !isnothing(psd_uri)
                        write(
                            io,
                            "        <div class=\"fig-box\"><img src=\"$psd_uri\" alt=\"Power spectral density for $(_html_escape(rec.stem))\"><div class=\"cap\">Spectral density vs. cohort mean (dotted guides at 10 Hz alpha &amp; 50 Hz line).</div></div>\n",
                        )
                    end
                    if !isnothing(raw_clean_uri)
                        write(
                            io,
                            "        <div class=\"fig-box\"><img src=\"$raw_clean_uri\" alt=\"Pre- vs post-cleaning overlay for $(_html_escape(rec.stem))\"><div class=\"cap\">Pre-ICA raw vs. cleaned frontal ROI and mean GFP.</div></div>\n",
                        )
                    end
                    write(io, "      </div>\n")
                end
                if !isnothing(ica_uri)
                    write(
                        io,
                        "      <div class=\"fig-box\"><img src=\"$ica_uri\" alt=\"Removed ICA topographies for $(_html_escape(rec.stem))\"><div class=\"cap\">Subtracted ICA artifact components ($(_html_escape(_format_ic_breakdown(rec.ic_counts)))).</div></div>\n",
                    )
                end
                write(io, "    </div>\n")
            end
        end

        write(io, "  </div>\n</details>\n")
    end
    write(
        io,
        """</section>
</div>
<script>
function toggleFlaggedOnly() {
  const active = document.body.classList.toggle("flagged-only");
  document.querySelectorAll(".flag-filter-btn").forEach(btn => btn.classList.toggle("active", active));
}
function copyMethodsText(btn) {
  const el = document.getElementById("methods-boilerplate");
  if (!el || !navigator.clipboard) return;
  navigator.clipboard.writeText(el.innerText).then(() => {
    const prev = btn.textContent;
    btn.textContent = "Copied!";
    setTimeout(() => { btn.textContent = prev; }, 1600);
  });
}
(function() {
  function openHashTarget() {
    if (!location.hash || location.hash.length < 2) return;
    const id = decodeURIComponent(location.hash.slice(1));
    const el = document.getElementById(id);
    if (el && el.tagName === "DETAILS") {
      document.body.classList.remove("flagged-only");
      document.querySelectorAll(".flag-filter-btn").forEach(btn => btn.classList.remove("active"));
      el.open = true;
      el.scrollIntoView({ behavior: "smooth", block: "start" });
    }
  }
  window.addEventListener("DOMContentLoaded", openHashTarget);
  window.addEventListener("hashchange", openHashTarget);
})();
</script>
</body>
</html>
""",
    )
    return String(take!(io))
end
