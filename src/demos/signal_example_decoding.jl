"""
    signal_example_decoding()

Interactive MVPA Decoding Playground — Teaching Time-Resolved Classification.

Demonstrates how a Support Vector Machine (SVM) classifier can learn to
distinguish two experimental conditions from multi-channel EEG data, one
time point at a time. Synthetic data is generated internally — no real
data required.

## Plots

| Row | Plot | Description |
|-----|------|-------------|
| 1 | ERP Waveforms | Mean waveform per condition (blue vs red) with ±1 SE shading |
| 2 | Decoding Accuracy | Time-resolved classification accuracy with 50% chance line and SE band |

## Controls

| Control | Range | Description |
|---------|-------|-------------|
| Signal Strength | 0–3 | Amplitude of the condition difference (larger = easier to decode) |
| Noise Level | 0.1–5 | Background noise amplitude |
| N Trials | 10–200 | Trials per condition — more trials = better estimates |
| N Channels | 1–10 | Number of EEG channels (features for the classifier) |
| Effect Onset | −0.2–0.8 s | When the condition difference begins |
| Effect Duration | 0.05–0.8 s | How long the condition difference lasts |
| [Decode!] | button | Runs the SVM classification (takes a few seconds) |

## Teaching Notes

**MVPA (Multi-Variate Pattern Analysis)** asks: at each time point, can a
classifier distinguish condition A from condition B using the pattern of
activity across channels?

Unlike a t-test (which tests one channel at a time), MVPA considers the
*joint* pattern across all channels — exploiting correlations that univariate
methods miss.

**Cross-validation** prevents overfitting. The data is split into training
and test sets; the classifier only sees training data during learning. If it
performs above chance on the held-out test data, the effect is genuine.

**Key insights to discover:**

- *Signal strength* directly controls accuracy — this is the effect size.
- *More channels* help: the "multi" in multivariate adds information.
- *More trials* produce smoother, more reliable accuracy estimates.
- *Effect onset/duration* maps directly to the accuracy peak — the classifier
  finds the information exactly where it exists.
- *Zero signal strength* → accuracy fluctuates around 50% (chance). This is
  what a null result looks like.

## See Also

- Grootswagers, T., Wardle, S. G., & Carlson, T. A. (2017). Decoding dynamic brain
  patterns from evoked responses: a tutorial on multivariate pattern analysis
  applied to time series neuroimaging data. *Journal of Cognitive Neuroscience*,
  29(4), 677–697.
- King, J.R. & Dehaene, S. (2014). Characterizing the dynamics of mental
  representations: the temporal generalization method. *Trends in Cognitive
  Sciences*, 18(4), 203-210.

# Examples
```julia
using EegFun
EegFun.signal_example_decoding()
```

# Returns
- `fig::Figure`: The Makie figure object
- `axes::Tuple`: `(ax_erp, ax_acc)` — the two axis objects
"""
function signal_example_decoding()

    FS     = 250.0    # sample rate (Hz) — typical EEG
    T_PRE  = -0.2     # pre-stimulus (s)
    T_POST = 1.0      # post-stimulus (s)
    t_arr  = collect(range(T_PRE, T_POST, step = 1.0 / FS))
    n_time = length(t_arr)

    fig = Figure(size = (1400, 750), figure_padding = (8, 8, 8, 8), title = "MVPA Decoding Playground")

    # ── Font sizing ──────────────────────────────────────────────────────────
    lbl_sz   = Observable(17)
    tick_sz  = Observable(14)
    title_sz = Observable(19)
    ctrl_sz  = Observable(17)

    on(fig.scene.viewport) do area
        sf         = area.widths[1] / 3000
        lbl_sz[]   = max(12, round(Int, 17 * sf))
        tick_sz[]  = max(10, round(Int, 14 * sf))
        title_sz[] = max(14, round(Int, 19 * sf))
        ctrl_sz[]  = max(12, round(Int, 17 * sf))
    end

    # ── Axes ─────────────────────────────────────────────────────────────────
    ax_erp = Axis(
        fig[1, 1],
        title          = "ERP Waveforms (channel average)",
        ylabel         = "Amplitude (μV)",
        titlesize      = title_sz,
        xlabelsize     = lbl_sz,
        ylabelsize     = lbl_sz,
        xticklabelsize = tick_sz,
        yticklabelsize = tick_sz,
        xticklabelsvisible = false,
        xgridvisible   = true,
        ygridvisible   = true,
    )

    ax_acc = Axis(
        fig[2, 1],
        title          = "Decoding Accuracy (SVM, k-fold CV)",
        xlabel         = "Time (s)",
        ylabel         = "Accuracy",
        titlesize      = title_sz,
        xlabelsize     = lbl_sz,
        ylabelsize     = lbl_sz,
        xticklabelsize = tick_sz,
        yticklabelsize = tick_sz,
        xgridvisible   = true,
        ygridvisible   = true,
    )

    # ── Parameter observables ────────────────────────────────────────────────
    sig_strength   = Observable(1.0)
    noise_level    = Observable(1.0)
    n_trials       = Observable(50)
    n_channels     = Observable(3)
    effect_onset   = Observable(0.1)
    effect_duration = Observable(0.4)

    # ── Data observables ─────────────────────────────────────────────────────
    erp_a_mean = Observable(zeros(n_time))
    erp_b_mean = Observable(zeros(n_time))
    erp_a_lo   = Observable(zeros(n_time))
    erp_a_hi   = Observable(zeros(n_time))
    erp_b_lo   = Observable(zeros(n_time))
    erp_b_hi   = Observable(zeros(n_time))

    acc_mean   = Observable(fill(0.5, n_time))
    acc_lo     = Observable(fill(0.5, n_time))
    acc_hi     = Observable(fill(0.5, n_time))

    status_text = Observable("Press [Decode!] to classify")

    # ── Generate synthetic data ──────────────────────────────────────────────
    function generate_data()
        n_ch  = n_channels[]
        n_tr  = n_trials[]
        sig   = sig_strength[]
        noise = noise_level[]
        onset = effect_onset[]
        dur   = effect_duration[]

        # Gaussian bump centred in the effect window
        centre = onset + dur / 2.0
        sigma  = dur / 4.0  # ~95% of bump within the window

        bump = sig .* exp.(-((t_arr .- centre) .^ 2) ./ (2 * sigma^2))

        # Generate trials: [channels × time × trials]
        data_a = noise .* randn(n_ch, n_time, n_tr)
        data_b = noise .* randn(n_ch, n_time, n_tr)

        # Add signal: condition A = +bump, condition B = -bump (across all channels)
        for tr = 1:n_tr
            for ch = 1:n_ch
                # Each channel gets slightly different bump amplitude, centred around 1.0
                ch_scale = 0.7 + 0.6 * ((ch - 0.5) / n_ch)  # range ~0.7–1.3, symmetric
                @views data_a[ch, :, tr] .+= bump .* ch_scale
                @views data_b[ch, :, tr] .-= bump .* ch_scale
            end
        end

        # Compute ERPs (channel-averaged across channels, then mean ± SE across trials)
        erp_a_trials_2d = dropdims(mean(data_a, dims=1), dims=1)  # [time × trials]
        erp_b_trials_2d = dropdims(mean(data_b, dims=1), dims=1)  # [time × trials]

        a_m = vec(mean(erp_a_trials_2d, dims=2))
        a_se = vec(std(erp_a_trials_2d, dims=2) ./ sqrt(n_tr))
        b_m = vec(mean(erp_b_trials_2d, dims=2))
        b_se = vec(std(erp_b_trials_2d, dims=2) ./ sqrt(n_tr))

        erp_a_mean[] = a_m
        erp_a_lo[]   = a_m .- a_se
        erp_a_hi[]   = a_m .+ a_se
        erp_b_mean[] = b_m
        erp_b_lo[]   = b_m .- b_se
        erp_b_hi[]   = b_m .+ b_se

        # Auto-scale ERP axis
        all_vals = vcat(erp_a_lo[], erp_a_hi[], erp_b_lo[], erp_b_hi[])
        y_range = maximum(abs.(all_vals)) * 1.2
        ylims!(ax_erp, -y_range, y_range)

        return data_a, data_b
    end

    # ── Lightweight decoding ─────────────────────────────────────────────────
    is_decoding = Ref(false)  # guard against overlapping runs

    function run_decoding(data_a, data_b)
        n_ch   = size(data_a, 1)
        n_tr_a = size(data_a, 3)
        n_tr_b = size(data_b, 3)
        n_tr   = min(n_tr_a, n_tr_b)  # equalize

        n_iterations = 10
        n_folds      = 3
        n_per_fold   = n_tr ÷ n_folds

        # Pre-compute CV splits once (avoids per-timepoint allocation)
        total_trials = 2 * n_tr
        cv_splits = Vector{Tuple{Vector{Int}, Vector{Int}}}(undef, n_folds)
        for fold = 1:n_folds
            test_start_a = (fold - 1) * n_per_fold + 1
            test_end_a   = fold * n_per_fold
            test_start_b = n_tr + test_start_a
            test_end_b   = n_tr + test_end_a
            test_idx  = vcat(test_start_a:test_end_a, test_start_b:test_end_b)
            train_idx = setdiff(1:total_trials, test_idx)
            cv_splits[fold] = (train_idx, test_idx)
        end

        # Preallocate
        accuracies = zeros(n_iterations, n_time)
        x_all  = Matrix{Float64}(undef, total_trials, n_ch)
        labels = Vector{Int}(undef, total_trials)

        for iter = 1:n_iterations
            # Shuffle trials
            perm_a = shuffle(1:n_tr)
            perm_b = shuffle(1:n_tr)

            for t = 1:n_time
                # Extract features at this time point
                for i = 1:n_tr
                    @inbounds for ch = 1:n_ch
                        x_all[i, ch]      = data_a[ch, t, perm_a[i]]
                        x_all[n_tr+i, ch] = data_b[ch, t, perm_b[i]]
                    end
                    labels[i]      = 1
                    labels[n_tr+i] = 2
                end

                # k-fold cross-validation
                fold_accs = zeros(n_folds)
                for fold = 1:n_folds
                    train_idx, test_idx = cv_splits[fold]

                    X_train = @view x_all[train_idx, :]
                    y_train = labels[train_idx]
                    X_test  = @view x_all[test_idx, :]
                    y_test  = @view labels[test_idx]

                    y_pred = libsvm_classifier(X_train, y_train, X_test)

                    n_correct = 0
                    @inbounds for i in eachindex(y_test)
                        n_correct += (y_pred[i] == y_test[i])
                    end
                    fold_accs[fold] = n_correct / length(y_test)
                end

                accuracies[iter, t] = mean(fold_accs)
            end
        end

        # Mean ± SE across iterations
        a_m  = vec(mean(accuracies, dims=1))
        a_se = vec(std(accuracies, dims=1) ./ sqrt(n_iterations))

        acc_mean[] = a_m
        acc_lo[]   = a_m .- a_se
        acc_hi[]   = a_m .+ a_se

        # Auto-scale accuracy axis
        y_lo = max(0.3, minimum(acc_lo[]) - 0.05)
        y_hi = min(1.0, maximum(acc_hi[]) + 0.05)
        ylims!(ax_acc, y_lo, y_hi)
    end

    # ── Initial data + ERP plot ──────────────────────────────────────────────
    current_data = Ref{Tuple{Array{Float64,3}, Array{Float64,3}}}(generate_data())

    # ERP plots
    band!(ax_erp, t_arr, erp_a_lo, erp_a_hi, color = (:royalblue, 0.2))
    lines!(ax_erp, t_arr, erp_a_mean, color = :royalblue, linewidth = 2, label = "Condition A")
    band!(ax_erp, t_arr, erp_b_lo, erp_b_hi, color = (:crimson, 0.2))
    lines!(ax_erp, t_arr, erp_b_mean, color = :crimson, linewidth = 2, label = "Condition B")
    vlines!(ax_erp, [0.0], color = :black, linewidth = 1, linestyle = :dash)
    hlines!(ax_erp, [0.0], color = :black, linewidth = 0.5)
    axislegend(ax_erp, position = :rt)

    # Accuracy plot
    hlines!(ax_acc, [0.5], color = :grey50, linewidth = 1.5, linestyle = :dash, label = "Chance (50%)")
    band!(ax_acc, t_arr, acc_lo, acc_hi, color = (:black, 0.15))
    lines!(ax_acc, t_arr, acc_mean, color = :black, linewidth = 2.5, label = "Accuracy")
    vlines!(ax_acc, [0.0], color = :black, linewidth = 1, linestyle = :dash)
    axislegend(ax_acc, position = :rt)

    # Status text (top-left of accuracy panel)
    text!(
        ax_acc,
        0.02,
        0.97,
        text = status_text,
        space = :relative,
        fontsize = @lift(max(14, round(Int, $lbl_sz * 1.1))),
        color = :grey40,
        align = (:left, :top),
    )

    xlims!(ax_erp, T_PRE, T_POST)
    xlims!(ax_acc, T_PRE, T_POST)
    ylims!(ax_acc, 0.35, 0.85)

    # ── Controls ─────────────────────────────────────────────────────────────
    ctrl = GridLayout(fig[3, 1], colgap = 16)

    function labelled_slider(parent, col, header, range_vals, startval, fmt)
        Label(parent[1, col], header, fontsize = ctrl_sz, halign = :center)
        sl  = Slider(parent[2, col], range = range_vals, startvalue = startval)
        lbl = Label(parent[3, col], fmt(startval), fontsize = ctrl_sz, halign = :center)
        return sl, lbl
    end

    sl_ss, lbl_ss = labelled_slider(ctrl, 1, "Signal Strength", 0.0:0.1:3.0, 1.0, v -> "$(round(v, digits=1))")
    sl_nl, lbl_nl = labelled_slider(ctrl, 2, "Noise Level", 0.1:0.1:5.0, 1.0, v -> "$(round(v, digits=1))")
    sl_nt, lbl_nt = labelled_slider(ctrl, 3, "N Trials", 10:10:200, 50, v -> "$(v)")
    sl_nc, lbl_nc = labelled_slider(ctrl, 4, "N Channels", 1:1:10, 3, v -> "$(v)")
    sl_eo, lbl_eo = labelled_slider(ctrl, 5, "Effect Onset", -0.2:0.05:0.8, 0.1, v -> "$(round(v, digits=2)) s")
    sl_ed, lbl_ed = labelled_slider(ctrl, 6, "Effect Duration", 0.05:0.05:0.8, 0.4, v -> "$(round(v, digits=2)) s")

    # Decode button
    btn_decode = Button(ctrl[2, 7], label = "Decode!", fontsize = ctrl_sz)
    Label(ctrl[1, 7], "", fontsize = ctrl_sz, halign = :center)
    Label(ctrl[3, 7], "", fontsize = ctrl_sz, halign = :center)

    # Wire sliders → update ERP immediately (fast)
    on(sl_ss.value) do v
        sig_strength[] = v
        lbl_ss.text = "$(round(v, digits=1))"
        current_data[] = generate_data()
        status_text[] = "Parameters changed — press [Decode!]"
    end
    on(sl_nl.value) do v
        noise_level[] = v
        lbl_nl.text = "$(round(v, digits=1))"
        current_data[] = generate_data()
        status_text[] = "Parameters changed — press [Decode!]"
    end
    on(sl_nt.value) do v
        n_trials[] = v
        lbl_nt.text = "$(v)"
        current_data[] = generate_data()
        status_text[] = "Parameters changed — press [Decode!]"
    end
    on(sl_nc.value) do v
        n_channels[] = v
        lbl_nc.text = "$(v)"
        current_data[] = generate_data()
        status_text[] = "Parameters changed — press [Decode!]"
    end
    on(sl_eo.value) do v
        effect_onset[] = v
        lbl_eo.text = "$(round(v, digits=2)) s"
        current_data[] = generate_data()
        status_text[] = "Parameters changed — press [Decode!]"
    end
    on(sl_ed.value) do v
        effect_duration[] = v
        lbl_ed.text = "$(round(v, digits=2)) s"
        current_data[] = generate_data()
        status_text[] = "Parameters changed — press [Decode!]"
    end

    # Button → run decoding (slow, with guard against overlapping runs)
    on(btn_decode.clicks) do _
        is_decoding[] && return  # ignore if already running
        is_decoding[] = true
        status_text[] = "Decoding... (this takes a few seconds)"
        @async begin
            try
                data_a, data_b = current_data[]
                run_decoding(data_a, data_b)
                status_text[] = "Done! Adjust parameters and decode again."
            catch e
                status_text[] = "Error: $(sprint(showerror, e))"
            finally
                is_decoding[] = false
            end
        end
    end

    # Force control columns to fill the full figure width
    for col = 1:7
        colsize!(ctrl, col, Relative(1 / 7))
    end

    # ── Row sizing ───────────────────────────────────────────────────────────
    rowsize!(fig.layout, 1, Relative(0.35))
    rowsize!(fig.layout, 2, Relative(0.45))
    rowsize!(fig.layout, 3, Relative(0.20))

    display(fig)
    return fig, (ax_erp, ax_acc)
end
