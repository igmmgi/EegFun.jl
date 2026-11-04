"""
    signal_example_decoding()

Interactive MVPA decoding playground demonstrating time-resolved k-fold cross-validated SVM
classification of two synthetic multi-channel EEG conditions.

# Examples
```julia
using GLMakie, EegFun
EegFun.signal_example_decoding()
```

# Returns
- `fig::Figure`: The Makie figure object
- `axes::Tuple`: `(ax_erp, ax_acc)` — the two axis objects
"""
function signal_example_decoding()
    _check_interactive_demo_backend("MVPA Decoding demo") || return nothing

    FS = 250.0
    T_PRE = -0.2
    T_POST = 1.0
    t_arr = collect(range(T_PRE, T_POST, step = 1.0 / FS))
    n_time = length(t_arr)

    fig = _eegfun_figure(size = (1400, 750), figure_padding = (8, 8, 8, 8), title = "MVPA Decoding Playground")
    lbl_sz, tick_sz, title_sz, ctrl_sz = _demo_adaptive_fonts(fig, 1400, (15, 11), (12, 10), (18, 13), (15, 11))

    ax_erp = Axis(
        fig[1, 1],
        title = "ERP Waveforms (channel average)",
        ylabel = "Amplitude (μV)",
        titlesize = title_sz,
        xlabelsize = lbl_sz,
        ylabelsize = lbl_sz,
        xticklabelsize = tick_sz,
        yticklabelsize = tick_sz,
        xticklabelsvisible = false,
        xgridvisible = true,
        ygridvisible = true,
    )
    ax_acc = Axis(
        fig[2, 1],
        title = "Decoding Accuracy (SVM, k-fold CV)",
        xlabel = "Time (s)",
        ylabel = "Accuracy",
        titlesize = title_sz,
        xlabelsize = lbl_sz,
        ylabelsize = lbl_sz,
        xticklabelsize = tick_sz,
        yticklabelsize = tick_sz,
        xgridvisible = true,
        ygridvisible = true,
    )

    sig_strength = Observable(1.0)
    noise_level = Observable(1.0)
    n_trials = Observable(50)
    n_channels = Observable(3)
    effect_onset = Observable(0.1)
    effect_duration = Observable(0.4)

    erp_a_mean = Observable(zeros(n_time))
    erp_b_mean = Observable(zeros(n_time))
    erp_a_lo = Observable(zeros(n_time))
    erp_a_hi = Observable(zeros(n_time))
    erp_b_lo = Observable(zeros(n_time))
    erp_b_hi = Observable(zeros(n_time))

    acc_mean = Observable(fill(0.5, n_time))
    acc_lo = Observable(fill(0.5, n_time))
    acc_hi = Observable(fill(0.5, n_time))
    status_text = Observable("Press [Decode!] to classify")

    function generate_data()
        n_ch = n_channels[]
        n_tr = n_trials[]
        sig = sig_strength[]
        noise = noise_level[]
        onset = effect_onset[]
        dur = effect_duration[]

        centre = onset + dur / 2.0
        sigma = dur / 4.0
        bump = @. sig * exp(-((t_arr - centre)^2) / (2.0 * sigma^2))

        data_a = noise .* randn(n_ch, n_time, n_tr)
        data_b = noise .* randn(n_ch, n_time, n_tr)
        ch_scales = [0.7 + 0.6 * ((ch - 0.5) / n_ch) for ch = 1:n_ch]

        @inbounds for tr = 1:n_tr, t = 1:n_time
            bval = bump[t]
            for ch = 1:n_ch
                delta = bval * ch_scales[ch]
                data_a[ch, t, tr] += delta
                data_b[ch, t, tr] -= delta
            end
        end

        inv_ch = 1.0 / n_ch
        inv_sqrt_tr = 1.0 / sqrt(n_tr)
        a_m = Vector{Float64}(undef, n_time)
        b_m = Vector{Float64}(undef, n_time)
        a_lo = Vector{Float64}(undef, n_time)
        a_hi = Vector{Float64}(undef, n_time)
        b_lo = Vector{Float64}(undef, n_time)
        b_hi = Vector{Float64}(undef, n_time)
        tr_buf_a = Vector{Float64}(undef, n_tr)
        tr_buf_b = Vector{Float64}(undef, n_tr)

        @inbounds for t = 1:n_time
            for tr = 1:n_tr
                sa = 0.0
                sb = 0.0
                for ch = 1:n_ch
                    sa += data_a[ch, t, tr]
                    sb += data_b[ch, t, tr]
                end
                tr_buf_a[tr] = sa * inv_ch
                tr_buf_b[tr] = sb * inv_ch
            end
            ma, se_a = mean(tr_buf_a), std(tr_buf_a) * inv_sqrt_tr
            mb, se_b = mean(tr_buf_b), std(tr_buf_b) * inv_sqrt_tr
            a_m[t], a_lo[t], a_hi[t] = ma, ma - se_a, ma + se_a
            b_m[t], b_lo[t], b_hi[t] = mb, mb - se_b, mb + se_b
        end

        erp_a_mean[] = a_m
        erp_a_lo[] = a_lo
        erp_a_hi[] = a_hi
        erp_b_mean[] = b_m
        erp_b_lo[] = b_lo
        erp_b_hi[] = b_hi

        y_range = max(maximum(abs, a_lo), maximum(abs, a_hi), maximum(abs, b_lo), maximum(abs, b_hi), 1e-3) * 1.2
        ylims!(ax_erp, -y_range, y_range)

        return data_a, data_b
    end

    is_decoding = Ref(false)

    function run_decoding(data_a::Array{Float64,3}, data_b::Array{Float64,3})
        n_ch = size(data_a, 1)
        n_tr = min(size(data_a, 3), size(data_b, 3))
        n_iterations = 10
        n_folds = 3
        n_per_fold = n_tr ÷ n_folds
        total_trials = 2 * n_tr

        labels = Vector{Int}(undef, total_trials)
        fill!(@view(labels[1:n_tr]), 1)
        fill!(@view(labels[(n_tr+1):total_trials]), 2)

        cv_splits = Vector{Tuple{Vector{Int},Vector{Int},Vector{Int}}}(undef, n_folds)
        for fold = 1:n_folds
            test_start_a = (fold - 1) * n_per_fold + 1
            test_end_a = fold * n_per_fold
            test_idx = vcat(test_start_a:test_end_a, (n_tr+test_start_a):(n_tr+test_end_a))
            train_idx = setdiff(1:total_trials, test_idx)
            cv_splits[fold] = (train_idx, test_idx, labels[train_idx])
        end

        accuracies = Matrix{Float64}(undef, n_iterations, n_time)
        x_all = Matrix{Float64}(undef, total_trials, n_ch)

        for iter = 1:n_iterations
            perm_a = randperm(n_tr)
            perm_b = randperm(n_tr)

            for t = 1:n_time
                @inbounds for ch = 1:n_ch, i = 1:n_tr
                    x_all[i, ch] = data_a[ch, t, perm_a[i]]
                    x_all[n_tr+i, ch] = data_b[ch, t, perm_b[i]]
                end

                acc_sum = 0.0
                for fold = 1:n_folds
                    train_idx, test_idx, y_train = cv_splits[fold]
                    X_train = @view x_all[train_idx, :]
                    X_test = @view x_all[test_idx, :]
                    y_test = @view labels[test_idx]
                    y_pred = libsvm_classifier(X_train, y_train, X_test)

                    n_correct = 0
                    @inbounds for i in eachindex(y_test)
                        n_correct += (y_pred[i] == y_test[i])
                    end
                    acc_sum += n_correct / length(y_test)
                end
                accuracies[iter, t] = acc_sum / n_folds
            end
        end

        a_m = vec(mean(accuracies, dims = 1))
        a_se = vec(std(accuracies, dims = 1) ./ sqrt(n_iterations))

        acc_mean[] = a_m
        acc_lo[] = a_m .- a_se
        acc_hi[] = a_m .+ a_se

        ylims!(ax_acc, max(0.3, minimum(acc_lo[]) - 0.05), min(1.0, maximum(acc_hi[]) + 0.05))
    end

    current_data = Ref{Tuple{Array{Float64,3},Array{Float64,3}}}(generate_data())

    band!(ax_erp, t_arr, erp_a_lo, erp_a_hi, color = (:royalblue, 0.2))
    lines!(ax_erp, t_arr, erp_a_mean, color = :royalblue, linewidth = 2, label = "Condition A")
    band!(ax_erp, t_arr, erp_b_lo, erp_b_hi, color = (:crimson, 0.2))
    lines!(ax_erp, t_arr, erp_b_mean, color = :crimson, linewidth = 2, label = "Condition B")
    vlines!(ax_erp, [0.0], color = :black, linewidth = 1, linestyle = :dash)
    hlines!(ax_erp, [0.0], color = :black, linewidth = 0.5)
    axislegend(ax_erp, position = :rt)

    hlines!(ax_acc, [0.5], color = :grey50, linewidth = 1.5, linestyle = :dash, label = "Chance (50%)")
    band!(ax_acc, t_arr, acc_lo, acc_hi, color = (:black, 0.15))
    lines!(ax_acc, t_arr, acc_mean, color = :black, linewidth = 2.5, label = "Accuracy")
    vlines!(ax_acc, [0.0], color = :black, linewidth = 1, linestyle = :dash)
    axislegend(ax_acc, position = :rt)

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

    ctrl = GridLayout(fig[3, 1], colgap = 16)
    sl_ss, lbl_ss = _demo_labelled_slider(ctrl, 1, "Signal Strength", 0.0:0.1:3.0, 1.0, v -> "$(round(v, digits=1))", ctrl_sz)
    sl_nl, lbl_nl = _demo_labelled_slider(ctrl, 2, "Noise Level", 0.1:0.1:5.0, 1.0, v -> "$(round(v, digits=1))", ctrl_sz)
    sl_nt, lbl_nt = _demo_labelled_slider(ctrl, 3, "N Trials", 10:10:200, 50, string, ctrl_sz)
    sl_nc, lbl_nc = _demo_labelled_slider(ctrl, 4, "N Channels", 1:1:10, 3, string, ctrl_sz)
    sl_eo, lbl_eo = _demo_labelled_slider(ctrl, 5, "Effect Onset", -0.2:0.05:0.8, 0.1, v -> "$(round(v, digits=2)) s", ctrl_sz)
    sl_ed, lbl_ed = _demo_labelled_slider(ctrl, 6, "Effect Duration", 0.05:0.05:0.8, 0.4, v -> "$(round(v, digits=2)) s", ctrl_sz)

    btn_decode = Button(ctrl[2, 7], label = "Decode!", fontsize = ctrl_sz)
    Label(ctrl[1, 7], "", fontsize = ctrl_sz, halign = :center)
    Label(ctrl[3, 7], "", fontsize = ctrl_sz, halign = :center)

    function _bind_param!(sl, obs, lbl, fmt)
        on(sl.value) do v
            obs[] = v
            lbl.text = fmt(v)
            current_data[] = generate_data()
            status_text[] = "Parameters changed — press [Decode!]"
        end
    end

    _bind_param!(sl_ss, sig_strength, lbl_ss, v -> "$(round(v, digits=1))")
    _bind_param!(sl_nl, noise_level, lbl_nl, v -> "$(round(v, digits=1))")
    _bind_param!(sl_nt, n_trials, lbl_nt, string)
    _bind_param!(sl_nc, n_channels, lbl_nc, string)
    _bind_param!(sl_eo, effect_onset, lbl_eo, v -> "$(round(v, digits=2)) s")
    _bind_param!(sl_ed, effect_duration, lbl_ed, v -> "$(round(v, digits=2)) s")

    on(btn_decode.clicks) do _
        is_decoding[] && return
        is_decoding[] = true
        status_text[] = "Decoding... (this takes a few seconds)"
        @async begin
            try
                events(fig).window_open[] || return
                data_a, data_b = current_data[]
                run_decoding(data_a, data_b)
                events(fig).window_open[] && (status_text[] = "Done! Adjust parameters and decode again.")
            catch e
                events(fig).window_open[] && (status_text[] = "Error: $(sprint(showerror, e))")
            finally
                is_decoding[] = false
            end
        end
    end

    for col = 1:7
        colsize!(ctrl, col, Relative(1 / 7))
    end

    rowsize!(fig.layout, 1, Relative(0.35))
    rowsize!(fig.layout, 2, Relative(0.45))
    rowsize!(fig.layout, 3, Relative(0.20))

    _display_demo_figure(fig)
    return fig, (ax_erp, ax_acc)
end
