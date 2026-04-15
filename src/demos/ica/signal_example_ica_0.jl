"""
    signal_example_ica_0()

Interactive Demo — Part 0: Matrix Math Basics (Row-by-Column).

A foundational sandbox demonstrating exactly *how* matrix multiplication
(U = W * X) operates on time-series data using a 3-channel setup.

## Educational Value
By freezing time and letting the user manually slide a time cursor across
the data, this demo exposes the hidden "row-by-column" dot product mechanism.
Users literally watch the matrix sweep across the data columns step-by-step
to draw the new component waveforms.

## See Also
- [`signal_example_ica_1`](@ref) — Part 1: What is a Mixture?

# Examples
```julia
using EegFun
EegFun.signal_example_ica_0()
```
"""
function signal_example_ica_0()
    N = 800
    Random.seed!(42)

    fig = Figure(size = (1500, 1000), figure_padding = (15, 15, 15, 15))

    lbl_sz   = Observable(16)
    title_sz = Observable(20)
    ctrl_sz  = Observable(18)
    math_sz  = Observable(22)

    # ── Source Generation ─────────────────────────────────────────────────────
    t = range(0, 4π, length = N)

    s1 = sin.(5 .* t) .+ 0.2 .* randn(N)                 # High freq
    s2 = 1.5 .* sin.(1.2 .* t) .+ 0.2 .* randn(N)        # Slow wave
    s3 = 0.8 .* sin.(15 .* t) .+ 1.5 .* exp.(-(t .- 6) .^ 2) # Blink-like transient + noise

    M = [
        1.0  0.5 -0.2
        0.5 -1.0  0.6
        0.1  0.8  1.0
    ]

    X_raw = M * vcat(s1', s2', s3')

    eeg1 = X_raw[1, :]
    eeg2 = X_raw[2, :]
    eeg3 = X_raw[3, :]

    # ── Observables ───────────────────────────────────────────────────────────
    time_idx = Observable(1)

    # 3x3 array for cleaner processing loop
    w_obs = [
        Observable(1.0) Observable(1.0) Observable(1.0)
        Observable(1.0) Observable(1.0) Observable(1.0)
        Observable(1.0) Observable(1.0) Observable(1.0)
    ]

    # Arrays for plotting the "drawn" components
    u1_full = Observable(zeros(N))
    u2_full = Observable(zeros(N))
    u3_full = Observable(zeros(N))

    # ── Layout: Matrix Math Board (Top) ───────────────────────────────────────
    Label(fig[1, 1], "Matrix Math (U = WX)", fontsize = @lift($title_sz + 8), font = :bold, color = :black)

    math_grid = GridLayout(fig[2, 1], alignmode = Outside(10))
    # Titles
    Label(math_grid[1, 1], "Unmixing Matrix [W]", fontsize = ctrl_sz, font = :bold, color = :black)
    Label(math_grid[1, 2], "", fontsize = ctrl_sz) # Spacer
    Label(math_grid[1, 3], "Data Column [X_t]", fontsize = ctrl_sz, font = :bold, color = :black)
    Label(math_grid[1, 4], "", fontsize = ctrl_sz) # Spacer
    Label(math_grid[1, 5], "Row × Column Math", fontsize = ctrl_sz, font = :bold, color = :black)
    Label(math_grid[1, 6], "", fontsize = ctrl_sz) # Spacer
    Label(math_grid[1, 7], "Output [U_t]", fontsize = ctrl_sz, font = :bold, color = :black)

    # Concept: Colors map to the Input Columns. 
    # Col 1 = Blue, Col 2 = Red, Col 3 = Purple
    colors = [:blue, :red, :purple]

    # W Matrix Box
    w_box = GridLayout(math_grid[2, 1])
    Box(w_box[1:3, 1:3], color = (:transparent, 0.0), strokewidth = 3, strokecolor = :black)

    sl_grid = Matrix{Any}(undef, 3, 3)

    for r = 1:3
        for c = 1:3
            start_val = 1.0
            sl = Slider(w_box[r, c], range = -2.0:0.05:2.0, startvalue = start_val, width = 120)
            sl.color_active = (:transparent, 0.0)
            sl.color_inactive = (:transparent, 0.0)

            # Color relates to the Column (c) it multiplies!
            Label(
                w_box[r, c],
                @lift(Printf.@sprintf("%.2f", $(sl.value))),
                fontsize = math_sz,
                font = :bold,
                color = colors[c],
                padding = (15, 15, 10, 10),
            )

            sl_grid[r, c] = sl
            on(sl.value) do val
                w_obs[r, c][] = val
                update_calcs()
            end
        end
    end

    colgap!(w_box, 15)
    rowgap!(w_box, 10)

    Label(math_grid[2, 2], " × ", fontsize = @lift($math_sz + 10), font = :bold)

    # X Column Vector
    x_box = GridLayout(math_grid[2, 3])
    Box(x_box[1:3, 1], color = (:transparent, 0.0), strokewidth = 3, strokecolor = :black)
    lbl_x1 = Label(
        x_box[1, 1],
        @lift(Printf.@sprintf("%.2f", eeg1[$time_idx])),
        fontsize = math_sz,
        font = :bold,
        color = colors[1],
        padding = (10, 10, 10, 10),
    )
    lbl_x2 = Label(
        x_box[2, 1],
        @lift(Printf.@sprintf("%.2f", eeg2[$time_idx])),
        fontsize = math_sz,
        font = :bold,
        color = colors[2],
        padding = (10, 10, 10, 10),
    )
    lbl_x3 = Label(
        x_box[3, 1],
        @lift(Printf.@sprintf("%.2f", eeg3[$time_idx])),
        fontsize = math_sz,
        font = :bold,
        color = colors[3],
        padding = (10, 10, 10, 10),
    )

    Label(math_grid[2, 4], " = ", fontsize = @lift($math_sz + 10), font = :bold)

    # The Math Execution Box
    exec_box = GridLayout(math_grid[2, 5])
    Box(exec_box[1:3, 1], color = (:transparent, 0.0), strokewidth = 3, strokecolor = :black)

    # Create the math strings dynamically with explicit color blocks per component
    for r = 1:3
        row_grid = GridLayout(exec_box[r, 1], alignmode = Inside())

        Label(row_grid[1, 1], "(", fontsize = math_sz, font = :bold, color = :black)
        Label(
            row_grid[1, 2],
            @lift(Printf.@sprintf("%.2f × %.2f", $(sl_grid[r, 1].value), eeg1[$time_idx])),
            fontsize = math_sz,
            font = :bold,
            color = colors[1],
        )
        Label(row_grid[1, 3], ")  +  (", fontsize = math_sz, font = :bold, color = :black)
        Label(
            row_grid[1, 4],
            @lift(Printf.@sprintf("%.2f × %.2f", $(sl_grid[r, 2].value), eeg2[$time_idx])),
            fontsize = math_sz,
            font = :bold,
            color = colors[2],
        )
        Label(row_grid[1, 5], ")  +  (", fontsize = math_sz, font = :bold, color = :black)
        Label(
            row_grid[1, 6],
            @lift(Printf.@sprintf("%.2f × %.2f", $(sl_grid[r, 3].value), eeg3[$time_idx])),
            fontsize = math_sz,
            font = :bold,
            color = colors[3],
        )
        Label(row_grid[1, 7], ")", fontsize = math_sz, font = :bold, color = :black)

        # Add slight padding to the row to space them out vertically
        rowgap!(exec_box, 15)
        colgap!(row_grid, 0)
    end

    Label(math_grid[2, 6], " = ", fontsize = @lift($math_sz + 10), font = :bold)

    # U Result Vector (Black, since it's the sum of the colors)
    u_box = GridLayout(math_grid[2, 7])
    Box(u_box[1:3, 1], color = (:transparent, 0.0), strokewidth = 3, strokecolor = :black)
    Label(
        u_box[1, 1],
        @lift(Printf.@sprintf("%.2f", u1_full[][$time_idx])),
        fontsize = math_sz,
        font = :bold,
        color = :black,
        padding = (15, 15, 10, 10),
    )
    Label(
        u_box[2, 1],
        @lift(Printf.@sprintf("%.2f", u2_full[][$time_idx])),
        fontsize = math_sz,
        font = :bold,
        color = :black,
        padding = (15, 15, 10, 10),
    )
    Label(
        u_box[3, 1],
        @lift(Printf.@sprintf("%.2f", u3_full[][$time_idx])),
        fontsize = math_sz,
        font = :bold,
        color = :black,
        padding = (15, 15, 10, 10),
    )


    # ── Calculations ──────────────────────────────────────────────────────────
    function update_calcs()
        u1_full[] = (w_obs[1, 1][] .* eeg1) .+ (w_obs[1, 2][] .* eeg2) .+ (w_obs[1, 3][] .* eeg3)
        u2_full[] = (w_obs[2, 1][] .* eeg1) .+ (w_obs[2, 2][] .* eeg2) .+ (w_obs[2, 3][] .* eeg3)
        u3_full[] = (w_obs[3, 1][] .* eeg1) .+ (w_obs[3, 2][] .* eeg2) .+ (w_obs[3, 3][] .* eeg3)
        notify(time_idx) # Force refresh
    end

    # ── Time Slider ───────────────────────────────────────────────────────────
    time_grid = GridLayout(fig[3, 1])
    Label(time_grid[1, 1], "Time ", fontsize = ctrl_sz, font = :bold)
    sl_time = Slider(time_grid[1, 2], range = 1:N, startvalue = 1)
    on(sl_time.value) do val
        time_idx[] = val
    end

    # ── Layout: Plot Columns ──────────────────────────────────────────────────
    plot_grid = GridLayout(fig[4, 1])

    # Column 1: EEG Data (X) - Colored by their respective column!
    axes_x = []
    eeg_data = [eeg1, eeg2, eeg3]
    for i = 1:3
        ax = Axis(plot_grid[i, 1], title = "EEG $i [ X_$i ]", titlesize = title_sz)
        lines!(ax, 1:N, eeg_data[i], color = colors[i], linewidth = 2)
        vlines!(ax, time_idx, color = :lightgray, linewidth = 3)
        ylims!(ax, -5, 5)
        hidedecorations!(ax)
        push!(axes_x, ax)
    end

    # Column 2: Components (U) - Black, as they are the output sum
    t_hist = @lift(1:($time_idx))
    axes_u = []
    comp_full = [u1_full, u2_full, u3_full]

    for i = 1:3
        c_hist = @lift($(comp_full[i])[1:($time_idx)])
        ax = Axis(plot_grid[i, 2], title = "Component $i [ U_$i ]", titlesize = title_sz)
        lines!(ax, t_hist, c_hist, color = :black, linewidth = 3)
        scatter!(ax, time_idx, @lift($(comp_full[i])[$time_idx]), color = :lightgray, markersize = 16)
        xlims!(ax, 1, N)
        ylims!(ax, -10, 10)
        hidedecorations!(ax)
        push!(axes_u, ax)
    end

    # Sizing
    rowsize!(fig.layout, 1, Relative(0.05))
    rowsize!(fig.layout, 2, Relative(0.20))
    rowsize!(fig.layout, 3, Relative(0.05))
    rowsize!(fig.layout, 4, Relative(0.70))

    update_calcs() # bootstrap

    display(fig)
    return fig
end
