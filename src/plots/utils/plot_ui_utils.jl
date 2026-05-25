"""
    _create_menu(fig, options, default, label; kwargs...)

Create a dropdown menu with an adjacent label.
"""
function _create_menu(fig, options, default, label; kwargs...)
    menu = Menu(fig, options = options, default = default, direction = :down, fontsize = 18, width = Auto(), tellwidth = false, kwargs...)
    return hcat(menu, Label(fig, label, fontsize = 22, halign = :left, tellwidth = false))
end

"""
    _build_checkbox_grid!(grid, items, cols, checked_fn; fontsize = 14, label_fn = string)

Build a labelled checkbox grid into `grid`. Returns `Vector{Checkbox}`.
"""
function _build_checkbox_grid!(grid, items, cols, checked_fn; fontsize = 14, label_fn = string)
    checkboxes = Checkbox[]
    for (i, item) in enumerate(items)
        row = ((i - 1) ÷ cols) + 1
        col = ((i - 1) % cols) + 1
        col_base = (col - 1) * 2
        cb = Checkbox(grid[row, col_base + 1], checked = checked_fn(item, i), width = 16, height = 16)
        push!(checkboxes, cb)
        Label(grid[row, col_base + 2], label_fn(item), fontsize = fontsize, halign = :left)
    end
    actual_cols = min(cols, length(items))
    for c in 1:actual_cols
        col_base = (c - 1) * 2
        colgap!(grid, col_base + 1, 5)
        c < actual_cols && colgap!(grid, col_base + 2, 20)
    end
    rowgap!(grid, 5)
    return checkboxes
end

"""
    _popup_button(fig, label, popup_fn)

Return a Button that opens `popup_fn` on click. Used to build control panels with popup windows.
"""
function _popup_button(fig, label, popup_fn)
    btn = Button(fig, label = label, width = Auto(), fontsize = 18)
    on(btn.clicks) do _ popup_fn() end
    return btn
end

"""
    _add_group_buttons!(fig, grid_row, cbs, items, groups)

Add a row of group-selection buttons that set checkboxes to match `predicate`.
"""
function _add_group_buttons!(fig, grid_row, cbs, items, groups)
    area = fig[grid_row, 1] = GridLayout()
    for (j, (lbl, pred)) in enumerate(groups)
        btn = Button(area[1, j], label = lbl)
        on(btn.clicks) do _
            for (i, cb) in enumerate(cbs)
                cb.checked[] = pred(items[i])
            end
        end
    end
end

"""
    _radio_buttons!(btns::Vector, selected_obs::Observable, values::Vector)

Wire buttons as a radio group: clicking one sets `selected_obs` and highlights it.
"""
function _radio_buttons!(btns::Vector, selected_obs::Observable, values::Vector)
    for (btn, val) in zip(btns, values)
        on(btn.clicks) do _
            selected_obs[] = val
            for (b, v) in zip(btns, values)
                b.buttoncolor[] = v == val ? :lightblue : :white
            end
        end
    end
    # Highlight the current default
    for (b, v) in zip(btns, values)
        b.buttoncolor[] = v == selected_obs[] ? :lightblue : :white
    end
end
