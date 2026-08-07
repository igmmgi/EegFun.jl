# === NOTE ON SPELLING: COLOR vs COLOUR ===
# This codebase uses "color" (American spelling) throughout to maintain
# consistency with Makie.jl's API, which uses "color" in all its attributes
# (e.g., colormap, colorbar, strokecolor, etc.). While "colour" may be the
# preferred spelling in some regions, we follow Makie's convention for
# consistency and to avoid confusion when working with Makie's API.
# =============================================================================

"""
    _find_continuous_regions(mask::BitVector, time_points::Vector{Float64})

Find continuous regions where mask is true, returning start and end times.
"""
function _find_continuous_regions(mask::BitVector, time_points::Vector{Float64})

    regions = Vector{Tuple{Float64,Float64}}()
    if isempty(mask) || !any(mask)
        return regions
    end

    in_region = false
    start_idx = 0

    for (i, is_sig) in enumerate(mask)
        if is_sig && !in_region # Start of a new region
            in_region = true
            start_idx = i
        elseif !is_sig && in_region # End of a region
            in_region = false
            push!(regions, (time_points[start_idx], time_points[i-1]))
        end
    end

    # Handle case where region extends to end
    if in_region
        push!(regions, (time_points[start_idx], time_points[end]))
    end

    return regions
end

function _open_save_settings_dialog(fig)
    fig_settings = Figure(size = (550, 350), figure_padding = 20)
    rowgap!(fig_settings.layout, 12)
    colgap!(fig_settings.layout, 15)

    Label(fig_settings[1, 1:2], "Export Figure Settings", fontsize = 18, halign = :center, font = :bold)

    Label(fig_settings[2, 1], "Format:", halign = :left, font = :bold)
    format_menu = Menu(fig_settings[2, 2], options = ["PNG", "PDF", "SVG", "JPG"], default = "PNG", halign = :left, width = 220)

    Label(fig_settings[3, 1], "Resolution (DPI):", halign = :left, font = :bold)
    scale_menu =
        Menu(fig_settings[3, 2], options = ["72 DPI", "150 DPI", "300 DPI", "600 DPI"], default = "72 DPI", halign = :left, width = 220)
    scale_na_label = Label(fig_settings[3, 2], "N/A (Vector Graphics are resolution-independent)", halign = :left, font = :italic)
    scale_na_label.blockscene.visible[] = false

    Label(fig_settings[4, 1], "Output Size (W × H):", halign = :left, font = :bold)
    size_grid = GridLayout(fig_settings[4, 2], halign = :left)
    w_orig, h_orig = size(fig.scene)
    w_orig_cm = round(w_orig * 2.54 / 72.0, digits = 2)
    h_orig_cm = round(h_orig * 2.54 / 72.0, digits = 2)

    w_box = Textbox(size_grid[1, 1], validator = Float64, width = 80)
    Label(size_grid[1, 2], "×", font = :bold)
    h_box = Textbox(size_grid[1, 3], validator = Float64, width = 80)
    Label(size_grid[1, 4], "cm", font = :bold)

    lock_toggle = Toggle(size_grid[1, 5], active = true, halign = :left)
    Label(size_grid[1, 6], "Lock Aspect", font = :bold)
    colgap!(size_grid, 6)

    # Initial values
    w_box.displayed_string[] = string(w_orig_cm)
    w_box.stored_string[] = string(w_orig_cm)
    h_box.displayed_string[] = string(h_orig_cm)
    h_box.stored_string[] = string(h_orig_cm)

    Label(fig_settings[5, 1], "Output Pixels:", halign = :left, font = :bold)
    pixel_label = Label(fig_settings[5, 2], "", halign = :left)

    Label(fig_settings[6, 1], "Transparent Background:", halign = :left, font = :bold)
    transparent_toggle = Toggle(fig_settings[6, 2], active = false, halign = :left)

    export_button = Button(fig_settings[7, 1:2], label = "Select File & Save...", font = :bold)

    function update_pixel_display()
        fmt = format_menu.selection[]
        if fmt in ("PDF", "SVG")
            pixel_label.text[] = "N/A (Vector Format)"
            return
        end

        w_cm = tryparse(Float64, w_box.displayed_string[])
        if isnothing(w_cm) || w_cm <= 0.0
            w_cm = tryparse(Float64, w_box.stored_string[])
        end
        if isnothing(w_cm) || w_cm <= 0.0
            w_cm = w_orig_cm
        end

        h_cm = tryparse(Float64, h_box.displayed_string[])
        if isnothing(h_cm) || h_cm <= 0.0
            h_cm = tryparse(Float64, h_box.stored_string[])
        end
        if isnothing(h_cm) || h_cm <= 0.0
            h_cm = h_orig_cm
        end

        scale_str = scale_menu.selection[]
        dpi_val = 72
        if !isnothing(scale_str)
            parsed_dpi = tryparse(Int, first(split(scale_str)))
            if !isnothing(parsed_dpi)
                dpi_val = parsed_dpi
            end
        end

        w_px = round(Int, (w_cm / 2.54) * dpi_val)
        h_px = round(Int, (h_cm / 2.54) * dpi_val)

        pixel_label.text[] = "$(w_px) × $(h_px) px"
    end

    on(scale_menu.selection) do _
        update_pixel_display()
    end

    on(format_menu.selection) do fmt
        is_vector = fmt in ("PDF", "SVG")
        scale_menu.blockscene.visible[] = !is_vector
        scale_na_label.blockscene.visible[] = is_vector
        update_pixel_display()
    end

    on(w_box.stored_string) do val_str
        w_val = tryparse(Float64, val_str)
        if !isnothing(w_val) && w_val > 0.0
            if lock_toggle.active[]
                h_curr = tryparse(Float64, h_box.stored_string[])
                h_new = round(w_val * (h_orig / w_orig), digits = 2)
                if isnothing(h_curr) || abs(h_curr - h_new) > 0.01
                    h_box.displayed_string[] = string(h_new)
                    h_box.stored_string[] = string(h_new)
                end
            end
            update_pixel_display()
        end
    end

    on(h_box.stored_string) do val_str
        h_val = tryparse(Float64, val_str)
        if !isnothing(h_val) && h_val > 0.0
            if lock_toggle.active[]
                w_curr = tryparse(Float64, w_box.stored_string[])
                w_new = round(h_val * (w_orig / h_orig), digits = 2)
                if isnothing(w_curr) || abs(w_curr - w_new) > 0.01
                    w_box.displayed_string[] = string(w_new)
                    w_box.stored_string[] = string(w_new)
                end
            end
            update_pixel_display()
        end
    end

    on(lock_toggle.active) do locked
        if locked
            w_val = tryparse(Float64, w_box.stored_string[])
            if !isnothing(w_val) && w_val > 0.0
                h_new = round(w_val * (h_orig / w_orig), digits = 2)
                h_box.displayed_string[] = string(h_new)
                h_box.stored_string[] = string(h_new)
                update_pixel_display()
            end
        end
    end

    update_pixel_display()

    on(export_button.clicks) do _
        fmt = format_menu.selection[]
        trans = transparent_toggle.active[]

        # Determine dimensions in cm
        w_str = w_box.displayed_string[]
        w_cm = tryparse(Float64, w_str)
        if isnothing(w_cm) || w_cm <= 0.0
            w_cm = tryparse(Float64, w_box.stored_string[])
        end
        if isnothing(w_cm) || w_cm <= 0.0
            w_cm = w_orig_cm
        end

        h_str = h_box.displayed_string[]
        h_cm = tryparse(Float64, h_str)
        if isnothing(h_cm) || h_cm <= 0.0
            h_cm = tryparse(Float64, h_box.stored_string[])
        end
        if isnothing(h_cm) || h_cm <= 0.0
            h_cm = h_orig_cm
        end

        scale_str = scale_menu.selection[]
        dpi_val = 72
        if !isnothing(scale_str)
            parsed_dpi = tryparse(Int, first(split(scale_str)))
            if !isnothing(parsed_dpi)
                dpi_val = parsed_dpi
            end
        end

        # Calculate final size in points
        w_target_pt = (w_cm / 2.54) * 72.0
        h_target_pt = (h_cm / 2.54) * 72.0

        ext_map = Dict("PNG" => ".png", "PDF" => ".pdf", "SVG" => ".svg", "JPG" => ".jpg")
        filter_ext = ext_map[fmt]

        @async begin
            try
                # Close settings window
                screen = Makie.getscreen(fig_settings.scene)
                if !isnothing(screen)
                    close(screen)
                end

                # Open save file dialog with specific filter
                path = NativeFileDialog.save_file(pwd(); filterlist = lowercase(fmt))
                if !isempty(path)
                    # Append correct extension if missing
                    if !endswith(lowercase(path), filter_ext)
                        path = path * filter_ext
                    end

                    # Handle transparency
                    orig_bg = fig.scene.backgroundcolor[]
                    orig_content_bgs = Dict()
                    if trans
                        fig.scene.backgroundcolor[] = RGBAf(0, 0, 0, 0)
                        for content in fig.content
                            if hasproperty(content, :backgroundcolor) && content.backgroundcolor isa Observable
                                orig_content_bgs[content] = content.backgroundcolor[]
                                content.backgroundcolor[] = RGBAf(0, 0, 0, 0)
                            end
                        end
                    end

                    # Store original size and resize figure temporarily
                    w_orig_sz, h_orig_sz = size(fig.scene)
                    resize!(fig.scene, (w_target_pt, h_target_pt))

                    try
                        if filter_ext in (".pdf", ".svg")
                            Makie.save(path, fig; backend = CairoMakie, pt_per_unit = 1.0)
                        else
                            Makie.save(path, fig; px_per_unit = dpi_val / 72.0)
                        end
                        @info "Figure saved to: $path"
                    finally
                        # Restore original figure size
                        resize!(fig.scene, (w_orig_sz, h_orig_sz))

                        if trans
                            fig.scene.backgroundcolor[] = orig_bg
                            for (content, bg) in orig_content_bgs
                                content.backgroundcolor[] = bg
                            end
                        end
                    end
                end
            catch err
                @warn "Failed to save figure: $err"
            end
        end
    end

    _display_popup(fig_settings)
end



function _display_figure(fig)
    # Register keyboard shortcut for saving the figure
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press && event.key == Keyboard.s
            _open_save_settings_dialog(fig)
        end
    end

    display(fig)
end

"""
    _display_popup(fig)

Display a figure in a new, independent popup window if the backend supports it (e.g., GLMakie).
Otherwise, falls back to the default `display(fig)`.
"""
function _display_popup(fig)
    backend = Makie.current_backend()
    if !isnothing(backend) && string(backend) == "GLMakie"
        screen = backend.Screen()
        display(screen, fig)
    else
        display(fig)
    end
end

"""
    _set_window_title(title::String)

Set the window title for GLMakie. Does nothing for CairoMakie (which doesn't support window titles).
"""
function _set_window_title(title::String)
    if string(Makie.current_backend()) == "GLMakie"
        Makie.current_backend().activate!(title = title)
    end
    # Note: CairoMakie doesn't support window titles.
end

"""
    _create_figure_with_axis(data; title_suffix::String = "", figure_kwargs...)

Create a Figure and Axis with window title set from data.

# Arguments
- `data`: Data object to generate window title from (passed to `_generate_window_title`)
- `title_suffix::String`: Optional suffix to append to the generated title
- `figure_kwargs...`: Keyword arguments passed to `Figure()` constructor

# Returns
- `fig::Figure`: The created figure
- `ax::Axis`: The created axis at position `[1, 1]`
"""
function _create_figure_with_axis(data; title_suffix::String = "", figure_kwargs...)
    title = _generate_window_title(data)
    if !isempty(title_suffix)
        title = title * title_suffix
    end
    _set_window_title(title)
    fig = Figure(; figure_kwargs...)
    ax = Axis(fig[1, 1])
    return (fig = fig, axes = [ax])
end

"""
    _get_colorbar_defaults()

Get all default values for Colorbar attributes by creating a single Colorbar instance.
Returns a dictionary mapping attribute names to their default values.
"""
function _get_colorbar_defaults()
    # Create a minimal figure
    fig = Figure()
    cb = Colorbar(fig)

    # Get all attribute values at once
    defaults = Dict{Symbol,Any}()
    for attr in propertynames(Colorbar)
        defaults[attr] = getproperty(cb, attr)
    end

    return defaults
end

# Cache the colorbar defaults
const COLORBAR_DEFAULTS = _get_colorbar_defaults()

"""
    _get_legend_defaults()

Get all default values for Legend attributes by creating a temporary Legend instance.
Returns a dictionary mapping attribute names to their default values.
"""
function _get_legend_defaults()
    # Create a minimal figure and axis with a dummy plot that has labels
    # This is necessary because axislegend() requires plots with labels
    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax, [1, 2], [1, 2], label = "dummy")  # Create a plot with a label
    leg = axislegend(ax)

    # Get all attribute values at once
    defaults = Dict{Symbol,Any}()
    for attr in propertynames(Legend)
        defaults[attr] = getproperty(leg, attr)
    end

    return defaults
end

# Cache the legend defaults
const LEGEND_DEFAULTS = _get_legend_defaults()

"""
    _extract_colorbar_kwargs!(plot_kwargs::Dict{Symbol, Any})

Extract all colorbar-related parameters from plot_kwargs and return a clean dictionary
suitable for passing to Colorbar constructor.

# Arguments
- `plot_kwargs`: Dictionary of plot parameters (modified in-place)

# Returns
- `Dict{Symbol, Any}`: Cleaned colorbar parameters with invalid attributes removed
"""
function _extract_colorbar_kwargs!(plot_kwargs::Dict{Symbol,Any})
    colorbar_kwargs = Dict{Symbol,Any}()
    colorbar_attrs = propertynames(Colorbar)

    for attr in colorbar_attrs
        colorbar_key = Symbol("colorbar_$(attr)")
        if haskey(plot_kwargs, colorbar_key)
            value = pop!(plot_kwargs, colorbar_key)
            if !isnothing(value)  # Only add if not the default nothing
                colorbar_kwargs[attr] = value
            end
        end
    end

    # Note: These attributes conflict with internal Makie colorbar creation if passed directly.
    # These cannot be passed to colorbar kwargs
    pop!(colorbar_kwargs, :colormap, nothing)
    pop!(colorbar_kwargs, :limits, nothing)
    pop!(colorbar_kwargs, :highclip, nothing)
    pop!(colorbar_kwargs, :lowclip, nothing)

    return colorbar_kwargs
end

"""
    _extract_legend_kwargs(plot_kwargs::Dict{Symbol, Any}; exclude_positioning::Bool=false)

Extract legend-related parameters from plot_kwargs and return a new dictionary
suitable for passing to axislegend().

Does not mutate plot_kwargs - only reads from it.

# Arguments
- `plot_kwargs`: Dictionary of plot parameters (read-only)
- `exclude_positioning`: If true, exclude positioning attributes (halign, valign, alignmode) 
  that conflict with the `position` parameter in axislegend()

# Returns
- `Dict{Symbol, Any}`: New dictionary containing extracted legend parameters (with legend_ prefix removed)
"""
function _extract_legend_kwargs(plot_kwargs::Dict{Symbol,Any}; exclude_positioning::Bool = false)
    legend_kwargs = Dict{Symbol,Any}()
    legend_attrs = propertynames(Legend)

    for attr in legend_attrs
        legend_key = Symbol("legend_$(attr)")
        if haskey(plot_kwargs, legend_key)
            value = plot_kwargs[legend_key]
            if !isnothing(value)
                legend_kwargs[attr] = value
            end
        end
    end

    # axislegend() uses its own `position` parameter for placement. If halign/valign/alignmode
    # are also present, they conflict and `position` is silently ignored. Remove them so that
    # legend_position works as expected.
    for attr in (:halign, :valign, :alignmode)
        pop!(legend_kwargs, attr, nothing)
    end

    return legend_kwargs
end

"""
    _extract_layout_kwargs(plot_kwargs::Dict{Symbol, Any})

Extract layout-related parameters from plot_kwargs and return a new dictionary
suitable for passing to create_layout().

Does not mutate plot_kwargs - only reads from it.

# Arguments
- `plot_kwargs`: Dictionary of plot parameters (read-only)

# Returns
- `Dict{Symbol, Any}`: New dictionary containing extracted layout parameters (with layout_ prefix removed)
"""
function _extract_layout_kwargs(plot_kwargs::Dict{Symbol,Any})
    layout_kwargs = Dict{Symbol,Any}()

    # Get all layout parameter names from LAYOUT_KWARGS
    # These are the base names (keys in LAYOUT_KWARGS, e.g., :topo_plot_width, :grid_rowgap)
    layout_param_names = keys(LAYOUT_KWARGS)

    # For each known layout parameter, check if it exists with layout_ prefix
    for param_name in layout_param_names
        layout_key = Symbol("layout_$(param_name)")
        if haskey(plot_kwargs, layout_key)
            value = plot_kwargs[layout_key]
            if !isnothing(value)
                layout_kwargs[param_name] = value
            end
        end
    end

    return layout_kwargs
end

# === AXIS STYLING FUNCTIONS ===

"""
    _set_axis_grid!(ax; xgrid = false, ygrid = false, xminorgrid = false, yminorgrid = false)

Apply grid settings to the axis.

# Arguments
- `ax`: Makie Axis object

# Keyword Arguments
- `xgrid`: Whether to show x-axis grid
- `ygrid`: Whether to show y-axis grid  
- `xminorgrid`: Whether to show x-axis minor grid
- `yminorgrid`: Whether to show y-axis minor grid
"""
function _set_axis_grid!(ax; xgrid = false, ygrid = false, xminorgrid = false, yminorgrid = false)
    ax.xgridvisible = xgrid
    ax.ygridvisible = ygrid
    ax.xminorgridvisible = xminorgrid
    ax.yminorgridvisible = yminorgrid
end

"""
    _set_axis_properties!(ax; xlim = nothing, ylim = nothing, xlabel = "", ylabel = "", yreversed = false)

Apply axis limits, labels, and direction to the axis.

# Arguments
- `ax`: Makie Axis object

# Keyword Arguments
- `xlim`: X-axis limits as (min, max) tuple or nothing for auto-scaling
- `ylim`: Y-axis limits as (min, max) tuple or nothing for auto-scaling
- `xlabel`: Label for x-axis (default: empty string)
- `ylabel`: Label for y-axis (default: empty string)
- `yreversed`: Whether to reverse the y-axis (default: false)
"""
function _set_axis_properties!(
    ax;
    xlim = nothing,
    ylim = nothing,
    xlabel = "",
    ylabel = "",
    yreversed = false,
    scale_x_value = nothing,
    scale_y_value = nothing,
)

    # Set axis labels
    ax.xlabel = xlabel
    ax.ylabel = ylabel
    ax.yreversed = yreversed

    # Set exact tick spacing if provided by the user
    if !isnothing(scale_x_value)
        ax.xticks = (vmin, vmax) -> begin
            start_x = ceil(vmin / scale_x_value) * scale_x_value
            end_x = floor(vmax / scale_x_value) * scale_x_value
            ticks = collect(start_x:scale_x_value:(end_x+1e-9))
            labels = string.(round.(ticks, sigdigits = 4))
            return ticks, labels
        end
    end

    if !isnothing(scale_y_value)
        ax.yticks = (vmin, vmax) -> begin
            start_y = ceil(vmin / scale_y_value) * scale_y_value
            end_y = floor(vmax / scale_y_value) * scale_y_value
            ticks = collect(start_y:scale_y_value:(end_y+1e-9))
            labels = string.(round.(ticks, sigdigits = 4))
            return ticks, labels
        end
    end

    # Set axis limits
    !isnothing(xlim) && xlims!(ax, xlim[1], xlim[2])
    !isnothing(ylim) && ylims!(ax, ylim[1], ylim[2])

end

"""
    _set_origin_lines!(ax; add_xy_origin = true, color = :gray, linewidth = 0.5, alpha = 0.7)

Add origin lines at x=0 and y=0 to the axis.

# Arguments
- `ax`: Makie Axis object

# Keyword Arguments
- `add_xy_origin`: Whether to add origin lines at x=0 and y=0
- `color`: Color of the origin lines
- `linewidth`: Line width of the origin lines
- `alpha`: Transparency of the origin lines
"""
function _set_origin_lines!(ax; add_xy_origin = true, color = :gray, linewidth = 1, alpha = 0.8)
    if add_xy_origin
        hlines!(ax, 0, color = color, linewidth = linewidth, alpha = alpha)
        vlines!(ax, 0, color = color, linewidth = linewidth, alpha = alpha)
    end
end

"""
    _add_origin_scale_indicator!(ax; show_scale_indicator = false, scale_x_value = 0.1, scale_y_value = 10.0, scale_x_label = "s", scale_y_label = "μV")

Replace outer spines and ticks with crosshair axes passing through the origin.
"""
function _add_origin_scale_indicator!(
    ax;
    axis_type = :standard,
    scale_x_value = nothing,
    scale_y_value = nothing,
    xlabel = "s",
    ylabel = "μV",
)
    if axis_type == :origin
        # Hide standard spines and decorations
        hidedecorations!(ax)
        hidespines!(ax)

        # Draw the thick scale bars starting from (0,0)?
        # No, for crosshair axes, the origin lines themselves ARE the axes.
        # We can just keep the hlines! and vlines! from _set_origin_lines! (which are already drawn)
        # But we'll add our own thick crosshair lines if we want them to stand out
        hlines!(ax, 0, color = :black, linewidth = 1, overdraw = true)
        vlines!(ax, 0, color = :black, linewidth = 1, overdraw = true)

        # Compute tick positions dynamically based on scale_x_value and scale_y_value
        x_ticks_obs = lift(ax.finallimits) do lims
            min_x, max_x = lims.origin[1], lims.origin[1] + lims.widths[1]
            start_x = ceil(min_x / scale_x_value) * scale_x_value
            end_x = floor(max_x / scale_x_value) * scale_x_value
            # Floating point issues can make range collection problematic, so we use a small epsilon
            collect(start_x:scale_x_value:(end_x+1e-9))
        end

        y_ticks_obs = lift(ax.finallimits) do lims
            min_y, max_y = lims.origin[2], lims.origin[2] + lims.widths[2]
            start_y = ceil(min_y / scale_y_value) * scale_y_value
            end_y = floor(max_y / scale_y_value) * scale_y_value
            collect(start_y:scale_y_value:(end_y+1e-9))
        end

        # Format labels
        x_labels_obs = lift(x_ticks_obs) do ticks
            [string(round(t, sigdigits = 3)) for t in ticks]
        end
        y_labels_obs = lift(y_ticks_obs) do ticks
            [string(round(t, sigdigits = 3)) for t in ticks]
        end

        # Tick markers (small lines crossing the axes)
        x_pts_obs = lift(x_ticks_obs) do ticks
            [Point2f(t, 0.0) for t in ticks if abs(t) > 1e-9]
        end
        y_pts_obs = lift(y_ticks_obs) do ticks
            [Point2f(0.0, t) for t in ticks if abs(t) > 1e-9]
        end
        scatter!(ax, x_pts_obs, marker = '|', markersize = 10, color = :black, overdraw = true, xautolimits = false, yautolimits = false)
        scatter!(ax, y_pts_obs, marker = '-', markersize = 10, color = :black, overdraw = true, xautolimits = false, yautolimits = false)

        # Text labels
        x_lbls = lift(x_labels_obs, x_ticks_obs) do labels, ticks
            [l for (l, t) in zip(labels, ticks) if abs(t) > 1e-9]
        end
        y_lbls = lift(y_labels_obs, y_ticks_obs) do labels, ticks
            [l for (l, t) in zip(labels, ticks) if abs(t) > 1e-9]
        end
        text!(
            ax,
            x_pts_obs,
            text = x_lbls,
            align = (:center, :top),
            offset = (0, -8),
            color = :black,
            fontsize = 12,
            overdraw = true,
            xautolimits = false,
            yautolimits = false,
        )
        text!(
            ax,
            y_pts_obs,
            text = y_lbls,
            align = (:right, :center),
            offset = (-8, 0),
            color = :black,
            fontsize = 12,
            overdraw = true,
            xautolimits = false,
            yautolimits = false,
        )

        # Axis labels (Time and Amplitude)
        # Place them at the end of the axes
        x_axis_label_pos = lift(ax.finallimits) do lims
            Point2f(lims.origin[1] + lims.widths[1], 0.0)
        end
        y_axis_label_pos = lift(ax.finallimits) do lims
            Point2f(0.0, lims.origin[2] + lims.widths[2])
        end
        text!(
            ax,
            x_axis_label_pos,
            text = string(" ", xlabel),
            align = (:left, :center),
            color = :black,
            fontsize = 14,
            overdraw = true,
            xautolimits = false,
            yautolimits = false,
        )
        text!(
            ax,
            y_axis_label_pos,
            text = string(ylabel, " "),
            align = (:center, :bottom),
            color = :black,
            fontsize = 14,
            overdraw = true,
            xautolimits = false,
            yautolimits = false,
        )
    end
end

# === WINDOW TITLE AND STRING ABBREVIATION UTILITIES ===

"""
    _split_into_parts(s::String)

Split a string into parts on uppercase letters, underscores, and trailing digits.
Returns a vector of parts, e.g., "ExampleCondition1" -> ["Example", "Condition", "1"]

# Arguments
- `s::String`: String to split

# Returns
- `Vector{String}`: Vector of parts
"""
function _split_into_parts(s::String)
    # Use regex to split: match words (uppercase or lowercase, with optional trailing digits), underscores, or standalone digits
    # Pattern: ([A-Z][a-z]*|[a-z]+)(\d*) - word with optional trailing digits, (_) - underscore, (\d+) - standalone digits
    parts = String[]
    pattern = r"([A-Z][a-z]*|[a-z]+)(\d*)|(_)|(\d+)"

    for m in eachmatch(pattern, s)
        if !isnothing(m.captures[3])  # Underscore
            push!(parts, "_")
        elseif !isnothing(m.captures[4])  # Standalone digits
            push!(parts, m.captures[4])
        else  # Word with optional digits
            word = m.captures[1]
            digits = m.captures[2]
            push!(parts, word)
            !isempty(digits) && push!(parts, digits)
        end
    end

    return parts
end

"""
    _generate_window_title(datasets; max_total_length::Int = 80, max_name_length::Int = 30)

Generate a window title from datasets, including file names and condition names.
Handles single and multiple conditions, with automatic shortening if needed.

Works with any type that has `.file` and `.condition_name` fields.

If all datasets share the same file name, uses format: "FileName: Cond1, Cond2, Cond3"
If datasets have different file names, uses format: "FileX:Cond1, FileY:Cond2"

# Arguments
- `datasets`: Vector of datasets (must have `.file` and `.condition_name` fields)
- `max_total_length::Int`: Maximum total length of the result string
- `max_name_length::Int`: Maximum length for each individual "file:condition" pair

# Returns
- `String`: Window title string
"""
function _generate_window_title(datasets::Vector{<:EegData}; max_total_length::Int = 80, max_name_length::Int = 20)

    isempty(datasets) && return ""
    length(datasets) == 1 && return "$(datasets[1].file):$(datasets[1].condition_name)"

    # Check if all datasets have the same file name
    first_file = datasets[1].file
    all_same_file = all(dataset.file == first_file for dataset in datasets)

    if all_same_file
        condition_names = [data.condition_name for data in datasets]
        condition_str = _shorten_condition_names(
            condition_names;
            max_total_length = max_total_length - length(first_file) - 2,
            max_name_length = max_name_length,
        )
        return "$first_file: $condition_str"
    else
        file_condition_pairs = ["$(data.file):$(data.condition_name)" for data in datasets]
        return _shorten_condition_names(file_condition_pairs; max_total_length = max_total_length, max_name_length = max_name_length)
    end
end

_generate_window_title(datasets::EegData) = "$(datasets.file):$(datasets.condition_name)"

"""
    _generate_window_title(datasets::ContinuousData)

Generate a window title from ContinuousData, using just the filename.

# Arguments
- `datasets::ContinuousData`: ContinuousData object

# Returns
- `String`: Window title string (just the filename)
"""
_generate_window_title(datasets::ContinuousData) = datasets.file

"""
    _generate_window_title(datasets::Vector{ContinuousData}; 
                          max_total_length::Int = 80,
                          max_name_length::Int = 30)

Generate a window title from a vector of ContinuousData objects.

If all datasets share the same file name, uses format: "FileName"
If datasets have different file names, uses format: "File1, File2, File3"

# Arguments
- `datasets::Vector{ContinuousData}`: Vector of ContinuousData objects
- `max_total_length::Int`: Maximum total length of the result string
- `max_name_length::Int`: Maximum length for each individual filename

# Returns
- `String`: Window title string
"""
function _generate_window_title(datasets::Vector{ContinuousData}; max_total_length::Int = 80, max_name_length::Int = 20)
    isempty(datasets) && return ""
    length(datasets) == 1 && return datasets[1].file

    # Check if all datasets have the same file name
    first_file = datasets[1].file
    all_same_file = all(dataset.file == first_file for dataset in datasets)

    all_same_file && return first_file

    file_names = [data.file for data in datasets]
    return _shorten_condition_names(file_names; max_total_length = max_total_length, max_name_length = max_name_length)

end

"""
    _abbreviate_name(name::String, common_prefix_parts::Vector{String})

Create an intelligent abbreviation of a name by:
- Abbreviating the common prefix parts (taking first letters)
- Preserving the unique suffix parts

# Arguments
- `name::String`: The name to abbreviate
- `common_prefix_parts::Vector{String}`: Common prefix parts found across all names

# Returns
- `String`: Abbreviated name
"""
function _abbreviate_name(name::String, common_prefix_parts::Vector{String})
    isempty(common_prefix_parts) && return name

    # abbreviate: take first 2-3 letters of each (or first letter if short)
    name_parts = _split_into_parts(name)
    abbrev_prefix = ""
    for part in common_prefix_parts
        if part == "_"
            abbrev_prefix *= "_"
        elseif !isempty(part)
            if isdigit(part[1])
                abbrev_prefix *= part  # Keep numbers as-is
            else
                # Take first 2-3 letters, or first letter if part is very short
                n = length(part) >= 3 ? 3 : (length(part) >= 2 ? 2 : 1)
                abbrev_prefix *= uppercase(part[1]) * part[2:min(n, length(part))]
            end
        end
    end

    # Get the unique suffix parts (everything after the common prefix)
    unique_suffix = join(name_parts[(length(common_prefix_parts)+1):end], "")

    return abbrev_prefix * unique_suffix
end

"""
    _find_common_prefix_parts(names::Vector{String})

Find common prefix parts across a list of names, splitting on uppercase letters and underscores.

# Arguments
- `names::Vector{String}`: Vector of names to analyze

# Returns
- `Vector{String}`: Common prefix parts
"""
function _find_common_prefix_parts(names::Vector{String})
    isempty(names) && return String[]
    length(names) == 1 && return String[]

    # Split all names into parts
    all_parts = [_split_into_parts(name) for name in names]

    # Find common prefix parts
    first_parts = all_parts[1]
    common_parts = String[]

    for (i, part) in enumerate(first_parts)
        if all(get(parts, i, nothing) == part for parts in all_parts)
            push!(common_parts, part)
        else
            break
        end
    end

    # Don't use if too short (less than 2 parts)
    return length(common_parts) >= 2 ? common_parts : String[]
end

"""
    _abbreviate_parts(parts::Vector{String})

Abbreviate a list of parts by taking first 2-3 letters of each word part.
"""
function _abbreviate_parts(parts::Vector{String})
    abbrev = ""
    for part in parts
        if part == "_"
            abbrev *= "_"
        elseif !isempty(part)
            if isdigit(part[1])
                abbrev *= part
            else
                n = min(3, length(part))
                abbrev *= uppercase(part[1]) * part[2:n]
            end
        end
    end
    return abbrev
end


"""
    _shorten_condition_names(condition_names::Vector{String}; 
                            max_total_length::Int = 40,
                            max_name_length::Int = 30,
                            separator::String = ", ",
                            show_ends::Int = 3)

Shorten a list of condition names while preserving identification.
Uses intelligent abbreviation to create shorter but still identifiable names.
If the total is still too long, shows first and last N conditions with "..." in between.

# Arguments
- `condition_names::Vector{String}`: Vector of condition names
- `max_total_length::Int`: Maximum total length of the result string
- `max_name_length::Int`: Maximum length for each individual name (after abbreviation)
- `separator::String`: Separator between names
- `show_ends::Int`: Number of conditions to show at start/end when truncating

# Returns
- `String`: Shortened condition name string
"""
function _shorten_condition_names(
    condition_names::Vector{String};
    max_total_length::Int = 80,
    max_name_length::Int = 30,
    separator::String = ", ",
    show_ends::Int = 3,
)
    isempty(condition_names) && return ""
    length(condition_names) == 1 && return condition_names[1]

    # Try intelligent abbreviation
    common_prefix_parts = _find_common_prefix_parts(condition_names)
    if isempty(common_prefix_parts)
        abbreviated_names = [_abbreviate_parts(_split_into_parts(name)) for name in condition_names]
    else
        abbreviated_names = [_abbreviate_name(name, common_prefix_parts) for name in condition_names]
    end

    # Use original if abbreviation didn't help (less than 20% shorter)
    avg_original = sum(length(name) for name in condition_names) / length(condition_names)
    avg_abbreviated = sum(length(name) for name in abbreviated_names) / length(abbreviated_names)
    if avg_abbreviated >= avg_original * 0.8
        abbreviated_names = condition_names
    end

    # Truncate individual names and join
    final_names = [length(n) > max_name_length ? n[1:max_name_length] * "…" : n for n in abbreviated_names]
    full_string = join(final_names, separator)

    # If still too long, show first N and last N
    if length(full_string) > max_total_length && length(condition_names) > 2 * show_ends
        first_part = join(final_names[1:show_ends], separator)
        last_part = join(final_names[(end-show_ends+1):end], separator)
        return first_part * separator * "…" * separator * last_part
    end

    return full_string
end

# === BASELINE HELPER FUNCTIONS ===

"""
    _extract_baseline_values(interval::Interval)

Extract numeric baseline values from interval using multiple dispatch.
Returns (start, stop) tuple or (nothing, nothing) if interval is nothing.
"""
_extract_baseline_values(::Nothing) = (nothing, nothing)
_extract_baseline_values(interval::AbstractRange) = (first(interval), last(interval))
_extract_baseline_values(interval::Tuple) = (interval[1], interval[2])

"""
    _parse_baseline_values(start_str::String, stop_str::String)

Parse numeric baseline values from textbox strings.
Returns (start, stop) tuple or (-Inf, Inf) for empty strings, or (nothing, nothing) on parse error.
"""
function _parse_baseline_values(start_str::String, stop_str::String)
    (start_str == " " || stop_str == " ") && return -Inf, Inf
    try
        return parse(Float64, start_str), parse(Float64, stop_str)
    catch e
        @minimal_warning "Invalid baseline values: $e"
        return nothing, nothing
    end
end

"""
    _create_baseline_textbox(layout, row, label, obs, placeholder, width)

Helper to create and connect a textbox to an observable.
Returns the created Textbox.
"""
function _create_baseline_textbox(layout, row, label, obs, placeholder, width)
    Label(layout[row, 1], label, width = 60)
    tb = Textbox(layout[row, 2], placeholder = placeholder, width = width)
    tb.stored_string[] = obs[]
    hasproperty(tb, :displayed_string) ? connect!(obs, tb.displayed_string) : connect!(obs, tb.stored_string)
    return tb
end
