# Shared plot kwargs dictionaries, Makie backend hooks, and axis/colormap utilities.

######################################
######### Types & Constants ##########
######################################

mutable struct MakieExtensionState
    glmakie_active::Bool
    create_screen::Function
    display_screen::Function
    save_vector::Function
    set_title::Function
end

const MAKIE_EXT_STATE = MakieExtensionState(
    false,
    (; size = nothing) -> @minimal_error("Opening interactive windows requires GLMakie. Please run `using GLMakie`."),
    (fig; size = nothing) -> begin
        try
            screen = MAKIE_EXT_STATE.create_screen(; size = size)
            display(screen, fig)
            return screen
        catch e
            @debug "Could not open interactive screen, falling back to display()" exception = e
            display(fig)
            return nothing
        end
    end,
    (path, fig; pt_per_unit = 1.0) -> begin
        if string(Makie.current_backend()) != "CairoMakie"
            @minimal_warning "Saving vector figures (.pdf/.svg) with optimal quality requires CairoMakie. Please run `using CairoMakie`. Falling back to default Makie save."
        end
        Makie.save(path, fig; pt_per_unit = pt_per_unit)
    end,
    (title::String) -> begin
        backend = Makie.current_backend()
        if string(backend) == "GLMakie" && hasproperty(backend, :activate!)
            try
                backend.activate!(title = title)
            catch e
                @debug "Could not set window title" exception = e
            end
        end
    end,
)

const _DIVERGING_COLORMAPS = Set{Symbol}([
    :coolwarm,
    :RdBu,
    :RdYlBu,
    :RdYlGn,
    :RdGy,
    :PuOr,
    :PRGn,
    :PiYG,
    :BrBG,
    :Spectral,
    :bwr,
    :seismic,
    :balance,
    :curl,
    :delta,
    :vik,
    :roma,
    :broc,
    :cork,
    :bam,
    :berlin,
    :tofino,
    :lisbon,
    :managua,
])

const _SPLIT_PARTS_RE = r"([A-Z][a-z]*|[a-z]+)(\d*)|(_)|(\d+)"

const _SIGNED_TF_BASELINE_METHODS = (:db, :percent, :relchange, :normchange, :zscore, :absolute)

# Shared kwargs for all topography plots (ICA, ERP, and time-frequency)
const PLOT_TOPOGRAPHY_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :display_plot => (true, "Whether to display the plot"),
    :figure_title => ("", "Title drawn at the top of the entire figure canvas"),
    :figure_title_fontsize => (24, "Font size for figure title"),
    :interactive => (true, "Whether to enable interactive features"),
    :theme_fontsize => (nothing, "Font size for theme"),
    :zoom_step => (0.2, "Fractional zoom step for arrow keys (e.g. 0.2 means 20% zoom in/out)"),
    :figure_padding => ((10, 30, 10, 10), "Padding around entire figure as (left, right, bottom, top) tuple (in pixels)"),
    :method => (
        :thin_plate,
        "Interpolation method: :multiquadratic, :inverse_multiquadratic, :gaussian, :inverse_quadratic, :thin_plate, :polyharmonic, :shepard, :nearest, :spherical_spline. See https://eljungsk.github.io/ScatteredInterpolation.jl/dev/methods/ for details on methods.",
    ),
    :colormap => (nothing, "Colormap for the topography"),
    :gridscale => (75, "Grid resolution for interpolation"),
    :dims => (nothing, "Grid dimensions (rows, cols). If nothing, calculates best square-ish grid"),
    :colorrange => (nothing, "Color scale limits as (min, max) tuple (nothing for auto)"),
    :num_levels => (20, "Number of contour levels (for ICA plots)."),
    :plot_title => (nothing, "Plot title"),
    :show_title => (true, "Whether to show the plot title"),
    :plot_title_position => (
        nothing,
        "Relative (x, y) coordinates for the plot title (e.g., (0.5, 0.95)). If provided, the title is drawn inside the axis.",
    ),
    :plot_title_align => ((:center, :top), "Alignment of the inner plot title"),
    :plot_title_fontsize => (16, "Font size for the title"),
    :time_unit =>
        (:s, "Time unit for display labels (:s or :ms). Only affects title strings — all intervals and selections remain in seconds."),
    :head_color => (nothing, "Color of the head shape outline."),
    :head_linewidth => (2, "Line width of the head shape outline."),
    :head_radius => (1.0, "Radius of the head shape in mm."),
    :point_plot => (true, "Whether to plot electrode points."),
    :point_marker => (:circle, "Marker style for electrode points."),
    :point_markersize => (12, "Size of electrode point markers."),
    :point_color => (nothing, "Color of electrode points."),
    :label_plot => (true, "Whether to plot electrode labels."),
    :label_fontsize => (20, "Font size for electrode labels."),
    :label_color => (nothing, "Color of electrode labels."),
    :label_xoffset => (0, "X-axis offset for electrode labels."),
    :label_yoffset => (0, "Y-axis offset for electrode labels."),
    :colorbar_kwargs => (
        NamedTuple(),
        "Additional kwargs passed directly to the Makie Colorbar block (see Makie's [Colorbar documentation](https://docs.makie.org/stable/reference/blocks/colorbar/) for available attributes).",
    ),
    :colorbar_plot => (true, "Whether to show the colorbar"),
    :colorbar_position => (:right, "Position of the colorbar (:right, :left, :top, :bottom, or tuple)"),
    :colorbar_label => ("", "Label for the colorbar"),
    :colorbar_plot_numbers => ([], "Plot indices for which to show colorbars. Empty list shows colorbars for all plots."),
    :highlight_channels => (
        nothing,
        "Highlight channel groups as a NamedTuple or Vector of NamedTuples. Each group: (channels=[:Cz, :Pz], color=:white, size=8, marker=:circle)",
    ),
    :use_global_scale => (false, "Do multiple ICA topoplots share the same color scale based on min/max across all components?"),
    :component_selection => (components(), "Function that returns boolean vector for component filtering"),
)

const PLOT_ICA_QUALITY_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :display_plot => (true, "Whether to display the plot"),
    :plot_title => (nothing, "Custom plot title override"),
    :show_title => (true, "Whether to show the plot title"),
    :legend => (true, "Whether to display the legend"),
    :z_threshold => (3.0, "Z-score threshold for identifying components"),
    :min_harmonic_power => (1.5, "Minimum harmonic power threshold for line noise detection"),
    :max_ibi_std_s => (0.2, "Maximum inter-beat interval standard deviation for ECG detection"),
    :min_peak_ratio => (0.7, "Minimum peak ratio for ECG detection"),
    :plot_title_fontsize => (16, "Font size for plot titles"),
    :label_fontsize => (14, "Font size for axis labels"),
    :tick_fontsize => (12, "Font size for tick labels"),
    :legend_fontsize => (12, "Font size for legend"),
    :linewidth => (2, "Line width for plot lines"),
    :line_alpha => (0.8, "Transparency for plot lines"),
    :threshold_style => (:dash, "Line style for threshold lines"),
    :zscore_threshold_color => (:red, "Color for z-score threshold lines"),
    :zscore_threshold_style => (:dash, "Line style for z-score threshold lines"),
    :marker_size => (8, "Size of markers"),
    :marker_alpha => (0.7, "Transparency of markers"),
    :correlation_color => (:steelblue, "Color for correlation plots"),
    :threshold_color => (:red, "Color for threshold lines"),
    :identified_color => (:green, "Color for identified components"),
    :rejected_color => (:red, "Color for rejected components"),
)

const PLOT_ERP_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Display the plot (true/false)"),
    :window_title => ("", "Title for the OS plot window. If empty, it's generated automatically."),
    :figure_title => ("", "Title drawn at the top of the entire figure canvas"),
    :figure_title_fontsize => (24, "Font size for figure title"),
    :interactive => (true, "Enable interactive features (true/false)"),
    :zoom_step => (0.2, "Fractional zoom step for arrow keys (e.g. 0.2 means 20% zoom in/out)"),
    :selection_color => (:blue, "Color for interactive selection rectangles"),
    :selection_alpha => (0.3, "Alpha (transparency) for interactive selection rectangles"),
    :theme_fontsize => (nothing, "Font size for theme"),

    # Axis limits and labels
    :xlim => (nothing, "X-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :ylim => (nothing, "Y-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :xlabel => ("Time (s)", "Label for x-axis"),
    :ylabel => ("μV", "Label for y-axis"),
    :yreversed => (false, "Whether to reverse the y-axis"),
    :xticks => (nothing, "X-axis tick positions (e.g., -0.1:0.1:0.8 or [0, 0.2, 0.4]). If nothing, automatically determined"),
    :yticks => (nothing, "Y-axis tick positions (e.g., -4:2:4 or [-2, 0, 2]). If nothing, automatically determined"),
    :time_unit => (
        :s,
        "Time unit for x-axis display (:s or :ms). Only affects axis labels and tick formatting — all intervals remain in seconds.",
    ),

    # Title
    :plot_title => (nothing, "Plot title"),
    :show_title => (true, "Show plot title (true/false)"),
    :plot_title_fontsize => (16, "Font size for plot titles"),
    :plot_title_position => (
        nothing,
        "Relative (x, y) coordinates for the plot title (e.g., (0.5, 0.95)). If provided, the title is drawn inside the axis.",
    ),
    :plot_title_align => ((:center, :top), "Alignment of the inner plot title"),

    # Line styling
    :linewidth => (nothing, "Line width for ERPs"),
    :color => (nothing, "Color for ERPs (single color or a vector of colors, one per dataset)"),
    :linestyle => (nothing, "Line style for ERPs (single style or a vector of styles, one per dataset)"),
    :colormap => (nothing, "Colormap for multi-channel plots"),

    # Plot configuration
    :average_channels => (false, "Average across channels (true/false)"),
    :error_bars => (:none, "Error bars type to plot: :none, :sem, :within_sem, :ci95"),

    # Legend parameters
    :legend_kwargs => (
        NamedTuple(),
        "Additional kwargs passed directly to the Makie Legend block (see Makie's [Legend documentation](https://docs.makie.org/stable/reference/blocks/legend/) for available attributes).",
    ),
    :legend => (true, "Show the legend (true/false)"),
    :legend_label => ("", "Title for the legend"),
    :legend_framevisible => (true, "Legend frame visible (true/false)"),
    :legend_position => (:lt, "Legend position (:lt, :rt, :lb, :rb, or tuple like (:left, :top), or (0.5, 0.5))"),
    :legend_channel => ([], "Which channel to put the legend on."),
    :legend_labels => ([], "Legend labels."),
    :legend_nbanks => (nothing, "Number of legend columns."),

    # Individual axes grids
    :xgrid => (false, "Show x-axis grid (true/false)"),
    :ygrid => (false, "Show y-axis grid (true/false)"),
    :xminorgrid => (false, "Show x-axis minor grid (true/false)"),
    :yminorgrid => (false, "Show y-axis minor grid (true/false)"),

    # Origin lines and scale indicator
    :add_xy_origin => (true, "Add origin lines at x=0 and y=0 (true/false)"),
    :axis_type => (:standard, "Type of axis to draw (:standard or :origin)"),
    :scale_x_value => (nothing, "X-axis scale value/step size (e.g. 0.1 for 100 ms)"),
    :scale_y_value => (nothing, "Y-axis scale value/step size (e.g. 5.0 for 5 μV)"),

    # Layout parameters
    :layout_kwargs => (
        NamedTuple(),
        "Additional parameters to configure grid or topographical layouts (e.g., grid_dims, topo_plot_width, grid_rowgap)",
    ),
    :figure_padding => ((10, 30, 10, 10), "Padding around entire figure as (left, right, bottom, top) tuple (in pixels)"),

    # Highlight regions
    :highlight_regions => (
        nothing,
        "Highlight regions as a NamedTuple or Vector of NamedTuples. Each region: (x1, x2, y1=-Inf, y2=Inf, color=:gray, alpha=0.3)",
    ),
)

######################################
######### Internal Functions #########
######################################

function _estimate_dt(time_points::AbstractVector{<:Real})::Float64
    n = length(time_points)
    n <= 1 && return 0.001
    time_points isa AbstractRange && return Float64(step(time_points))
    @inbounds d1 = Float64(time_points[2] - time_points[1])
    tol = 1e-9 * max(1.0, abs(d1))
    @inbounds for i = 2:(n-1)
        if abs(Float64(time_points[i+1] - time_points[i]) - d1) > tol
            return Float64(median(diff(time_points)))
        end
    end
    return d1
end

function _find_continuous_regions(mask::AbstractVector{Bool}, time_points::AbstractVector{<:Real})
    regions = Tuple{Float64,Float64}[]
    (isempty(mask) || !any(mask)) && return regions

    half_dt = _estimate_dt(time_points) * 0.5
    in_region = false
    start_idx = 0

    @inbounds for i in eachindex(mask)
        is_sig = mask[i]
        if is_sig && !in_region
            in_region = true
            start_idx = i
        elseif !is_sig && in_region
            in_region = false
            push!(regions, (Float64(time_points[start_idx]) - half_dt, Float64(time_points[i-1]) + half_dt))
        end
    end
    if in_region
        push!(regions, (Float64(time_points[start_idx]) - half_dt, Float64(time_points[end]) + half_dt))
    end
    return regions
end

function _parse_positive_cm(tb::Textbox, fallback::Float64)::Float64
    val = tryparse(Float64, tb.displayed_string[])
    (isnothing(val) || val <= 0.0) && (val = tryparse(Float64, tb.stored_string[]))
    return (isnothing(val) || val <= 0.0) ? fallback : val
end

function _parse_menu_dpi(scale_menu::Menu)::Int
    scale_str = scale_menu.selection[]
    isnothing(scale_str) && return 72
    parsed = tryparse(Int, first(eachsplit(scale_str)))
    return something(parsed, 72)
end

function _set_textbox_str!(tb::Textbox, s::String)
    tb.displayed_string[] = s
    tb.stored_string[] = s
end

function _open_save_settings_dialog(fig)
    fig_settings = _eegfun_figure(size = (550, 350), figure_padding = 20)
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

    _set_textbox_str!(w_box, string(w_orig_cm))
    _set_textbox_str!(h_box, string(h_orig_cm))

    Label(fig_settings[5, 1], "Output Pixels:", halign = :left, font = :bold)
    pixel_label = Label(fig_settings[5, 2], "", halign = :left)

    Label(fig_settings[6, 1], "Transparent Background:", halign = :left, font = :bold)
    transparent_toggle = Toggle(fig_settings[6, 2], active = false, halign = :left)

    export_button = Button(fig_settings[7, 1:2], label = "Select File & Save...", font = :bold)

    function update_pixel_display()
        if format_menu.selection[] in ("PDF", "SVG")
            pixel_label.text[] = "N/A (Vector Format)"
            return
        end
        w_cm = _parse_positive_cm(w_box, w_orig_cm)
        h_cm = _parse_positive_cm(h_box, h_orig_cm)
        dpi_val = _parse_menu_dpi(scale_menu)
        pixel_label.text[] = "$(round(Int, (w_cm / 2.54) * dpi_val)) × $(round(Int, (h_cm / 2.54) * dpi_val)) px"
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
                (isnothing(h_curr) || abs(h_curr - h_new) > 0.01) && _set_textbox_str!(h_box, string(h_new))
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
                (isnothing(w_curr) || abs(w_curr - w_new) > 0.01) && _set_textbox_str!(w_box, string(w_new))
            end
            update_pixel_display()
        end
    end

    on(lock_toggle.active) do locked
        if locked
            w_val = tryparse(Float64, w_box.stored_string[])
            if !isnothing(w_val) && w_val > 0.0
                _set_textbox_str!(h_box, string(round(w_val * (h_orig / w_orig), digits = 2)))
                update_pixel_display()
            end
        end
    end

    update_pixel_display()

    on(export_button.clicks) do _
        fmt = format_menu.selection[]
        trans = transparent_toggle.active[]
        w_target_pt = (_parse_positive_cm(w_box, w_orig_cm) / 2.54) * 72.0
        h_target_pt = (_parse_positive_cm(h_box, h_orig_cm) / 2.54) * 72.0
        dpi_val = _parse_menu_dpi(scale_menu)
        filter_ext = "." * lowercase(fmt)

        @async begin
            try
                screen = Makie.getscreen(fig_settings.scene)
                !isnothing(screen) && close(screen)

                path = NativeFileDialog.save_file(pwd(); filterlist = lowercase(fmt))
                if !isempty(path)
                    endswith(lowercase(path), filter_ext) || (path *= filter_ext)

                    orig_bg = fig.scene.backgroundcolor[]
                    orig_content_bgs = Dict{Any,Any}()
                    if trans
                        fig.scene.backgroundcolor[] = RGBAf(0, 0, 0, 0)
                        for content in fig.content
                            if hasproperty(content, :backgroundcolor) && content.backgroundcolor isa Observable
                                orig_content_bgs[content] = content.backgroundcolor[]
                                content.backgroundcolor[] = RGBAf(0, 0, 0, 0)
                            end
                        end
                    end

                    w_orig_sz, h_orig_sz = size(fig.scene)
                    resize!(fig.scene, (w_target_pt, h_target_pt))
                    try
                        if filter_ext in (".pdf", ".svg")
                            _save_vector_figure(path, fig; pt_per_unit = 1.0)
                        else
                            Makie.save(path, fig; px_per_unit = dpi_val / 72.0)
                        end
                        @info "Figure saved to: $path"
                    finally
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
                @minimal_warning "Failed to save figure: $err"
            end
        end
    end

    _display_in_screen(fig_settings)
end

_create_screen(; size = nothing) = MAKIE_EXT_STATE.create_screen(; size = size)

_display_in_screen(fig; size = nothing) = MAKIE_EXT_STATE.display_screen(fig; size = size)

_save_vector_figure(path, fig; pt_per_unit = 1.0) = MAKIE_EXT_STATE.save_vector(path, fig; pt_per_unit = pt_per_unit)

_is_glmakie_available() = MAKIE_EXT_STATE.glmakie_active

_set_window_title(title::String) = MAKIE_EXT_STATE.set_title(title)

function _check_interactive_demo_backend(demo_name::String = "Interactive demo")
    if isdefined(Main, :GLMakie)
        try
            Main.GLMakie.activate!(inline = false)
        catch
        end
    end
    backend = Makie.current_backend()
    if backend === missing || string(backend) == "CairoMakie" || (!_is_glmakie_available() && !isdefined(Main, :GLMakie))
        @minimal_warning """
        $demo_name requires GLMakie for live sliders, buttons, and animations.
        Currently active backend: $(backend)
        Please activate GLMakie in your session:
            using GLMakie
            GLMakie.activate!(inline = false)
        """
        return false
    end
    return true
end

function _display_demo_figure(fig)
    isnothing(fig) && return nothing
    backend = Makie.current_backend()
    (backend === missing || isnothing(backend) || string(backend) == "CairoMakie") && return nothing
    try
        display(fig)
    catch e
        @debug "Could not display demo figure" exception = e
    end
    return nothing
end

_has_focused_textbox(fig::Figure)::Bool = any(block -> block isa Textbox && block.focused[], fig.content)

function _display_figure(fig)
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press && event.key == Keyboard.s
            kb_state = events(fig).keyboardstate
            if (Keyboard.left_control in kb_state || Keyboard.right_control in kb_state) && !_has_focused_textbox(fig)
                _open_save_settings_dialog(fig)
            end
        end
    end
    _display_in_screen(fig)
end

_resolve_from_scene(scene, attr::Symbol, user_val, default_val) =
    isnothing(user_val) ? (haskey(Makie.theme(scene), attr) ? Makie.theme(scene)[attr][] : default_val) : user_val

_resolve_theme_colormap(obj, user_cmap, default_cmap = DEFAULT_COLORMAP) =
    _resolve_from_scene(obj.scene, :colormap, user_cmap, default_cmap)

_resolve_theme_linewidth(obj, user_lw, default_lw = 2) = _resolve_from_scene(obj.scene, :linewidth, user_lw, default_lw)

_resolve_theme_textcolor(obj, user_color, default_color = :black) = _resolve_from_scene(obj.scene, :textcolor, user_color, default_color)

_resolve_theme_linecolor(obj, user_color, default_color = :black) = _resolve_from_scene(obj.scene, :linecolor, user_color, default_color)

function _is_diverging_colormap(cmap)::Bool
    cmap isa Symbol && return cmap in _DIVERGING_COLORMAPS
    cmap isa AbstractString && return Symbol(cmap) in _DIVERGING_COLORMAPS
    cmap isa Makie.Reverse && return _is_diverging_colormap(cmap.data)
    return false
end

function _symmetric_colorrange(min_val::Real, max_val::Real; fallback::Float64 = 1.0)::Tuple{Float64,Float64}
    max_abs = Float64(max(abs(min_val), abs(max_val)))
    (!isfinite(max_abs) || max_abs == 0.0) && return (-fallback, fallback)
    return (-max_abs, max_abs)
end

@inline _kwargs_dict(kwargs)::Dict{Symbol,Any} = Dict{Symbol,Any}(Symbol(k) => v for (k, v) in pairs(kwargs))

function _get_colorbar_position(pos, row_span, col_span)
    pos isa Tuple && return pos
    pos == :right && return (row_span, last(col_span) + 1)
    pos == :left && return (row_span, first(col_span) - 1)
    pos == :bottom && return (last(row_span) + 1, col_span)
    pos == :top && return (first(row_span) - 1, col_span)
    return (row_span, last(col_span) + 1)
end

function _set_axis_grid!(ax; xgrid = false, ygrid = false, xminorgrid = false, yminorgrid = false)
    ax.xgridvisible = xgrid
    ax.ygridvisible = ygrid
    ax.xminorgridvisible = xminorgrid
    ax.yminorgridvisible = yminorgrid
end

function _make_step_ticks(step_val)
    return (vmin, vmax) -> begin
        start_v = ceil(vmin / step_val) * step_val
        end_v = floor(vmax / step_val) * step_val
        ticks = collect(start_v:step_val:(end_v+1e-9))
        return ticks, string.(round.(ticks, sigdigits = 4))
    end
end

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
    ax.xlabel = xlabel
    ax.ylabel = ylabel
    ax.yreversed = yreversed
    !isnothing(scale_x_value) && (ax.xticks = _make_step_ticks(scale_x_value))
    !isnothing(scale_y_value) && (ax.yticks = _make_step_ticks(scale_y_value))
    !isnothing(xlim) && xlims!(ax, xlim[1], xlim[2])
    !isnothing(ylim) && ylims!(ax, ylim[1], ylim[2])
end

function _set_origin_lines!(ax; add_xy_origin = true, color = :gray, linewidth = 1, alpha = 0.8)
    if add_xy_origin
        hlines!(ax, 0, color = color, linewidth = linewidth, alpha = alpha)
        vlines!(ax, 0, color = color, linewidth = linewidth, alpha = alpha)
    end
end

function _add_origin_scale_indicator!(
    ax;
    axis_type = :standard,
    scale_x_value = nothing,
    scale_y_value = nothing,
    xlabel = "s",
    ylabel = "μV",
)
    axis_type == :origin || return
    hidedecorations!(ax)
    hidespines!(ax)

    linecolor = _resolve_theme_linecolor(ax, nothing)
    textcolor = _resolve_theme_textcolor(ax, nothing)
    hlines!(ax, 0, color = linecolor, linewidth = 1, overdraw = true)
    vlines!(ax, 0, color = linecolor, linewidth = 1, overdraw = true)

    x_ticks_obs = lift(ax.finallimits) do lims
        min_x, max_x = lims.origin[1], lims.origin[1] + lims.widths[1]
        collect((ceil(min_x/scale_x_value)*scale_x_value):scale_x_value:(floor(max_x/scale_x_value)*scale_x_value+1e-9))
    end
    y_ticks_obs = lift(ax.finallimits) do lims
        min_y, max_y = lims.origin[2], lims.origin[2] + lims.widths[2]
        collect((ceil(min_y/scale_y_value)*scale_y_value):scale_y_value:(floor(max_y/scale_y_value)*scale_y_value+1e-9))
    end

    x_pts_obs = lift(ticks -> [Point2f(t, 0.0) for t in ticks if abs(t) > 1e-9], x_ticks_obs)
    y_pts_obs = lift(ticks -> [Point2f(0.0, t) for t in ticks if abs(t) > 1e-9], y_ticks_obs)
    x_lbls = lift(ticks -> [string(round(t, sigdigits = 3)) for t in ticks if abs(t) > 1e-9], x_ticks_obs)
    y_lbls = lift(ticks -> [string(round(t, sigdigits = 3)) for t in ticks if abs(t) > 1e-9], y_ticks_obs)

    scatter!(ax, x_pts_obs, marker = '|', markersize = 10, color = linecolor, overdraw = true, xautolimits = false, yautolimits = false)
    scatter!(ax, y_pts_obs, marker = '-', markersize = 10, color = linecolor, overdraw = true, xautolimits = false, yautolimits = false)

    text!(
        ax,
        x_pts_obs,
        text = x_lbls,
        align = (:center, :top),
        offset = (0, -8),
        color = textcolor,
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
        color = textcolor,
        fontsize = 12,
        overdraw = true,
        xautolimits = false,
        yautolimits = false,
    )

    x_axis_label_pos = lift(lims -> Point2f(lims.origin[1] + lims.widths[1], 0.0), ax.finallimits)
    y_axis_label_pos = lift(lims -> Point2f(0.0, lims.origin[2] + lims.widths[2]), ax.finallimits)
    text!(
        ax,
        x_axis_label_pos,
        text = string(" ", xlabel),
        align = (:left, :center),
        color = textcolor,
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
        color = textcolor,
        fontsize = 14,
        overdraw = true,
        xautolimits = false,
        yautolimits = false,
    )
end

function _split_into_parts(s::String)
    parts = String[]
    isempty(s) && return parts
    for m in eachmatch(_SPLIT_PARTS_RE, s)
        if !isnothing(m.captures[3])
            push!(parts, "_")
        elseif !isnothing(m.captures[4])
            push!(parts, String(m.captures[4]))
        else
            push!(parts, String(m.captures[1]))
            digits = m.captures[2]
            !isempty(digits) && push!(parts, String(digits))
        end
    end
    return parts
end

function _generate_window_title(datasets::Vector{<:EegData}; max_total_length::Int = 80, max_name_length::Int = 20)
    isempty(datasets) && return ""
    length(datasets) == 1 && return "$(datasets[1].file):$(datasets[1].condition_name)"

    first_file = datasets[1].file
    if all(d -> d.file == first_file, datasets)
        condition_names = [d.condition_name for d in datasets]
        condition_str = _shorten_condition_names(
            condition_names;
            max_total_length = max_total_length - length(first_file) - 2,
            max_name_length = max_name_length,
        )
        return "$first_file: $condition_str"
    end
    file_condition_pairs = ["$(d.file):$(d.condition_name)" for d in datasets]
    return _shorten_condition_names(file_condition_pairs; max_total_length = max_total_length, max_name_length = max_name_length)
end

_generate_window_title(datasets::EegData) = "$(datasets.file):$(datasets.condition_name)"

_generate_window_title(datasets::ContinuousData) = datasets.file

function _generate_window_title(datasets::Vector{ContinuousData}; max_total_length::Int = 80, max_name_length::Int = 20)
    isempty(datasets) && return ""
    length(datasets) == 1 && return datasets[1].file
    first_file = datasets[1].file
    all(d -> d.file == first_file, datasets) && return first_file
    return _shorten_condition_names([d.file for d in datasets]; max_total_length = max_total_length, max_name_length = max_name_length)
end

function _write_abbreviated_part!(io::IOBuffer, part::String)
    isempty(part) && return
    if part == "_" || isdigit(part[1])
        write(io, part)
    else
        write(io, uppercase(part[1]))
        n = min(3, length(part))
        n >= 2 && write(io, SubString(part, 2, n))
    end
end

function _abbreviate_parts(parts::Vector{String})
    io = IOBuffer()
    for part in parts
        _write_abbreviated_part!(io, part)
    end
    return String(take!(io))
end

function _abbreviate_parts_with_prefix(name_parts::Vector{String}, n_prefix::Int)
    io = IOBuffer()
    for i = 1:min(n_prefix, length(name_parts))
        _write_abbreviated_part!(io, name_parts[i])
    end
    for i = (n_prefix+1):length(name_parts)
        write(io, name_parts[i])
    end
    return String(take!(io))
end

function _common_prefix_from_parts(all_parts::Vector{Vector{String}})
    length(all_parts) <= 1 && return String[]
    first_parts = all_parts[1]
    n_common = 0
    for (i, part) in enumerate(first_parts)
        if all(parts -> i <= length(parts) && parts[i] == part, all_parts)
            n_common = i
        else
            break
        end
    end
    return n_common >= 2 ? first_parts[1:n_common] : String[]
end

function _shorten_condition_names(
    condition_names::Vector{String};
    max_total_length::Int = 80,
    max_name_length::Int = 30,
    separator::String = ", ",
    show_ends::Int = 3,
)
    isempty(condition_names) && return ""
    length(condition_names) == 1 && return condition_names[1]

    all_parts = [_split_into_parts(name) for name in condition_names]
    common_prefix = _common_prefix_from_parts(all_parts)
    n_prefix = length(common_prefix)
    abbreviated_names =
        n_prefix == 0 ? [_abbreviate_parts(parts) for parts in all_parts] :
        [_abbreviate_parts_with_prefix(parts, n_prefix) for parts in all_parts]

    orig_len = sum(length, condition_names)
    abbrev_len = sum(length, abbreviated_names)
    abbrev_len >= orig_len * 0.8 && (abbreviated_names = condition_names)

    final_names = [length(n) > max_name_length ? first(n, max_name_length) * "…" : n for n in abbreviated_names]
    full_string = join(final_names, separator)

    if length(full_string) > max_total_length && length(condition_names) > 2 * show_ends
        return join(@view(final_names[1:show_ends]), separator) *
               separator *
               "…" *
               separator *
               join(@view(final_names[(end-show_ends+1):end]), separator)
    end
    return full_string
end

_extract_baseline_values(::Nothing) = (nothing, nothing)

_extract_baseline_values(interval::AbstractRange) = (first(interval), last(interval))

_extract_baseline_values(interval::Tuple) = (interval[1], interval[2])

_extract_baseline_values(interval::AbstractVector{<:Real}) = (Float64(interval[1]), Float64(interval[2]))

_extract_baseline_values(interval::TimeSelection) = (Float64(interval.start), Float64(interval.stop))

function _parse_baseline_values(start_str::String, stop_str::String)
    (start_str == " " || stop_str == " ") && return -Inf, Inf
    start_val = tryparse(Float64, start_str)
    stop_val = tryparse(Float64, stop_str)
    if isnothing(start_val) || isnothing(stop_val)
        @minimal_warning "Invalid baseline values: \"$start_str\", \"$stop_str\""
        return nothing, nothing
    end
    return start_val, stop_val
end

function _create_baseline_textbox(layout, row, label, obs, placeholder, width)
    Label(layout[row, 1], label, width = 60)
    tb = Textbox(layout[row, 2], placeholder = placeholder, width = width)
    tb.stored_string[] = obs[]
    hasproperty(tb, :displayed_string) ? connect!(obs, tb.displayed_string) : connect!(obs, tb.stored_string)
    return tb
end
