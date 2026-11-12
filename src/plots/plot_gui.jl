const UI_STYLE = (label_font = 20, button_font = 18, textbox_font = 16, input_width = 160, input_height = 30)

Base.@kwdef mutable struct GUIState
    directory::Observable{String} = Observable("")
    filename::Observable{String} = Observable("")
    participant::Observable{String} = Observable("")
    condition::Observable{String} = Observable("")
    epoch::Observable{String} = Observable("")
    plottype::Observable{String} = Observable("select")
    layout::Observable{String} = Observable("select")
    layout_file::Observable{String} = Observable("")
    layout_object::Observable{Union{Nothing,Layout}} = Observable{Union{Nothing,Layout}}(nothing)
    electrodes::Observable{Vector{String}} = Observable(String[])
    xlim::Observable{Tuple{Union{Nothing,Float64},Union{Nothing,Float64}}} =
        Observable{Tuple{Union{Nothing,Float64},Union{Nothing,Float64}}}((nothing, nothing))
    ylim::Observable{Tuple{Union{Nothing,Float64},Union{Nothing,Float64}}} =
        Observable{Tuple{Union{Nothing,Float64},Union{Nothing,Float64}}}((nothing, nothing))
    zlim::Observable{Tuple{Union{Nothing,Float64},Union{Nothing,Float64}}} =
        Observable{Tuple{Union{Nothing,Float64},Union{Nothing,Float64}}}((nothing, nothing))
    baseline_start::Observable{Union{Nothing,Float64}} = Observable{Union{Nothing,Float64}}(nothing)
    baseline_end::Observable{Union{Nothing,Float64}} = Observable{Union{Nothing,Float64}}(nothing)
    baseline_type::Observable{String} = Observable("select")
    layout_type::Observable{String} = Observable("single")
    average_channels::Observable{Bool} = Observable(false)
    invert_y::Observable{Bool} = Observable(false)
    submenu_active::Observable{Bool} = Observable(false)
    submenu_type::Observable{String} = Observable("")
    channel_menu::Union{Nothing,Menu} = nothing
    ica_components::Observable{String} = Observable("")
end

_create_select_button(parent, label) = Button(
    parent,
    label = label,
    fontsize = UI_STYLE.button_font,
    width = UI_STYLE.input_width,
    height = UI_STYLE.input_height,
    buttoncolor = :lightgrey,
    buttoncolor_hover = :grey,
    buttoncolor_active = :green,
)

function _create_textbox(parent; width = nothing, placeholder = "", boxcolor_obs = nothing, halign = :center)
    Textbox(
        parent,
        placeholder = placeholder,
        fontsize = UI_STYLE.textbox_font,
        width = isnothing(width) ? UI_STYLE.input_width : width,
        height = UI_STYLE.input_height,
        halign = halign,
        boxcolor = isnothing(boxcolor_obs) ? Observable(:white) : boxcolor_obs,
        cornerradius = 0,
    )
end

_create_menu(parent; options = ["Select"], width = nothing, height = nothing) = Menu(
    parent,
    options = options,
    width = isnothing(width) ? UI_STYLE.input_width : width,
    height = isnothing(height) ? UI_STYLE.input_height : height,
    fontsize = UI_STYLE.textbox_font,
)

_create_label(parent, text; fontsize = UI_STYLE.label_font, color = :black) = Label(parent, text, fontsize = fontsize, color = color)

_truncate_path(path::String, max_length::Int = 20) = length(path) <= max_length ? path : "..." * path[(end-(max_length-4)):end]

"""
    plot_gui(; display_plot::Bool = true)

Launch the interactive GUI for quick EEG data exploration and visualization.

Provides file/directory browsing, plot-type selection (Data Browser, Epochs, ERP,
Topography, ERP Image, Time-Frequency, ICA, ERP Analysis, and Diagnostic plots),
participant/condition/epoch/channel filtering, and axis/baseline controls.

# Examples
```julia
EegFun.plot_gui()
```
"""
function plot_gui(; display_plot::Bool = true)
    _set_window_title("PLOT GUI")
    gui_fig = _eegfun_figure(size = (700, 600), title = "Plot GUI", backgroundcolor = :lightgrey, figure_padding = (20, 20, 20, 20))
    main_layout = GridLayout(gui_fig[1, 1:3], rowgap = 2, colgap = 4)

    channel_menu = _create_menu(main_layout[10, 2], options = ["Select"])
    gui_state = GUIState(; channel_menu = channel_menu)

    _build_file_column!(main_layout, gui_state)
    _build_filter_column!(main_layout, gui_state)
    _build_axis_column!(main_layout, gui_state)

    colsize!(main_layout, 1, Fixed(UI_STYLE.input_width))
    colsize!(main_layout, 2, Auto())
    colsize!(main_layout, 3, Auto())

    screen = display_plot ? _display_in_screen(gui_fig) : nothing
    _set_window_title("Makie")
    return (fig = gui_fig, axes = Axis[], screen = screen)
end

function _build_file_column!(main_layout, gui_state)
    directory_select_button = _create_select_button(main_layout[1, 1], "Select Directory")
    directory_label_text = Observable("Dir:")
    _create_label(main_layout[2, 1], directory_label_text, fontsize = UI_STYLE.textbox_font, color = :gray)

    file_select_button = _create_select_button(main_layout[3, 1], "Select File")
    file_pattern_input = _create_textbox(main_layout[4, 1], placeholder = "", halign = :left)

    layout_select_button = _create_select_button(main_layout[5, 1], "Select Layout")
    layout_label_text = Observable("File:")
    _create_label(main_layout[6, 1], layout_label_text, fontsize = UI_STYLE.textbox_font, color = :gray)

    _create_label(main_layout[7, 1], "Plot Type")
    plottype_options = [
        "Select",
        "Data Browser",
        "Epochs",
        "ERP",
        "Topography",
        "ERP Image",
        "Time-Frequency",
        "─────────────────",
        "ICA >",
        "ERP Analysis >",
        "Diagnostic Plots >",
    ]
    plottype_dropdown = _create_menu(main_layout[8, 1], options = plottype_options)
    submenu_dropdown = _create_menu(main_layout[9, 1], options = ["---"])

    component_bc = Observable(:lightgrey)
    component_input = _create_textbox(main_layout[10, 1], placeholder = "", boxcolor_obs = component_bc)

    on(directory_select_button.clicks) do _
        dir_path = pick_folder(_gui_dir(gui_state))
        if !isnothing(dir_path) && dir_path != ""
            gui_state.directory[] = dir_path
            directory_label_text[] = "." * _truncate_path(String(strip(dir_path)))
        end
    end

    on(file_select_button.clicks) do _
        filename = pick_file(_gui_dir(gui_state))
        if !isnothing(filename) && filename != ""
            file_pattern_input.displayed_string[] = ""
            file_pattern_input.displayed_string[] = basename(filename)
            gui_state.filename[] = filename
        end
    end

    _bind_textbox_validation!(file_pattern_input) do val
        if !isempty(val) && val != basename(gui_state.filename[])
            gui_state.filename[] = val
        end
    end

    on(layout_select_button.clicks) do _
        filename = pick_file(_gui_dir(gui_state))
        if !isnothing(filename) && filename != ""
            basename_only = basename(filename)
            gui_state.layout_file[] = filename
            gui_state.layout[] = basename_only
            layout_label_text[] = basename_only
            try
                layout = read_layout(filename)
                gui_state.layout_object[] = layout
                gui_state.channel_menu.options = vcat(["Select"], string.(channel_labels(layout)))
                @info "Loaded layout with $(length(channel_labels(layout))) channels"
            catch layout_error
                @minimal_warning "Error loading layout: $layout_error"
                gui_state.layout_object[] = nothing
            end
        end
    end

    updating_submenu = Ref(false)
    on(plottype_dropdown.selection) do selection
        (selection == "Select" || selection == "─────────────────") && return
        updating_submenu[] = true
        try
            if selection == "ICA >"
                submenu_dropdown.options = ["Select", "Components", "Activation"]
                gui_state.submenu_active[] = true
                gui_state.submenu_type[] = "ica"
                component_bc[] = :white
            elseif selection == "ERP Analysis >"
                submenu_dropdown.options = ["Select", "ERP GUI"]
                gui_state.submenu_active[] = true
                gui_state.submenu_type[] = "erp_analysis"
                component_bc[] = :lightgrey
            elseif selection == "Diagnostic Plots >"
                submenu_dropdown.options = ["Select", "Channel Summary", "Layout View"]
                gui_state.submenu_active[] = true
                gui_state.submenu_type[] = "diagnostic"
                component_bc[] = :lightgrey
            else
                gui_state.plottype[] = selection
                gui_state.submenu_active[] = false
                component_bc[] = :lightgrey
            end
        finally
            updating_submenu[] = false
        end
    end

    on(submenu_dropdown.selection) do selection
        (updating_submenu[] || isnothing(selection) || selection == "Select") && return
        gui_state.plottype[] = selection
    end

    _bind_textbox_validation!(component_input) do val
        gui_state.ica_components[] = val
    end
end

function _build_filter_column!(main_layout, gui_state)
    participant_bc = Observable(:white)
    condition_bc = Observable(:white)
    epoch_bc = Observable(:white)

    _create_label(main_layout[1, 2], "Participant")
    participant_input = _create_textbox(main_layout[2, 2], boxcolor_obs = participant_bc)
    _create_label(main_layout[3, 2], "Condition")
    condition_input = _create_textbox(main_layout[4, 2], boxcolor_obs = condition_bc)
    _create_label(main_layout[5, 2], "Epoch")
    epoch_input = _create_textbox(main_layout[6, 2], boxcolor_obs = epoch_bc)

    _create_label(main_layout[7, 2], "Layout")
    layout_dropdown = _create_menu(main_layout[8, 2], options = ["Single", "Single Avg", "Grid", "Topo"])

    _create_label(main_layout[9, 2], "Channel(s)")
    selected_channels_text = Observable("Channels: ")
    _create_label(main_layout[11, 2], selected_channels_text, fontsize = UI_STYLE.textbox_font, color = :gray)

    setup_multi_int_callback(participant_input, participant_bc, "Participant", gui_state.participant)
    setup_multi_int_callback(condition_input, condition_bc, "Condition", gui_state.condition)
    setup_multi_int_callback(epoch_input, epoch_bc, "Epoch", gui_state.epoch)

    on(layout_dropdown.selection) do selection
        isnothing(selection) && return
        if selection == "Single Avg"
            gui_state.layout_type[] = "single"
            gui_state.average_channels[] = true
        else
            gui_state.layout_type[] = lowercase(selection)
            gui_state.average_channels[] = false
        end
    end

    resetting_channel_menu = Ref(false)
    on(gui_state.channel_menu.selection) do selection
        (resetting_channel_menu[] || isnothing(selection)) && return
        if selection == "Select"
            gui_state.electrodes[] = String[]
        elseif selection in gui_state.electrodes[]
            gui_state.electrodes[] = filter(!=(selection), gui_state.electrodes[])
        else
            gui_state.electrodes[] = vcat(gui_state.electrodes[], selection)
        end
        selected_channels_text[] = isempty(gui_state.electrodes[]) ? "Selected: " : "Selected: " * join(gui_state.electrodes[], ", ")
        resetting_channel_menu[] = true
        try
            gui_state.channel_menu.i_selected[] = 0
        finally
            resetting_channel_menu[] = false
        end
    end
end

function _create_limit_pair!(main_layout, row, label_str)
    min_bc = Observable(:white)
    max_bc = Observable(:white)
    _create_label(main_layout[row, 3], label_str, fontsize = UI_STYLE.textbox_font)
    sub = GridLayout(main_layout[row+1, 3], tellwidth = false, colgap = 8)
    min_in = _create_textbox(sub[1, 1], width = 80, boxcolor_obs = min_bc)
    max_in = _create_textbox(sub[1, 2], width = 80, boxcolor_obs = max_bc)
    return min_in, max_in, min_bc, max_bc
end

function _build_axis_column!(main_layout, gui_state)
    xmin_input, xmax_input, xmin_bc, xmax_bc = _create_limit_pair!(main_layout, 1, "X Limits")
    ymin_input, ymax_input, ymin_bc, ymax_bc = _create_limit_pair!(main_layout, 3, "Y Limits")
    zmin_input, zmax_input, zmin_bc, zmax_bc = _create_limit_pair!(main_layout, 5, "Z Limits")
    baseline_start, baseline_end, bl_start_bc, bl_end_bc = _create_limit_pair!(main_layout, 7, "Baseline")

    _create_label(main_layout[9, 3], "Baseline Type TF", fontsize = UI_STYLE.textbox_font)
    baseline_type =
        _create_menu(main_layout[10, 3], options = ["Select", "db", "absolute", "relative", "relchange", "percent", "zscore", "normchange"])

    invert_y_layout = GridLayout(main_layout[11, 3], tellwidth = false, colgap = 8)
    _create_label(invert_y_layout[1, 1], "Invert Y axis", fontsize = UI_STYLE.textbox_font)
    invert_y_checkbox = Checkbox(invert_y_layout[1, 2], checked = false)

    plot_button = Button(
        main_layout[13, 3],
        label = "PLOT",
        fontsize = UI_STYLE.button_font,
        buttoncolor = :darkgrey,
        buttoncolor_hover = :grey,
        buttoncolor_active = :green,
        width = 100,
        height = 50,
    )

    x_check = () -> check_range(gui_state.xlim, xmin_bc, xmax_bc, "X Limits")
    y_check = () -> check_range(gui_state.ylim, ymin_bc, ymax_bc, "Y Limits")
    z_check = () -> check_range(gui_state.zlim, zmin_bc, zmax_bc, "Z Limits")
    bl_check = () -> begin
        s, e = gui_state.baseline_start[], gui_state.baseline_end[]
        if !isnothing(s) && !isnothing(e) && s > e
            bl_start_bc[] = :lightcoral
            bl_end_bc[] = :lightcoral
            @minimal_warning "Invalid Baseline range: start ($s) must be less than end ($e)"
        end
    end

    setup_numeric_callback(xmin_input, xmin_bc, "X min", v -> gui_state.xlim[] = (v, gui_state.xlim[][2]); range_check = x_check)
    setup_numeric_callback(xmax_input, xmax_bc, "X max", v -> gui_state.xlim[] = (gui_state.xlim[][1], v); range_check = x_check)
    setup_numeric_callback(ymin_input, ymin_bc, "Y min", v -> gui_state.ylim[] = (v, gui_state.ylim[][2]); range_check = y_check)
    setup_numeric_callback(ymax_input, ymax_bc, "Y max", v -> gui_state.ylim[] = (gui_state.ylim[][1], v); range_check = y_check)
    setup_numeric_callback(zmin_input, zmin_bc, "Z min", v -> gui_state.zlim[] = (v, gui_state.zlim[][2]); range_check = z_check)
    setup_numeric_callback(zmax_input, zmax_bc, "Z max", v -> gui_state.zlim[] = (gui_state.zlim[][1], v); range_check = z_check)
    setup_numeric_callback(baseline_start, bl_start_bc, "Baseline start", v -> gui_state.baseline_start[] = v; range_check = bl_check)
    setup_numeric_callback(baseline_end, bl_end_bc, "Baseline end", v -> gui_state.baseline_end[] = v; range_check = bl_check)

    on(baseline_type.selection) do selection
        !isnothing(selection) && selection != "Select" && (gui_state.baseline_type[] = selection)
    end
    on(invert_y_checkbox.checked) do is_checked
        gui_state.invert_y[] = is_checked
    end

    plot_dispatch = Dict{String,Function}(
        "Data Browser" => _plot_databrowser,
        "Epochs" => _plot_epochs,
        "ERP" => _plot_erp,
        "ERP Image" => _plot_erp_image,
        "Topography" => _plot_topography,
        "Time-Frequency" => _plot_time_frequency,
        "Components" => _plot_ica,
        "Activation" => _plot_ica_activation,
        "ERP GUI" => _plot_erp_measurement_gui,
        "Channel Summary" => _plot_channel_summary,
        "Layout View" => _plot_layout,
    )

    on(plot_button.clicks) do _
        plot_type = gui_state.plottype[]
        (plot_type == "Select" || startswith(plot_type, "─") || plot_type == "select") && return
        handler = get(plot_dispatch, plot_type, nothing)
        isnothing(handler) ? @minimal_warning("Unsupported plot type: $plot_type") : handler(gui_state)
    end
end

function _plot_databrowser(gui_state)
    fname = gui_state.filename[]
    file_ext = lowercase(splitext(fname)[2])
    if fname == "" || file_ext ∉ (".bdf", ".jld2")
        @minimal_warning "Requested plot settings incompatible: recheck!"
        return
    end

    _guarded("Data Browser") do
        if file_ext == ".jld2"
            _guarded_async("Data Browser") do
                dat = read_data(fname)
                isnothing(dat) && return
                dat = _apply_gui_filters(dat, gui_state)
                isnothing(dat) && return
                if dat isa Vector
                    plot_databrowser(dat)
                else
                    plot_databrowser(dat; screen = _create_screen())
                end
            end
        else
            layout = nothing
            if gui_state.layout_file[] != ""
                layout = read_layout(gui_state.layout_file[])
                polar_to_cartesian_xy!(layout)
            end
            raw = read_raw_data(fname)
            dat = isnothing(layout) ? create_eegfun_data(raw) : create_eegfun_data(raw, layout)
            gui_state.channel_menu.options = vcat(["Select"], string.(channel_labels(dat)))
            _guarded_async("Data Browser") do
                plot_databrowser(dat; screen = _create_screen())
            end
        end
    end
end

function _plot_epochs(gui_state)
    isnothing(_validate_file(gui_state)) && return
    fname = gui_state.filename[]
    cond_list = _parse_gui_int_field(gui_state.condition[])
    if isnothing(cond_list)
        @minimal_warning "Requested plot settings incompatible: recheck!"
        return
    end
    cond_sel = isempty(cond_list) ? conditions() : conditions(cond_list)
    selected_channels = _gui_selected_channels(gui_state)

    epoch_list = _parse_gui_int_field(gui_state.epoch[])
    if isnothing(epoch_list)
        @minimal_warning "Requested plot settings incompatible: recheck!"
        return
    end
    epoch_sel = isempty(epoch_list) ? epochs() : epochs(epoch_list)

    if _is_pattern(fname)
        res = _resolve_batch_files_helper(gui_state)
        isnothing(res) && return
        _, part_sel, _ = res
        _guarded_async("Epochs") do
            plot_epochs(
                fname;
                input_dir = _gui_dir(gui_state),
                participant_selection = part_sel,
                condition_selection = cond_sel,
                channel_selection = selected_channels,
                epoch_selection = epoch_sel,
                layout = Symbol(gui_state.layout_type[]),
                average_channels = gui_state.average_channels[],
                xlim = _gui_lim(gui_state.xlim[]),
                ylim = _gui_lim(gui_state.ylim[]),
                yreversed = gui_state.invert_y[],
            )
        end
        return
    end

    _guarded_async("Epochs") do
        plot_epochs(
            fname;
            channel_selection = selected_channels,
            condition_selection = cond_sel,
            epoch_selection = epoch_sel,
            layout = Symbol(gui_state.layout_type[]),
            average_channels = gui_state.average_channels[],
            xlim = _gui_lim(gui_state.xlim[]),
            ylim = _gui_lim(gui_state.ylim[]),
            yreversed = gui_state.invert_y[],
        )
    end
end

function _plot_erp(gui_state)
    isnothing(_validate_file(gui_state)) && return
    fname = gui_state.filename[]
    cond_list = _parse_gui_int_field(gui_state.condition[])
    if isnothing(cond_list)
        @minimal_warning "Requested plot settings incompatible: recheck!"
        return
    end
    cond_sel = isempty(cond_list) ? conditions() : conditions(cond_list)
    selected_channels = _gui_selected_channels(gui_state)
    baseline = _gui_lim((gui_state.baseline_start[], gui_state.baseline_end[]))

    if _is_pattern(fname)
        res = _resolve_batch_files_helper(gui_state)
        isnothing(res) && return
        _, part_sel, _ = res
        input_dir = _gui_dir(gui_state)
        _guarded_async("ERP") do
            plot_erp(
                fname;
                input_dir = input_dir,
                participant_selection = part_sel,
                layout = Symbol(gui_state.layout_type[]),
                channel_selection = selected_channels,
                condition_selection = cond_sel,
                baseline_interval = baseline,
                average_channels = gui_state.average_channels[],
                xlim = _gui_lim(gui_state.xlim[]),
                ylim = _gui_lim(gui_state.ylim[]),
                yreversed = gui_state.invert_y[],
            )
        end
        return
    end

    _guarded_async("ERP") do
        plot_erp(
            fname;
            layout = Symbol(gui_state.layout_type[]),
            condition_selection = cond_sel,
            channel_selection = selected_channels,
            baseline_interval = baseline,
            average_channels = gui_state.average_channels[],
            xlim = _gui_lim(gui_state.xlim[]),
            ylim = _gui_lim(gui_state.ylim[]),
            yreversed = gui_state.invert_y[],
        )
    end
end

function _plot_topography(gui_state)
    isnothing(_validate_file(gui_state)) && return
    fname = gui_state.filename[]
    selected_channels = _gui_selected_channels(gui_state)

    function _do_topo_plot(data)
        data = _apply_gui_filters(data, gui_state)
        isnothing(data) && return

        baseline_interval = _gui_lim((gui_state.baseline_start[], gui_state.baseline_end[]))
        !isnothing(baseline_interval) && baseline!(data, baseline_interval)

        x_lim = _gui_lim(gui_state.xlim[])
        interval_sel = !isnothing(x_lim) ? times(x_lim[1], x_lim[2]) : times()
        z_lim = _gui_lim(gui_state.zlim[])
        y_lim = _gui_lim(gui_state.ylim[])
        color_range = !isnothing(z_lim) ? z_lim : y_lim
        extra_kwargs = isnothing(color_range) ? (;) : (colorrange = color_range,)

        if data isa Vector{<:ErpData} || data isa ErpData
            plot_topography(data; channel_selection = selected_channels, interval_selection = interval_sel, extra_kwargs...)
        elseif data isa Vector{<:EpochData} || data isa EpochData
            epoch_str = strip(gui_state.epoch[])
            epoch_num = isempty(epoch_str) ? 1 : something(tryparse(Int, epoch_str), 1)
            plot_topography(data, epoch_num; channel_selection = selected_channels, interval_selection = interval_sel, extra_kwargs...)
        else
            @minimal_warning "Requested plot settings incompatible: recheck!"
        end
    end

    if _is_pattern(fname)
        if isnothing(_parse_gui_int_field(gui_state.condition[]))
            @minimal_warning "Requested plot settings incompatible: recheck!"
            return
        end
        res = _resolve_batch_files_helper(gui_state)
        isnothing(res) && return
        files, _, input_dir = res

        _guarded_async("Topography") do
            for file in sort(files, by = _natural_sort_key)
                @info "Plotting: $file"
                data = read_data(joinpath(input_dir, file))
                !isnothing(data) && _do_topo_plot(data)
            end
        end
        return
    end

    _guarded_async("Topography") do
        data = read_data(fname)
        isnothing(data) ? @minimal_warning("Requested plot settings incompatible: recheck!") : _do_topo_plot(data)
    end
end
