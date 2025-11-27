# =============================================================================
# NOTE ON SPELLING: COLOR vs COLOUR
# =============================================================================
# This codebase uses "color" (American spelling) throughout to maintain
# consistency with Makie.jl's API, which uses "color" in all its attributes
# (e.g., colormap, colorbar, strokecolor, etc.). While "colour" may be the
# preferred spelling in some regions, we follow Makie's convention for
# consistency and to avoid confusion when working with Makie's API.
# =============================================================================

function display_figure(fig)
    display(get_makie_screen(get_makie_backend()), fig)
end

# TODO: I must be doing something wrong here!!! 
# Must be a better way to do this?
function get_makie_backend()
    backend_str = string(Makie.current_backend())
    if occursin("GLMakie", backend_str)
        return :GLMakie
    elseif occursin("CairoMakie", backend_str)
        return :CairoMakie
    else
        @minimal_error "Makie backend not found"
    end
end

get_makie_screen(makie_backend::Symbol) = getfield(Main, makie_backend).Screen()

"""
    set_window_title(title::String)

Set the window title for GLMakie. Does nothing for CairoMakie (which doesn't support window titles).
"""
function set_window_title(title::String)
    backend = get_makie_backend()
    if backend == :GLMakie
        Makie.current_backend().activate!(title = title)
    end
    # TODO: CairoMakie doesn't support window titles?
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
    set_window_title(title)
    fig = Figure(; figure_kwargs...)
    ax = Axis(fig[1, 1])
    return fig, ax
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
            if value !== nothing  # Only add if not the default nothing
                colorbar_kwargs[attr] = value
            end
        end
    end

    # TODO: what did I need to add this?
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
        # attr in positioning_attrs && continue  # Skip positioning attributes if requested
        legend_key = Symbol("legend_$(attr)")
        if haskey(plot_kwargs, legend_key)
            value = plot_kwargs[legend_key]  # Read from plot_kwargs
            if value in [:legend_halign, :legend_valign, :legend_alignmode]
                println("Skipping positioning attribute: $attr")
                continue
            end
            if value !== nothing  # Only add if not the default nothing
                legend_kwargs[attr] = value  # Store in new dict without "legend_" prefix
            end
        end
    end
  
    # TODO: I do not know why this is needed? Bug here? Bug Makie axislegend?
    # Without it legend_position is ignored!
    # Remove positioning attributes that conflict with explicit position parameter
    [pop!(legend_kwargs, attr, nothing) for attr in [:halign, :valign, :alignmode]]

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
            if value !== nothing
                layout_kwargs[param_name] = value
            end
        end
    end
    
    return layout_kwargs
end

# =============================================================================
# AXIS STYLING FUNCTIONS
# =============================================================================

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
function _set_axis_properties!(ax; xlim = nothing, ylim = nothing, xlabel = "", ylabel = "", yreversed = false)

    # Set axis labels
    ax.xlabel = xlabel
    ax.ylabel = ylabel
    ax.yreversed = yreversed
    
    # Set axis limits
    xlim !== nothing && xlims!(ax, xlim[1], xlim[2])
    ylim !== nothing && ylims!(ax, ylim[1], ylim[2])

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

# =============================================================================
# WINDOW TITLE AND STRING ABBREVIATION UTILITIES
# =============================================================================

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
        if m.captures[3] !== nothing  # Underscore
            push!(parts, "_")
        elseif m.captures[4] !== nothing  # Standalone digits
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
function _generate_window_title( datasets::Vector{<:EegData}; max_total_length::Int = 80, max_name_length::Int = 20)

    isempty(datasets) && return ""
    length(datasets) == 1 && return "$(datasets[1].file):$(datasets[1].condition_name)" 

    # Check if all datasets have the same file name
    first_file = datasets[1].file
    all_same_file = all(dataset.file == first_file for dataset in datasets)
    
    if all_same_file 
        condition_names = [data.condition_name for data in datasets]
        condition_str = _shorten_condition_names(condition_names; 
                                                max_total_length = max_total_length - length(first_file) - 2, 
                                                max_name_length = max_name_length)
        return "$first_file: $condition_str"
    else
        file_condition_pairs = ["$(data.file):$(data.condition_name)" for data in datasets]
        return _shorten_condition_names(file_condition_pairs; 
                                       max_total_length = max_total_length, 
                                       max_name_length = max_name_length)
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
    unique_suffix = join(name_parts[length(common_prefix_parts)+1:end], "")
    
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
        last_part = join(final_names[(end - show_ends + 1):end], separator)
        return first_part * separator * "…" * separator * last_part
    end
    
    return full_string
end

# =============================================================================
# BASELINE HELPER FUNCTIONS
# =============================================================================

"""
    _extract_baseline_values(interval::BaselineInterval)

Extract numeric baseline values from interval.
Returns (start, stop) tuple or (nothing, nothing) if interval is nothing.
"""
function _extract_baseline_values(interval::BaselineInterval)
    if interval === nothing
        return nothing, nothing
    else
        return interval.start, interval.stop
    end
end

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
