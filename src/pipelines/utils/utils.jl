# Helper functions used during preprocess to make logging easier/prettier
function _applied_filters(filter_cfg::FilterConfig; filter_sections::Vector{Symbol} = [:highpass, :lowpass])
    applied_filters = []
    for name in filter_sections
        section = getproperty(filter_cfg, name)
        if section.apply
            push!(
                applied_filters,
                "$(section.freq) Hz $(section.type), $(section.method), $(section.func), order $(section.order)",
            )
        end
    end
    return join(applied_filters, "; ")
end

function _eog_config_string(eog_cfg::EogConfig)

    format_channels(channels) = begin
        left = join(channels[1], ", ")
        right = join(channels[2], ", ")
        ref = channels[3][1]
        return "$left - $right → $ref"
    end

    vEOG = format_channels(eog_cfg.vEOG_channels)
    hEOG = format_channels(eog_cfg.hEOG_channels)
return "vEOG: $vEOG ($(eog_cfg.vEOG_criterion) μV); hEOG: $hEOG ($(eog_cfg.hEOG_criterion) μV)"
end

_flag_symbol(base::AbstractString, criterion) = Symbol("$(base)_$(criterion)")

function _center_title(title::String, width::Int)
    total_dashes = width - length(title) - 2  # 2 spaces around title
    left_dashes = div(total_dashes, 2)
    right_dashes = total_dashes - left_dashes
    return "-" ^ left_dashes * " $title " * "-" ^ right_dashes
end

function section(title::String; width::Int = 80)
    dash_line = "-" ^ width
    middle_line = _center_title(title, width)
    return "\n$dash_line\n$middle_line\n$dash_line"
end

subsection(title::String; width::Int = 80) = "\n" * _center_title(title, width)
subsubsection(title::String) = "\n# " * title
