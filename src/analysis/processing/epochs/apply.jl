# Replay and application of recorded AnalysisSettings onto EEG datasets.

######################################
########### User Functions ###########
######################################

"""
    apply_analysis_settings(dat::EegData, settings::Union{AnalysisSettings, Observable{AnalysisSettings}}) -> EegData
    apply_analysis_settings(dat::EegData, ica::InfoIca, settings) -> EegData

Apply interactive `AnalysisSettings` (filtering, rereferencing, channel repair, region selection,
and optional ICA component subtraction) to a copy of `dat`, leaving `dat` unmodified.
"""
apply_analysis_settings(dat::EegData, ica::InfoIca, settings) = apply_analysis_settings!(copy(dat), copy(ica), settings)

apply_analysis_settings(dat::EegData, args...; kwargs...) = apply_analysis_settings!(copy(dat), args...; kwargs...)

"""
    apply_analysis_settings!(dat::EegData, settings::Union{AnalysisSettings, Observable{AnalysisSettings}}) -> EegData
    apply_analysis_settings!(dat::EegData, ica::InfoIca, settings) -> EegData

In-place mutating version of [`apply_analysis_settings`](@ref).
"""
function apply_analysis_settings!(dat::EegData, settings::AnalysisSettings)
    _apply_repair_and_reference!(dat, settings)
    _apply_filters_and_regions!(dat, settings)
    return dat
end

apply_analysis_settings!(dat::EegData, settings::Observable{AnalysisSettings}) = apply_analysis_settings!(dat, settings[])

function apply_analysis_settings!(dat::EegData, ica::InfoIca, settings::AnalysisSettings)
    _apply_repair_and_reference!(dat, settings)
    if !isempty(settings.removed_ica_components)
        subtract_ica_components!(dat, copy(ica); component_selection = components(settings.removed_ica_components))
    end
    _apply_filters_and_regions!(dat, settings)
    return dat
end

apply_analysis_settings!(dat::EegData, ica::InfoIca, settings::Observable{AnalysisSettings}) =
    apply_analysis_settings!(dat, ica, settings[])

######################################
######### Internal Functions #########
######################################

function _apply_repair_and_reference!(dat::EegData, settings::AnalysisSettings)
    !isempty(settings.repaired_channels) && repair_channels!(dat, settings.repaired_channels; method = settings.repair_method)
    if settings.reference != :none
        ref_target =
            (settings.reference == :custom && !isempty(settings.reference_channels)) ? settings.reference_channels : settings.reference
        rereference!(dat, ref_target)
    end
    return dat
end

function _apply_filters_and_regions!(dat::EegData, settings::AnalysisSettings)
    settings.hp_filter != 0 && highpass_filter!(dat, settings.hp_filter)
    settings.lp_filter != 0 && lowpass_filter!(dat, settings.lp_filter)
    !isempty(settings.selected_regions) && _add_selected_regions!(dat, settings.selected_regions)
    return dat
end
