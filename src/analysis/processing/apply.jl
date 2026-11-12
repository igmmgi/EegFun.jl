"""
    apply_analysis_settings!(dat::EegData, settings::Union{AnalysisSettings, Observable{AnalysisSettings}}) -> EegData
    apply_analysis_settings!(dat::EegData, ica::InfoIca, settings) -> EegData

Apply interactive `AnalysisSettings` (filtering, rereferencing, channel repair, region selection,
and optional ICA component subtraction) to `dat` in place.
"""
function apply_analysis_settings!(dat::EegData, settings::AnalysisSettings)
    settings.hp_filter != 0 && highpass_filter!(dat, settings.hp_filter)
    settings.lp_filter != 0 && lowpass_filter!(dat, settings.lp_filter)
    settings.reference != :none && rereference!(dat, settings.reference)
    !isempty(settings.repaired_channels) && repair_channels!(dat, settings.repaired_channels; method = settings.repair_method)
    !isempty(settings.selected_regions) && _add_selected_regions!(dat, settings.selected_regions)
    return dat
end

apply_analysis_settings!(dat::EegData, settings::Observable{AnalysisSettings}) = apply_analysis_settings!(dat, settings[])

function apply_analysis_settings!(dat::EegData, ica::InfoIca, settings::AnalysisSettings)
    apply_analysis_settings!(dat, settings)
    if !isempty(settings.removed_ica_components)
        subtract_ica_components!(dat, ica; component_selection = components(settings.removed_ica_components))
    end
    return dat
end

apply_analysis_settings!(dat::EegData, ica::InfoIca, settings::Observable{AnalysisSettings}) =
    apply_analysis_settings!(dat, ica, settings[])

"""
    apply_analysis_settings(dat::EegData, args...; kwargs...) -> EegData

Non-mutating version of [`apply_analysis_settings!`](@ref).
"""
apply_analysis_settings(dat::EegData, ica::InfoIca, settings) = apply_analysis_settings!(copy(dat), copy(ica), settings)
apply_analysis_settings(dat::EegData, args...; kwargs...) = apply_analysis_settings!(copy(dat), args...; kwargs...)
