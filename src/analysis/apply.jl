"""
    apply_analysis_settings!(data::EegData, settings::AnalysisSettings)

Apply analysis settings to data in-place.
"""
function apply_analysis_settings!(dat::EegData, settings::AnalysisSettings)

    settings.hp_filter != 0 && highpass_filter!(dat, settings.hp_filter)
    settings.lp_filter != 0 && lowpass_filter!(dat, settings.lp_filter)
    settings.reference != :none && rereference!(dat, settings.reference)
    !isempty(settings.repaired_channels) &&
        repair_channels!(dat, settings.repaired_channels, method = settings.repair_method)
    !isempty(settings.selected_regions) && _add_selected_regions!(dat, settings.selected_regions)

    return nothing
end
apply_analysis_settings!(dat::EegData, settings::Observable{AnalysisSettings}) =
    apply_analysis_settings!(dat, settings[])

"""
    apply_analysis_settings!(data::EegData, settings::AnalysisSettings, ica::InfoIca)

Apply analysis settings to data in-place, including ICA component removal.
"""
function apply_analysis_settings!(dat::EegData, ica::InfoIca, settings::AnalysisSettings)

    settings.hp_filter != 0 && highpass_filter!(dat, settings.hp_filter)
    settings.lp_filter != 0 && lowpass_filter!(dat, settings.lp_filter)
    settings.reference != :none && rereference!(dat, settings.reference)
    !isempty(settings.repaired_channels) &&
        repair_channels!(dat, settings.repaired_channels, method = settings.repair_method)
    !isempty(settings.selected_regions) && _add_selected_regions!(dat, settings.selected_regions)

    # ICA component removal if selected
    if !isempty(settings.removed_ica_components)
        remove_ica_components!(dat, ica, component_selection = components(settings.removed_ica_components))
    end

    return nothing
end
apply_analysis_settings!(dat::EegData, ica::InfoIca, settings::Observable{AnalysisSettings}) =
    apply_analysis_settings!(dat, ica, settings[])

@add_nonmutating apply_analysis_settings!
