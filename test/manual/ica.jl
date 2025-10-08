using eegfun
using GLMakie
using BenchmarkTools

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..",  "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
eegfun.polar_to_cartesian_xy!(layout_file)
dat = eegfun.read_bdf(data_file);
dat = eegfun.create_eeg_dataframe(dat, layout_file);
eegfun.rereference!(dat, :avg)
eegfun.filter_data!(dat, "hp", 1)
eegfun.is_extreme_value!(dat, 100);
ica_result = eegfun.run_ica(dat; sample_selection = eegfun.samples_not(:is_extreme_value_100))

# plot_ica_topoplot
eegfun.plot_ica_topoplot(ica_result)
eegfun.plot_ica_topoplot(ica_result, component_selection = eegfun.components(1:4), plot_labels = false);
eegfun.plot_ica_topoplot(ica_result, component_selection = eegfun.components(1:10), plot_colorbar = true);
eegfun.plot_ica_topoplot(ica_result, component_selection = eegfun.components(1:10), plot_colorbar = true, colorbar_plot_numbers = [1]);
eegfun.plot_ica_topoplot(ica_result, component_selection = eegfun.components(1:10), use_global_scale = true);
eegfun.plot_ica_topoplot(ica_result; component_selection = eegfun.components(1:10), method = :spherical_spline)
eegfun.plot_ica_topoplot(ica_result; component_selection = eegfun.components(1:10), method = :multiquadratic)

eegfun.plot_ica_topoplot(
    ica_result,
    component_selection = eegfun.components([1, 3, 5]);
    use_global_scale = true,
    plot_colorbar = true,
    colorbar_plot_numbers = [2],
)

eegfun.plot_ica_topoplot(
    ica_result,
    component_selection = eegfun.components([1, 3, 5, 7, 9]);
    dims = (2, 3),
    use_global_scale = true,
    plot_colorbar = true,
    colorbar_plot_numbers = [5],
)


# plot_ica_component_activation
eegfun.plot_ica_component_activation(dat, ica_result, method = :multiquadratic)



GLMakie.closeall()