using EegFun
using JLD2
dat = EegFun.Dataset_3D_Shape()

# run ica
ica_result = EegFun.run_ica(dat; sample_selection = EegFun.samples_not(:is_extreme_value_100), percentage_of_data = 25)

# simulate applying ICA component 1
state = EegFun.ContinuousDataBrowserState(dat, ica_result, Dict{Symbol,Any}())
EegFun._toggle_ica_component!(state, 1)

# Check the removed components
println("Removed: ", state.ica.removed_components)

# Compute what subtracted_data_obs would compute
col = state.channels.labels[1]
current_raw = EegFun.get_data(state.data.current[], 1:100, col)
original_raw = copy(current_raw)
for comp in state.ica.removed_components
    unmix_vec = state.ica.original.unmixing[comp, :]
    mix_weight = state.ica.original.mixing[findfirst(==(col), state.ica.original.layout.data.label), comp]
    activation = zeros(Float64, 100)
    for (ci, col_sym) in enumerate(state.ica.original.layout.data.label)
        norm_ch = (EegFun.get_data(state.data.original, 1:100, col_sym) .- state.ica.original.mean[ci]) ./ state.ica.original.scale
        activation .+= unmix_vec[ci] .* norm_ch
    end
    original_raw .+= mix_weight .* state.ica.original.scale .* activation
end
println("Activation max diff: ", maximum(abs.(original_raw .- current_raw)))

# Test toggle nothing and then 1
EegFun._toggle_ica_component!(state, nothing)
EegFun._toggle_ica_component!(state, 1)

current_raw = EegFun.get_data(state.data.current[], 1:100, col)
original_raw = copy(current_raw)
for comp in state.ica.removed_components
    unmix_vec = state.ica.original.unmixing[comp, :]
    mix_weight = state.ica.original.mixing[findfirst(==(col), state.ica.original.layout.data.label), comp]
    activation = zeros(Float64, 100)
    for (ci, col_sym) in enumerate(state.ica.original.layout.data.label)
        norm_ch = (EegFun.get_data(state.data.original, 1:100, col_sym) .- state.ica.original.mean[ci]) ./ state.ica.original.scale
        activation .+= unmix_vec[ci] .* norm_ch
    end
    original_raw .+= mix_weight .* state.ica.original.scale .* activation
end
println("Activation max diff after toggle nothing: ", maximum(abs.(original_raw .- current_raw)))
