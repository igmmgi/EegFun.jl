using BioSemiBDF
using CSV
using DataFrames
using DSP
using GLMakie
using LibGEOS
using LinearAlgebra
using OrderedCollections
using Random
using Statistics
using StatsBase: kurtosis
using FFTW
using Printf


"""
    IcaPrms

Structure containing ICA algorithm parameters.

# Fields
- `l_rate::Float64`: Initial learning rate (default: 0.001)
- `max_iter::Int`: Maximum number of iterations (default: 512)
- `w_change::Float64`: Change threshold for stopping (default: 1e-6)
- `anneal_deg::Float64`: Angle for learning rate reduction (default: 60.0)
- `anneal_step::Float64`: Learning rate reduction factor (default: 0.9)
- `blowup::Float64`: Maximum weight change allowed (default: 1e15)
- `blowup_fac::Float64`: Learning rate reduction on blowup (default: 0.8)
- `max_weight::Float64`: Maximum weight magnitude (default: 1e8)
- `restart_factor::Float64`: Learning rate factor on restart (default: 0.9)
- `degconst::Float64`: Degrees to radians conversion (default: 180.0 / π)
- `default_stop::Float64`: Default weight stop criterion (default: 1e-6)
"""
mutable struct IcaPrms
    l_rate::Float64
    max_iter::Int
    w_change::Float64
    anneal_deg::Float64
    anneal_step::Float64
    blowup::Float64
    blowup_fac::Float64
    max_weight::Float64
    restart_factor::Float64
    degconst::Float64
    default_stop::Float64
end

# Outer constructor with default values
function IcaPrms(;
    l_rate = 0.001,
    max_iter = 512,
    w_change = 1e-6,
    anneal_deg = 60.0,
    anneal_step = 0.9,
    blowup = 1e15,
    blowup_fac = 0.8,
    max_weight = 1e8,
    restart_factor = 0.9,
    degconst = 180.0 / π,
    default_stop = 1e-6,
)
    IcaPrms(
        l_rate,
        max_iter,
        w_change,
        anneal_deg,
        anneal_step,
        blowup,
        blowup_fac,
        max_weight,
        restart_factor,
        degconst,
        default_stop,
    )
end

"""
    InfoIca

Structure containing ICA decomposition results.

# Fields
- `unmixing::Matrix{Float64}`: Transforms data to ICs
- `mixing::Matrix{Float64}`: Transforms ICs back to data
- `scale::Float64`: Data scaling factor
- `variance::Vector{Float64}`: Variance explained by each component
- `kurtosis::Vector{Float64}`: Kurtosis of each component
- `spectra::Matrix{Float64}`: Power spectra of components (freqs × components)
- `frequencies::Vector{Float64}`: Frequency bins for spectra
- `ica_label::Vector{String}`: Component labels
- `data_label::Vector{String}`: Original data channel labels
"""
struct InfoIca
    unmixing::Matrix{Float64}
    mixing::Matrix{Float64}
    scale::Float64
    variance::Vector{Float64}
    ica_label::Vector{String}
    data_label::Vector{String}
end



"""
    infomax_ica(data::Matrix{Float64}; params::IcaPrms = IcaPrms(), extended::Bool = false, ext_params::Union{Nothing, ExtendedIcaPrms} = nothing, n_components::Union{Nothing,Int} = nothing)

Perform Infomax ICA decomposition on EEG data using the algorithm from EEGLAB's runica.

# Arguments
- `data::Matrix{Float64}`: Data matrix (channels × samples)
- `params::IcaPrms`: ICA parameters
- `extended::Bool`: Whether to use extended Infomax
- `ext_params::ExtendedIcaPrms`: Extended ICA parameters
- `n_components::Union{Nothing,Int}`: Number of components for PCA reduction (default: no reduction)

# Returns
- `InfoICA`: Structure containing:
  - `topo`: Topography matrix (mixing matrix)
  - `unmixing`: Unmixing matrix (weights)
  - `label`: Component labels
"""

function to_data_frame(dat::EpochData)
    return vcat(dat.data...)
end

function to_data_frame(dat::Vector{EpochData})
    return vcat([vcat(dat[idx].data[:]...) for idx in eachindex(dat)]...)
end

function create_ica_data_matrix(dat::DataFrame, channels; samples_to_include = nothing)
    if isnothing(samples_to_include)
        samples_to_include = dat.sample
    end
    dat = dat[samples_to_include, :]
    # need to make sure we have unique samples as with longer epochs there is potential for overlap
    dat = unique(dat, :sample)
    # select only the channels we want
    dat = dat[!, intersect(names(dat), channels)]
    return permutedims(Matrix(dat))
end




# function infomax_ica(
#     dat::ContinuousData,
#     data_labels;
#     n_components::Union{Nothing,Int} = nothing,
#     params::IcaPrms = IcaPrms(),
# )
#     # select actual eeg data columns
#     dat = permutedims(Matrix(dat.data[!, intersect(names(dat.data), data_labels)]))
#     return infomax_ica(dat, data_labels, params = params, n_components = n_components)
# end


# Pre-allocate all arrays in a single struct for better cache locality
mutable struct WorkArrays
    weights::Matrix{Float64}
    BI::Matrix{Float64}
    u::Matrix{Float64}
    y::Matrix{Float64}
    data_block::Matrix{Float64}
    oldweights::Matrix{Float64}
    startweights::Matrix{Float64}
    weights_temp::Matrix{Float64}
    y_temp::Matrix{Float64}
    bi_weights::Matrix{Float64}
    wu_term::Matrix{Float64}
    delta::Matrix{Float64}
    olddelta::Matrix{Float64}
end

function create_work_arrays(n_components::Int, block_size::Int)
    weights = Matrix{Float64}(I, n_components, n_components)  # Initialize as identity matrix
    return WorkArrays(
        weights,  # weights - start as identity matrix
        block_size * Matrix{Float64}(I, n_components, n_components),  # BI
        zeros(n_components, block_size),  # u
        zeros(n_components, block_size),  # y
        zeros(n_components, block_size),  # data_block
        copy(weights),  # oldweights - copy of identity matrix
        copy(weights),  # startweights - copy of identity matrix
        zeros(n_components, n_components),  # weights_temp
        zeros(n_components, block_size),  # y_temp
        zeros(n_components, n_components),  # bi_weights
        zeros(n_components, n_components),  # wu_term
        zeros(1, n_components^2),  # delta
        zeros(1, n_components^2)   # olddelta
    )
end


function infomax_ica(
    dat_ica::Matrix{Float64},
    data_labels;
    n_components::Union{Nothing,Int} = nothing,
    params::IcaPrms = IcaPrms(),
    sample_rate::Int = 256
)

    # demean and scale data
    dat_ica .-= mean(dat_ica, dims = 2)
    scale = sqrt(norm((dat_ica * dat_ica') / size(dat_ica, 2)))
    dat_ica ./= scale

    # PCA reduction
    if !isnothing(n_components)
        F = svd(dat_ica)
        # Ensure we're using the correct dimensions for PCA components
        pca_components = F.U[1:size(dat_ica,1), 1:n_components]  # Explicitly specify both dimensions
        dat_ica = pca_components' * dat_ica
    else
        n_components = size(dat_ica, 1)
        pca_components = nothing
    end

    # sphering 
    sphere = inv(sqrt(cov(dat_ica, dims = 2)))
    dat_ica = sphere * dat_ica

    # initialize
    n_channels = size(dat_ica, 1)
    n_samples = size(dat_ica, 2)
    block = min(Int(floor(sqrt(n_samples / 3.0))), 512)
    work = create_work_arrays(n_channels, block)
    lastt = block * div(n_samples, block)
    step = 0
    wts_blowup = false
    change = 0.0
    oldchange = 0.0
    angledelta = 0.0

    # pre-allocate permutation vector
    permute_indices = Vector{Int}(undef, n_samples)

    @inbounds while step < params.max_iter
        randperm!(permute_indices)

        for t = 1:block:lastt
            block_end = min(t + block - 1, n_samples)
            block_size = block_end - t + 1

            # extract data block
            copyto!(view(work.data_block, :, 1:block_size), view(dat_ica, :, view(permute_indices, t:block_end)))

            # forward pass
            mul!(work.u, work.weights, work.data_block)
            @. work.y = 1 / (1 + exp(-work.u))
            @. work.y_temp = 1 - 2 * work.y

            # update weights 
            mul!(work.wu_term, work.y_temp, transpose(work.u))
            work.bi_weights .= work.BI .+ work.wu_term
            mul!(work.weights_temp, work.bi_weights, work.weights)
            @. work.weights += params.l_rate * work.weights_temp

            # boom?
            if maximum(abs, work.weights) > params.max_weight
                wts_blowup = true
                change = NaN
                break
            end
        end

        if !wts_blowup
            work.oldweights .-= work.weights
            step += 1
            work.delta .= reshape(work.oldweights, 1, :)
            change = dot(work.delta, work.delta)
        end

        if wts_blowup || isnan(change) || isinf(change)
            step = 0
            change = NaN
            wts_blowup = false
            params.l_rate *= params.restart_factor
            work.weights .= work.startweights
            work.oldweights .= work.startweights
            continue
        end

        if step > 2
            angledelta = acos(clamp(dot(work.delta, work.olddelta) / sqrt(change * oldchange), -1, 1))
            if params.degconst * angledelta > params.anneal_deg
                params.l_rate *= params.anneal_step
                work.olddelta .= work.delta
                oldchange = change
            end
        elseif step == 1
            work.olddelta .= work.delta
            oldchange = change
        end

        work.oldweights .= work.weights

        if step > 2 && change < params.w_change
            break
        elseif change > params.blowup
            params.l_rate *= params.blowup_fac
        end

        @info "Step $step, change = $change, lrate = $(params.l_rate), angle = $((params.degconst) * angledelta)"

    end

    # Final calculations
    if !isnothing(pca_components)
        work.weights = work.weights * sphere * pca_components'
    else
        work.weights = work.weights * sphere
    end

    mixing = pinv(work.weights)

    # Calculate total variance explained
    meanvar = vec(sum(abs2, mixing, dims = 1) .* sum(abs2, dat_ica, dims = 2)' ./ (n_components * n_samples - 1))
    meanvar_normalized = meanvar ./ sum(meanvar)
    order = sortperm(meanvar_normalized, rev = true)

    return InfoIca(
        work.weights[order, :],
        mixing[:, order],
        scale,
        meanvar_normalized[order],
        ["IC$i" for i = 1:size(work.weights, 1)],
        data_labels,
    )
end





function plot_ica_component_activation(
    dat::ContinuousData,
    ica_result::InfoIca,
    n_visible_components::Int = 10,  # Number of components visible at once
)
    # convert ContinuousData to appropriate matrix
    dat_matrix = permutedims(Matrix(dat.data[!, ica_result.data_label]))

    # Scale dat matrix the same way as in ICA
    dat_matrix .-= mean(dat_matrix, dims = 2)
    dat_matrix ./= ica_result.scale

    # Transform data to component space
    components = ica_result.unmixing * dat_matrix
    total_components = size(components, 1)

    # Create figure
    fig = Figure()

    # Create main layout
    # main_layout = fig[1,1] = GridLayout()
    
    # # Create plot area (90% of space)
    # plot_area = main_layout[1,1] = GridLayout()
    # rowsize!(main_layout, 1, Relative(0.9))
    
    # Create slider area (10% of space)
    # slider_area = main_layout[2,1] = GridLayout(tellheight=false)
    # rowsize!(main_layout, 2, Auto())
    
    # Create plot grids
    # gl = plot_area[1,1] = GridLayout()

    # Set up observables for interactive plotting
    window_size = 2000
    xrange = Observable(1:window_size)  # Initial window
    xlims = @lift((dat.data.time[first($xrange)], dat.data.time[last($xrange)]))
    
    # Initialize y-axis limits based on initial data window
    initial_range = maximum(abs.(extrema(components[1:n_visible_components, 1:window_size])))
    ylims = Observable((-initial_range, initial_range))
    
    # Observable for component range
    comp_start = Observable(1)
    
    # Create all subplots at once
    axs = []  # Store time series axes
    lines_obs = []  # Store line observables
    v_eogs = []  # Store line observables
    h_eogs = []  # Store line observables
    topo_axs = []  # Store topography axes
    
    for i in 1:n_visible_components
        # Create subplot for topography
        ax_topo = Axis(
            fig[i, 1],
            width = 150,
            height = 150,
            title = @sprintf("IC %d (%.1f%%)", i, ica_result.variance[i] * 100)
        )
        push!(topo_axs, ax_topo)
        
        # Initial topography plot
        plot_ica_topoplot(fig, ax_topo, ica_result, i, layout, colorbar_kwargs = Dict(:plot_colorbar => false))
        
        hidexdecorations!(ax_topo)

        # Create subplot for time series
        ax_time = Axis(fig[i, 2], ylabel = "ICA Amplitude")
        push!(axs, ax_time)

        # Set initial limits
        xlims!(ax_time, xlims[])
        ylims!(ax_time, ylims[])

        # Create observable for component data
        line_obs = Observable(components[i, :])
        v_eog = Observable(dat.data.vEOG)
        h_eog = Observable(dat.data.hEOG)
        lines!(ax_time, @lift(dat.data.time[$xrange]), @lift($line_obs[$xrange]))
        lines!(ax_time, @lift(dat.data.time[$xrange]), @lift($v_eog[$xrange]))
        lines!(ax_time, @lift(dat.data.time[$xrange]), @lift($h_eog[$xrange]))
        push!(lines_obs, line_obs)
        push!(v_eogs, v_eog)
        push!(h_eogs, h_eog)

        # Hide x-axis decorations for all but the last plot
        if i != n_visible_components
            hidexdecorations!(ax_time, grid = false)
        end
    end

    # Link all time series axes
    linkaxes!(axs...)

    # Adjust layout
    colsize!(fig.layout, 1, Auto(150))  # Fixed width for topo plots

    # Function to update component data
    function update_components(start_idx)
        for i in 1:n_visible_components
            comp_idx = start_idx + i - 1
            if comp_idx <= total_components
                lines_obs[i][] = components[comp_idx, :]
                
                # Clear and redraw topography
                empty!(topo_axs[i])
                plot_ica_topoplot(fig, topo_axs[i], ica_result, comp_idx, layout, colorbar_kwargs = Dict(:plot_colorbar => false))
                
                # Update title
                topo_axs[i].title = @sprintf("IC %d (%.1f%%)", comp_idx, ica_result.variance[comp_idx] * 100)
            end
        end
    end

    # # Add navigation buttons below topo plots
    topo_nav = GridLayout(fig[end+1,1])
    prev_topo = Button(topo_nav[1,1], label = "◄ Previous")
    next_topo = Button(topo_nav[1,2], label = "Next ►")

    # Connect topo navigation buttons
    on(prev_topo.clicks) do _
        new_start = max(1, comp_start[] - n_visible_components)
        comp_start[] = new_start
        update_components(new_start)
    end

    on(next_topo.clicks) do _
        new_start = min(total_components - n_visible_components + 1, comp_start[] + n_visible_components)
        comp_start[] = new_start
        update_components(new_start)
    end

             # Add keyboard controls
    on(events(fig).keyboardbutton) do event
        if event.action in (Keyboard.press,)
            if event.key == Keyboard.left || event.key == Keyboard.right
                # Handle x-axis scrolling
                current_range = xrange[]
                if event.key == Keyboard.left
                    new_start = max(1, first(current_range) - window_size)
                    xrange[] = new_start:(new_start + window_size - 1)
                else  # right
                    new_start = min(size(components, 2) - window_size + 1, first(current_range) + window_size)
                    xrange[] = new_start:(new_start + window_size - 1)
                end
                
                # Update x-axis limits for all axes
                new_xlims = (dat.data.time[first(xrange[])], dat.data.time[last(xrange[])])
                for ax in axs
                    xlims!(ax, new_xlims)
                end
                
            elseif event.key == Keyboard.up || event.key == Keyboard.down
                # Handle y-axis scaling
                current_range = ylims[][2]  # Just take the positive limit since it's symmetric
                if event.key == Keyboard.up
                    # Zoom in - decrease range by 20%
                    new_range = current_range * 0.8
                else  # down
                    # Zoom out - increase range by 20%
                    new_range = current_range * 1.2
                end
                
                # Keep centered on zero
                new_ylims = (-new_range, new_range)
                ylims[] = new_ylims
                
                # Update y-axis limits for all axes
                for ax in axs
                    ylims!(ax, new_ylims)
                end
                
            elseif event.key == Keyboard.page_up || event.key == Keyboard.page_down
                # Handle component scrolling
                current_start = comp_start[]
                if event.key == Keyboard.page_up
                    new_start = max(1, current_start - n_visible_components)
                else  # page_down
                    new_start = min(total_components - n_visible_components + 1, current_start + n_visible_components)
                end
                
                if new_start != current_start
                    comp_start[] = new_start
                    update_components(new_start)
                end
            end
        end
    end

    return fig
end


# function plot_ica_component_activation(
#     dat::ContinuousData,
#     ica_result::InfoIca,
#     n_visible_components::Int = 10,  # Number of components visible at once
# )
#     # convert ContinuousData to appropriate matrix
#     dat_matrix = permutedims(Matrix(dat.data[!, ica_result.data_label]))
# 
#     # Scale dat matrix the same way as in ICA
#     dat_matrix .-= mean(dat_matrix, dims = 2)
#     dat_matrix ./= ica_result.scale
# 
#     # Transform data to component space
#     components = ica_result.unmixing * dat_matrix
#     total_components = size(components, 1)
# 
#     # Create figure
#     fig = Figure()
# 
#     # Create main layout
#     # main_layout = fig[1,1] = GridLayout()
#     
#     # # Create plot area (90% of space)
#     # plot_area = main_layout[1,1] = GridLayout()
#     # rowsize!(main_layout, 1, Relative(0.9))
#     
#     # Create slider area (10% of space)
#     # slider_area = main_layout[2,1] = GridLayout(tellheight=false)
#     # rowsize!(main_layout, 2, Auto())
#     
#     # Create plot grids
#     # gl = plot_area[1,1] = GridLayout()
# 
#     # Set up observables for interactive plotting
#     window_size = 2000
#     xrange = Observable(1:window_size)  # Initial window
#     xlims = @lift((dat.data.time[first($xrange)], dat.data.time[last($xrange)]))
#     
#     # Initialize y-axis limits based on initial data window
#     initial_range = maximum(abs.(extrema(components[1:n_visible_components, 1:window_size])))
#     ylims = Observable((-initial_range, initial_range))
#     
#     # Observable for component range
#     comp_start = Observable(1)
#     
#     # Create all subplots at once
#     axs = []  # Store time series axes
#     lines_obs = []  # Store line observables
#     v_eogs = []  # Store line observables
#     h_eogs = []  # Store line observables
#     topo_axs = []  # Store topography axes
#     
#     for i in 1:n_visible_components
#         # Create subplot for topography
#         ax_topo = Axis(
#             fig[i, 1],
#             width = 150,
#             height = 150,
#             title = @sprintf("IC %d (%.1f%%)", i, ica_result.variance[i] * 100)
#         )
#         push!(topo_axs, ax_topo)
#         
#         # Initial topography plot
#         plot_ica_topoplot(fig, ax_topo, ica_result, i, layout, colorbar_kwargs = Dict(:plot_colorbar => false))
#         
#         hidexdecorations!(ax_topo)
# 
#         # Create subplot for time series
#         ax_time = Axis(fig[i, 2], ylabel = "ICA Amplitude")
#         push!(axs, ax_time)
# 
#         # Set initial limits
#         xlims!(ax_time, xlims[])
#         ylims!(ax_time, ylims[])
# 
#         # Create observable for component data
#         line_obs = Observable(components[i, :])
#         v_eog = Observable(dat.data.vEOG)
#         h_eog = Observable(dat.data.hEOG)
#         lines!(ax_time, @lift(dat.data.time[$xrange]), @lift($line_obs[$xrange]))
#         lines!(ax_time, @lift(dat.data.time[$xrange]), @lift($v_eog[$xrange]))
#         lines!(ax_time, @lift(dat.data.time[$xrange]), @lift($h_eog[$xrange]))
#         push!(lines_obs, line_obs)
#         push!(v_eogs, v_eog)
#         push!(h_eogs, h_eog)
# 
#         # Hide x-axis decorations for all but the last plot
#         if i != n_visible_components
#             hidexdecorations!(ax_time, grid = false)
#         end
#     end
# 
#     # Link all time series axes
#     linkaxes!(axs...)
# 
#     # Adjust layout
#     colsize!(fig.layout, 1, Auto(150))  # Fixed width for topo plots
# 
#     # Function to update component data
#     function update_components(start_idx)
#         for i in 1:n_visible_components
#             comp_idx = start_idx + i - 1
#             if comp_idx <= total_components
#                 lines_obs[i][] = components[comp_idx, :]
#                 
#                 # Clear and redraw topography
#                 empty!(topo_axs[i])
#                 plot_ica_topoplot(fig, topo_axs[i], ica_result, comp_idx, layout, colorbar_kwargs = Dict(:plot_colorbar => false))
#                 
#                 # Update title
#                 topo_axs[i].title = @sprintf("IC %d (%.1f%%)", comp_idx, ica_result.variance[comp_idx] * 100)
#             end
#         end
#     end
# 
#     # # Add navigation buttons below topo plots
#     topo_nav = GridLayout(fig[end+1,1])
#     prev_topo = Button(topo_nav[1,1], label = "◄ Previous")
#     next_topo = Button(topo_nav[1,2], label = "Next ►")
# 
#     # Connect topo navigation buttons
#     on(prev_topo.clicks) do _
#         new_start = max(1, comp_start[] - n_visible_components)
#         comp_start[] = new_start
#         update_components(new_start)
#     end
# 
#     on(next_topo.clicks) do _
#         new_start = min(total_components - n_visible_components + 1, comp_start[] + n_visible_components)
#         comp_start[] = new_start
#         update_components(new_start)
#     end
# 
#              # Add keyboard controls
#     on(events(fig).keyboardbutton) do event
#         if event.action in (Keyboard.press,)
#             if event.key == Keyboard.left || event.key == Keyboard.right
#                 # Handle x-axis scrolling
#                 current_range = xrange[]
#                 if event.key == Keyboard.left
#                     new_start = max(1, first(current_range) - window_size)
#                     xrange[] = new_start:(new_start + window_size - 1)
#                 else  # right
#                     new_start = min(size(components, 2) - window_size + 1, first(current_range) + window_size)
#                     xrange[] = new_start:(new_start + window_size - 1)
#                 end
#                 
#                 # Update x-axis limits for all axes
#                 new_xlims = (dat.data.time[first(xrange[])], dat.data.time[last(xrange[])])
#                 for ax in axs
#                     xlims!(ax, new_xlims)
#                 end
#                 
#             elseif event.key == Keyboard.up || event.key == Keyboard.down
#                 # Handle y-axis scaling
#                 current_range = ylims[][2]  # Just take the positive limit since it's symmetric
#                 if event.key == Keyboard.up
#                     # Zoom in - decrease range by 20%
#                     new_range = current_range * 0.8
#                 else  # down
#                     # Zoom out - increase range by 20%
#                     new_range = current_range * 1.2
#                 end
#                 
#                 # Keep centered on zero
#                 new_ylims = (-new_range, new_range)
#                 ylims[] = new_ylims
#                 
#                 # Update y-axis limits for all axes
#                 for ax in axs
#                     ylims!(ax, new_ylims)
#                 end
#                 
#             elseif event.key == Keyboard.page_up || event.key == Keyboard.page_down
#                 # Handle component scrolling
#                 current_start = comp_start[]
#                 if event.key == Keyboard.page_up
#                     new_start = max(1, current_start - n_visible_components)
#                 else  # page_down
#                     new_start = min(total_components - n_visible_components + 1, current_start + n_visible_components)
#                 end
#                 
#                 if new_start != current_start
#                     comp_start[] = new_start
#                     update_components(new_start)
#                 end
#             end
#         end
#     end
# 
#     return fig
# end



