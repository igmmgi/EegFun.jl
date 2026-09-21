"""
    plot_lmm_heatmap(result::LmmStatsResult; coef_idx::Int=2, threshold_p=nothing, kwargs...)

Plot a 2D Time x Channels heatmap of LMM statistical values.
"""
function plot_lmm_heatmap(result::LmmStatsResult; coef_idx::Int=2, threshold_p=nothing, kwargs...)
    plot_kwargs = Dict{Symbol, Any}(kwargs)
    
    # Initialize figure
    fontsize_kw = haskey(plot_kwargs, :theme_fontsize) && !isnothing(plot_kwargs[:theme_fontsize]) ? (; fontsize=plot_kwargs[:theme_fontsize]) : (;)
    fig = Figure(; fontsize_kw...)
    ax = Axis(fig[1, 1], xlabel="Time (s)", ylabel="Channels")
    
    t_vals = result.t_values[:, :, coef_idx]
    
    if !isnothing(threshold_p)
        p_vals = result.p_values[:, :, coef_idx]
        mask = p_vals .> threshold_p
        t_vals_plot = copy(t_vals)
        t_vals_plot[mask] .= NaN
    else
        t_vals_plot = t_vals
    end
    
    # Compute limits for a symmetric diverging colormap
    max_abs_t = maximum(abs, filter(!isnan, t_vals_plot); init=0.1)
    clims = get(plot_kwargs, :colorrange, (-max_abs_t, max_abs_t))
    cmap = get(plot_kwargs, :colormap, :RdBu)
    
    hm = heatmap!(ax, result.time_points, 1:length(result.channels), t_vals_plot'; 
                  colormap=cmap, colorrange=clims)
                  
    ax.yticks = (1:length(result.channels), string.(result.channels))
    
    coef_name = result.coefficients[coef_idx]
    ax.title = get(plot_kwargs, :plot_title, "LMM Heatmap: \$(coef_name)")
    
    # Add colorbar
    Colorbar(fig[1, 2], hm, label="t-value")
    
    return fig
end

"""
    plot_lmm_topomap(result::LmmStatsResult, time_point::Real; coef_idx::Int=2, kwargs...)

Plot a topographic map of LMM statistical values at a specific time point.
Delegates to the native `EegFun.plot_topography` using a temporary EEG object.
"""
function plot_lmm_topomap(result::LmmStatsResult, time_point::Real; coef_idx::Int=2, kwargs...)
    # Find closest time point index
    t_idx = argmin(abs.(result.time_points .- time_point))
    
    # Extract t-values for all channels
    t_vals = result.t_values[:, t_idx, coef_idx]
    
    # Construct a temporary DataFrame
    df = DataFrame(time = [result.time_points[t_idx]])
    for (i, ch) in enumerate(result.channels)
        df[!, ch] = [t_vals[i]]
    end
    
    # Construct a temporary ContinuousData object to reuse native plotting
    temp_eeg = ContinuousData("LMM_Topomap", df, result.epochs.layout, result.epochs.sample_rate, AnalysisInfo())
    
    plot_kwargs = Dict{Symbol, Any}(kwargs)
    if !haskey(plot_kwargs, :colormap)
        plot_kwargs[:colormap] = :RdBu
    end
    
    coef_name = result.coefficients[coef_idx]
    if !haskey(plot_kwargs, :plot_title)
        plot_kwargs[:plot_title] = "LMM Topomap: \$(coef_name) at \$(round(result.time_points[t_idx], digits=3))s"
    end
    
    if !haskey(plot_kwargs, :ylim)
        max_abs_t = maximum(abs, t_vals; init=0.1)
        plot_kwargs[:ylim] = (-max_abs_t, max_abs_t)
    end
    
    return plot_topography(temp_eeg; plot_kwargs...)
end
