using GLMakie


#########################################
# 2D head shape
function head_shape_2d(f, ax, layout; linewidth=2, plot_points=true, plot_labels=true, fontsize=20, markersize=12, label_x_offset=0, label_y_offset=0)

  if (:x2 ∉ names(layout) || :y2 ∉ names(layout))
    polar_to_cartesian_xy!(layout)
  end

  # head shape
  radius = 88 # mm
  arc!(ax, Point2f(0), radius * 2, -π, π, color=:black, linewidth=linewidth) # head
  arc!(Point2f(radius * 2, 0), radius * 2 / 7, -π / 2, π / 2, color=:black, linewidth=linewidth) # ear right
  arc!(Point2f(-radius * 2, 0), -radius * 2 / 7, π / 2, -π / 2, color=:black, linewidth=linewidth) # ear left
  lines!(ax, Point2f[(-0.05, 0.5), (0.0, 0.6), (0.05, 0.5)] .* radius * 4, color=:black, linewidth=linewidth) # nose

  # points
  if plot_points
    scatter!(ax, layout[!, :x2], layout[!, :y2], marker=:circle, markersize=markersize, color=:black)
  end

  if plot_labels
    for label in eachrow(layout)
      text!(ax, fontsize=fontsize, position=(label.x2 + label_x_offset, label.y2 + label_y_offset), label.label)
    end
  end

  # hide some plot stuff
  hidexdecorations!(ax; label=true, ticklabels=true, ticks=true, grid=true, minorgrid=true, minorticks=true)
  hideydecorations!(ax; label=true, ticklabels=true, ticks=true, grid=true, minorgrid=true, minorticks=true)
  hidespines!(ax, :t, :r, :l, :b)

  return f

end

function head_shape_2d(layout; kwargs...)
  f = Figure()
  ax = GLMakie.Axis(f[1, 1])
  head_shape_2d(f, ax, layout; kwargs...)
end


#########################################
# 3D head shape
function head_shape_3d(f, ax, layout; linewidth=2, plot_points=true, plot_labels=true, fontsize=20, markersize=12, label_x_offset=0, label_y_offset=0, label_z_offset=0)

  if (:x3 ∉ names(layout) || :y3 ∉ names(layout) || :z3 ∉ names(layout))
    polar_to_cartesian_xyz!(layout)
  end

  # points
  if plot_points
    scatter!(ax, layout[!, :x3], layout[!, :y3], layout[!, :z3], marker=:circle, markersize=markersize, color=:black)
  end

  if plot_labels
    for label in eachrow(layout)
      text!(ax, fontsize=fontsize, position=(label.x3 + label_x_offset, label.y3 + label_y_offset, label.z3 + label_z_offset), label.label)
    end
  end

  # hide some plot stuff
  hidedecorations!(ax)
  hidespines!(ax)

  return f

end

function head_shape_3d(layout; kwargs...)
  f = Figure()
  ax = GLMakie.Axis3(f[1, 1])
  head_shape_3d(f, ax, layout; kwargs...)
end


##########################################
# 2D topographic plot
function plot_topoplot(dat; ylim=nothing, grid_scale=300, plot_points=true, plot_labels=true, label_x_offset=0, label_y_offset=0)

  radius = 88 # mm

  points = Matrix(dat.layout[!, [:x2, :y2]])'
  data = data_interpolation_topo(Vector(dat.data[1000, dat.layout.label]), points)

  if isnothing(ylim)
    ylim = minimum(data[.!isnan.(data)]), maximum(data[.!isnan.(data)])
  end

  f = Figure()
  ax = GLMakie.Axis(f[1, 1])
  co = contourf!(range(-radius, radius, length=grid_scale), range(-radius, radius, length=grid_scale), data, levels=100, colormap=:jet)
  Colorbar(f[1, 2], co)

  # head shape
  head_shape_2d(f, ax, dat.layout, plot_points=plot_points, plot_labels=plot_labels, label_x_offset=label_x_offset, label_y_offset=label_y_offset)

  return f

end


struct ToggleButton
  label::String
  fun::Function
end

##################################################################
# Data Browser: Continuous Data
function plot_databrowser(dat::ContinuousData, channel_labels::Vector{<:AbstractString})

  function butterfly_plot(active)
    clear_axes()
    if active
      ycentre!()
      draw(plot_labels=false)
    elseif !active
      yoffset!()
      draw(plot_labels=true)
    end
  end

  function apply_lp_filter(active)
    clear_axes()
    if active
      data_filtered = filter_data(data, channel_labels, "lp", slider_lp_filter.value.val, 6, dat.sample_rate)
    elseif !active
      data_filtered = nothing
    end
    draw(plot_labels=true)
  end

  function plot_trigger_lines(active)
    trigger_lines.visible = active
    trigger_text.visible = active
    text_pos = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in trigger_data_time.time]
    trigger_text.position = text_pos
  end

  function plot_vEOG_lines(active) 
    vEOG_lines.visible = active
    vEOG_text.visible = active
    text_pos = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in vEOG_data_time.time]
    vEOG_text.position = text_pos
  end

  function plot_hEOG_lines(active) 
    hEOG_lines.visible = active
    hEOG_text.visible = active
    text_pos = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in hEOG_data_time.time]
    hEOG_text.position = text_pos
  end

  function plot_extreme_lines(active)
    extreme_spans.visible = active
  end

  function toggle_button_group(fig, labels)
    toggles = []
    push!(toggles, ToggleButton("Butterfly Plot", butterfly_plot))
    "triggers" in labels && push!(toggles, ToggleButton("Trigger", plot_trigger_lines))
    "is_vEOG" in labels && push!(toggles, ToggleButton("vEOG", plot_vEOG_lines))
    "is_hEOG" in labels && push!(toggles, ToggleButton("hEOG", plot_hEOG_lines))
    "is_extreme" in labels && push!(toggles, ToggleButton("extreme", plot_extreme_lines))
    push!(toggles, ToggleButton("LP-Filter On/Off", apply_lp_filter))

    toggle_buttons = [Toggle(fig, active=false) for toggle in toggles]
    toggle_labels = [Label(fig, toggle.label, fontsize=30, halign=:left) for toggle in toggles]
    toggle_functions = [toggle.fun for toggle in toggles]

    return hcat(toggle_buttons, toggle_labels, toggle_functions)

  end

  function step_back()
    xrange.val[1] - 200 < 1 && return
    xrange[] = xrange.val .- 200
    xlims!(ax, data.time[xrange.val[1]], data.time[xrange.val[end]])
  end

  function step_forward()
    xrange.val[1] + 200 > nrow(data) && return
    xrange[] = xrange.val .+ 200
    xlims!(ax, data.time[xrange.val[1]], data.time[xrange.val[end]])
  end

  function chans_less()
    (yrange.val[1] + 100 >= 0 || yrange.val[end] - 100 <= 0) && return
    yrange.val = yrange[][1]+100:yrange[][end]-100
    ylims!(ax, yrange.val[1], yrange.val[end])
  end

  function chans_more()
    yrange.val = yrange[][1]-100:yrange[][end]+100
    ylims!(ax, yrange.val[1], yrange.val[end])
  end

  function yoffset!()
    for (idx, col) in enumerate(channel_labels)
      data[!, col] .+= offset.val[idx]
    end
  end

  function ycentre!()
    for (idx, col) in enumerate(channel_labels)
      data[!, col] .-= offset.val[idx]
    end
  end

  function clear_axes()
    [delete!(ax, value) for (key, value) in channel_data_original]
    [delete!(ax, value) for (key, value) in channel_data_filtered]
    [delete!(ax, value) for (key, value) in channel_data_labels]
  end

  # Makie Figure
  fig = Figure()
  ax = GLMakie.Axis(fig[1, 1])

  # controls
  # interactions(ax)
  # deregister_interaction!(ax, :rectanglezoom)
  # deregister_interaction!(ax, :dragpan)
  # deregister_interaction!(ax, :scrollzoom)
  # deregister_interaction!(ax, :limitreset)

  # data to plot
  data = copy(dat.data)
  data_filtered = nothing

  channel_data_original = Dict()
  channel_data_filtered = Dict()
  channel_data_labels = Dict()

  # default xrange/yrange
  xlimit = 10000
  xrange = GLMakie.Observable(1:xlimit)
  yrange = GLMakie.Observable(-1500:1500) 
  nchannels = length(channel_labels)
  channel_labels_original = channel_labels
  channels_to_plot = GLMakie.Observable(channel_labels)

  if nchannels > 1
    offset = GLMakie.Observable(LinRange(yrange.val[end] * 0.9, yrange.val[1] * 0.9, nchannels))
  else # just centre
    offset = GLMakie.Observable(0.0)
  end

  @lift xlims!(ax, data.time[$xrange[1]], data.time[$xrange[end]])
  ylims!(ax, yrange.val[1], yrange.val[end])
  ax.xlabel = "Time (S)"
  ax.ylabel = "Amplitude (mV)"

  # toggle buttons for showing events (triggers, vEOG/hEOG, extreme values ...)
  toggles = toggle_button_group(fig, names(data))
  for t in 1:length(toggles[:,1])
    on(toggles[t,1].active) do _
      toggles[t,3](toggles[t,1].active.val)
    end
  end

  # menu for electrode/channel selection
  menu = hcat(Menu(fig, options=vcat(["All", "Left", "Right", "Central"], dat.layout.label), default="All", direction=:down), Label(fig, "Labels", fontsize=30, halign=:left))
  on(menu[1].selection) do s
    channel_labels = [s]
    if s == "All"
      channel_labels = channel_labels_original
    elseif s == "Left"
      channel_labels = channel_labels_original[findall(occursin.(r"\d*[13579]$", channel_labels_original))]
    elseif s == "Right"
      channel_labels = channel_labels_original[findall(occursin.(r"\d*[24680]$", channel_labels_original))]
    elseif s == "Central"
      channel_labels = channel_labels_original[findall(occursin.(r"z$", channel_labels_original))]
    end
    nchannels = length(channel_labels)
    clear_axes()

    data = copy(dat.data)
    if !isnothing(data_filtered)
      apply_lp_filter(true)
    end
    if nchannels > 1
      offset = GLMakie.Observable(LinRange(yrange.val[end] * 0.9, yrange.val[1] * 0.9, nchannels))
    else # just centre
      offset = GLMakie.Observable(0.0)
    end
    yoffset!()
    if !isnothing(data_filtered)
      apply_lp_filter(true)
    else
      draw(plot_labels=true)
    end
  end

  slider_extreme = Slider(fig[1, 2], range=0:5:200, startvalue=200, width=100)
  slider_lp_filter = Slider(fig[1, 2], range=10:5:100, startvalue=30, width=100)

  crit_val = lift(slider_extreme.value) do x
     x
  end

  slider_range = Slider(fig[3, 1], range=100:30000, startvalue=xlimit, snap=false)
  slider_x = Slider(fig[2, 1], range=slider_range.value.val:nrow(data), startvalue=slider_range.value, snap=false)

  on(slider_x.value) do x
    xrange[] = max(1, x - slider_range.value.val):x-1
    xlims!(ax, data.time[xrange.val[1]], data.time[xrange.val[end]])
  end

  on(slider_range.value) do x
    xrange.val = max(1, slider_x.value.val - x):slider_x.value.val-1
    xlims!(ax, data.time[xrange.val[1]], data.time[xrange.val[end]])
  end

  # keyboard events
  on(events(fig).keyboardbutton) do event
    if event.action in (Keyboard.press, Keyboard.repeat)
      event.key == Keyboard.left && step_back()
      event.key == Keyboard.right && step_forward()
      event.key == Keyboard.down && chans_less()
      event.key == Keyboard.up && chans_more()
    end
    return Consume(false)
  end

  # position GUI controls
  fig[1, 2] = grid!(vcat(toggles[:,1:2],
                         hcat(slider_lp_filter, Label(fig, @lift("LP-Filter Value: $($(slider_lp_filter.value))"), fontsize=30, halign=:left)),
                         hcat(slider_extreme, Label(fig, @lift("Extreme Values: $($(slider_extreme.value)) μV"), fontsize=30)),
                         menu),
                    tellheight=false)
  colsize!(fig.layout, 2, Relative(1 / 8))


  #################### Triggers/Events ###############################
  trigger_data_time = @views data[findall(x -> x != 0, data[!, :].triggers), [:time, :triggers]]
  trigger_lines = vlines!(trigger_data_time.time, color=:grey, linewidth=1, visible = false)
  text_pos = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in trigger_data_time.time]
  trigger_text = text!(string.(trigger_data_time.triggers), position=text_pos, space=:data, align=(:center, :center), fontsize=30, visible=false)

  ################### vEOG/hEOG ###############################
  if ("is_vEOG" in names(dat.data) && "is_hEOG" in names(dat.data))
    vEOG_data_time = @views data[findall(x -> x != 0, data[!, :].is_vEOG), [:time, :is_vEOG]]
    vEOG_lines = vlines!(vEOG_data_time.time, color=:grey, linewidth=1, visible = false)
    text_pos = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in vEOG_data_time.time]
    vEOG_text = text!(repeat(["v"], nrow(vEOG_data_time)), position=text_pos, space=:data, align=(:center, :center), fontsize=30, visible=false)

    hEOG_data_time = @views data[findall(x -> x != 0, data[!, :].is_hEOG), [:time, :is_hEOG]]
    hEOG_lines = vlines!(hEOG_data_time.time, color=:grey, linewidth=1, visible = false)
    text_pos = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in hEOG_data_time.time]
    hEOG_text = text!(repeat(["h"], nrow(hEOG_data_time)), position=text_pos, space=:data, align=(:center, :center), fontsize=30, visible=false)
  end

  ################### Extreme Values ###############################
  if ("is_extreme" in names(data))
    extreme = @views splitgroups(findall(x -> x != 0, data[!, :].is_extreme))
    extreme_spans = vspan!(ax, data[extreme[1], :time], data[extreme[2], :time], color="LightGrey", alpha=0.5, visible=false)
  end

  function draw(; plot_labels=true)
    alpha_orig = 1
    linewidth_orig = 2
    plot_filtered_data = !isnothing(data_filtered)
    if plot_filtered_data
      alpha_orig = 0.5
      linewidth_orig = 1
    end
    for col in channel_labels
      # original data
      channel_data_original[col] = lines!(ax, data[!, :time], data[!, col], color=@lift(abs.(dat.data[!, col]) .>= $crit_val), colormap=[:darkgrey, :darkgrey, :red], linewidth=linewidth_orig, alpha=alpha_orig)
      if plot_labels
        channel_data_labels[col] = text!(ax, @lift(data[$xrange, :time][1]), @lift(data[$xrange, col][1]), text=col, align=(:left, :center), fontsize=20)
      end
      # also show filtered data
      if plot_filtered_data
        channel_data_filtered[col] = lines!(ax, data_filtered[!, :time], data_filtered[!, col], color=:blue, linewidth=2)
      end
    end
  end

  # plot theme adjustments
  fontsize_theme = Theme(fontsize = 24)
  update_theme!(fontsize_theme)

  hideydecorations!(ax, label=true)
  yoffset!()
  draw(plot_labels=true)
  display(fig)
  # DataInspector(fig)

end

function plot_databrowser(dat::ContinuousData)
  plot_databrowser(dat, dat.layout.label)
end

function plot_databrowser(dat::ContinuousData, channel_labels::Union{<:AbstractString,Vector})
  plot_databrowser(dat, [channel_labels])
end



###########################################################

function plot_databrowser(dat::EpochData, channel_labels::Vector{<:AbstractString})

  function butterfly_plot(active)
    clear_axes()
    if active
      ycentre!()
      draw(plot_labels=false)
    elseif !active
      yoffset!()
      draw(plot_labels=true)
    end
  end

  function apply_lp_filter(active)
    clear_axes()
    if active
      data_filtered = filter_data(data, channel_labels, "lp", slider_lp_filter.value.val, 6, dat.sample_rate)
    elseif !active
      data_filtered = nothing
    end
    draw(plot_labels=true)
  end

  function plot_trigger_lines(active)
    trigger_lines.visible = active
    trigger_text.visible = active
    text_pos = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in trigger_data_time.val.time]
    trigger_text.position = text_pos
  end

  function plot_vEOG_lines(active) 
    println(vEOG_data_time.val)
    println(vEOG_lines)
    vEOG_lines.visible = active
    vEOG_text.visible = active
    vEOG_lines[1].val = 1.0
    #text_pos = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in vEOG_data_time.val.time]
    #vEOG_text.position = text_pos
  end

  function plot_hEOG_lines(active) 
    hEOG_lines.visible = active
    hEOG_text.visible = active
    text_pos = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in hEOG_data_time.val.time]
    hEOG_text.position = text_pos
  end

  function plot_extreme_lines(active)
    if length(extreme_spans) > 0
      extreme_spans[1].visible = active
    end
  end

  function toggle_button_group(fig, labels)
    toggles = []
    push!(toggles, ToggleButton("Butterfly Plot", butterfly_plot))
    "triggers" in labels && push!(toggles, ToggleButton("Trigger", plot_trigger_lines))
    "is_vEOG" in labels && push!(toggles, ToggleButton("vEOG", plot_vEOG_lines))
    "is_hEOG" in labels && push!(toggles, ToggleButton("hEOG", plot_hEOG_lines))
    "is_extreme" in labels && push!(toggles, ToggleButton("extreme", plot_extreme_lines))
    push!(toggles, ToggleButton("LP-Filter On/Off", apply_lp_filter))

    toggle_buttons = [Toggle(fig, active=false) for toggle in toggles]
    toggle_labels = [Label(fig, toggle.label, fontsize=30, halign=:left) for toggle in toggles]
    toggle_functions = [toggle.fun for toggle in toggles]

    return hcat(toggle_buttons, toggle_labels, toggle_functions)

  end

  function step_epoch_forward()
    clear_axes()
    trial[] = min(length(dat.data), trial.val[1] + 1)
    ax.title = "Epoch $(trial.val)/$(length(dat.data))"
    update_extreme_spans!()
    draw(trial.val)
  end

  function step_epoch_backward()
    clear_axes()
    trial[] = max(1, trial.val[1] - 1)
    ax.title = "Epoch $(trial.val)/$(length(dat.data))"
    update_extreme_spans!()
    draw(trial.val)
  end

  function chans_less()
    (yrange.val[1] + 100 >= 0 || yrange.val[end] - 100 <= 0) && return
    yrange.val = yrange[][1]+100:yrange[][end]-100
    ylims!(ax, yrange.val[1], yrange.val[end])
  end

  function chans_more()
    yrange.val = yrange[][1]-100:yrange[][end]+100
    ylims!(ax, yrange.val[1], yrange.val[end])
  end

  function yoffset!()
    for t in eachindex(data)
      for (idx, col) in enumerate(channel_labels)
        data[t][!, col] .+= offset.val[idx]
      end
    end
  end

  function ycentre!()
    for t in eachindex(data)
      for (idx, col) in enumerate(channel_labels)
        data[t][!, col] .-= offset.val[idx]
      end
    end
  end

  function clear_axes()
    [delete!(ax, value) for (key, value) in channel_data_original]
    [delete!(ax, value) for (key, value) in channel_data_filtered]
    [delete!(ax, value) for (key, value) in channel_data_labels]
  end

  # Makie Figure
  fig = Figure()
  ax = GLMakie.Axis(fig[1, 1])

  # controls
  # interactions(ax)
  # deregister_interaction!(ax, :rectanglezoom)
  # deregister_interaction!(ax, :dragpan)
  # deregister_interaction!(ax, :scrollzoom)
  # deregister_interaction!(ax, :limitreset)

  # data to plot
  data = deepcopy(dat.data)
  data_filtered = nothing

  channel_data_original = Dict()
  channel_data_filtered = Dict()
  channel_data_labels = Dict()

  # default xrange/yrange
  xlimit = nrow(dat.data[1]) 
  xrange = GLMakie.Observable(1:xlimit)
  trial = GLMakie.Observable(1) # first trial
  yrange = GLMakie.Observable(-1500:1500) 
  nchannels = length(channel_labels)
  channel_labels_original = channel_labels
  channels_to_plot = GLMakie.Observable(channel_labels)

  if nchannels > 1
    offset = GLMakie.Observable(LinRange(yrange.val[end] * 0.9, yrange.val[1] * 0.9, nchannels))
  else # just centre
    offset = GLMakie.Observable(0.0)
  end

  xlims!(ax, data[1].time[xrange.val[1]], data[1].time[xrange.val[end]])
  ylims!(ax, yrange.val[1], yrange.val[end])
  ax.title = "Epoch $(trial.val)/$(length(dat.data))"
  ax.xlabel = "Time (S)"
  ax.ylabel = "Amplitude (mV)"

  # toggle buttons for showing events (triggers, vEOG/hEOG, extreme values ...)
  toggles = toggle_button_group(fig, names(data[1]))
  for t in 1:length(toggles[:,1])
    on(toggles[t,1].active) do _
      toggles[t,3](toggles[t,1].active.val)
    end
  end

  # menu for electrode/channel selection
  menu = hcat(Menu(fig, options=vcat(["All", "Left", "Right", "Central"], dat.layout.label), default="All", direction=:down), Label(fig, "Labels", fontsize=30, halign=:left))
  on(menu[1].selection) do s
    channel_labels = [s]
    if s == "All"
      channel_labels = channel_labels_original
    elseif s == "Left"
      channel_labels = channel_labels_original[findall(occursin.(r"\d*[13579]$", channel_labels_original))]
    elseif s == "Right"
      channel_labels = channel_labels_original[findall(occursin.(r"\d*[24680]$", channel_labels_original))]
    elseif s == "Central"
      channel_labels = channel_labels_original[findall(occursin.(r"z$", channel_labels_original))]
    end
    nchannels = length(channel_labels)
    clear_axes()

    data = deepcopy(dat.data)
    if !isnothing(data_filtered)
      apply_lp_filter(true)
    end
    if nchannels > 1
      offset = GLMakie.Observable(LinRange(yrange.val[end] * 0.9, yrange.val[1] * 0.9, nchannels))
    else # just centre
      offset = GLMakie.Observable(0.0)
    end
    yoffset!()
    if !isnothing(data_filtered)
      apply_lp_filter(true)
    else
      draw(plot_labels=true)
    end
  end

  slider_extreme = Slider(fig[1, 2], range=0:5:200, startvalue=200, width=100)
  slider_lp_filter = Slider(fig[1, 2], range=10:5:100, startvalue=30, width=100)

  crit_val = lift(slider_extreme.value) do x
     x
  end

  # keyboard events
  on(events(fig).keyboardbutton) do event
    if event.action in (Keyboard.press, Keyboard.repeat)
      event.key == Keyboard.left && step_epoch_backward()
      event.key == Keyboard.right && step_epoch_forward()
      event.key == Keyboard.down && chans_less()
      event.key == Keyboard.up && chans_more()
    end
    return Consume(false)
  end

  # position GUI controls
  fig[1, 2] = grid!(vcat(toggles[:,1:2],
                         hcat(slider_lp_filter, Label(fig, @lift("LP-Filter Value: $($(slider_lp_filter.value))"), fontsize=30, halign=:left)),
                         hcat(slider_extreme, Label(fig, @lift("Extreme Values: $($(slider_extreme.value)) μV"), fontsize=30)),
                         menu),
                    tellheight=false)
  colsize!(fig.layout, 2, Relative(1 / 8))


  #################### Triggers/Events ###############################
  trigger_data_time = @lift data[$trial][findall(x -> x != 0, data[$trial][!, :].triggers), [:time, :triggers]]
  trigger_lines = vlines!(trigger_data_time.val.time, color=:grey, linewidth=1, visible = false)
  text_pos = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in trigger_data_time.val.time]
  trigger_text = text!(string.(trigger_data_time.val.triggers), position=text_pos, space=:data, align=(:center, :center), fontsize=30, visible=false)

  ################### vEOG/hEOG ###############################
  if ("is_vEOG" in names(dat.data[1]) && "is_hEOG" in names(dat.data[1]))
    vEOG_data_time = @lift data[$trial][findall(x -> x != 0, data[$trial][!, :].is_vEOG), [:time, :is_vEOG]]
    vEOG_lines = vlines!(vEOG_data_time.val.time, color=:grey, linewidth=1, visible = false)
    text_pos = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in vEOG_data_time.val.time]
    vEOG_text = text!(repeat(["v"], nrow(vEOG_data_time.val)), position=text_pos, space=:data, align=(:center, :center), fontsize=30, visible=false)

    hEOG_data_time = @lift data[$trial][findall(x -> x != 0, data[$trial][!, :].is_hEOG), [:time, :is_hEOG]]
    hEOG_lines = vlines!(hEOG_data_time.val.time, color=:grey, linewidth=1, visible = false)
    text_pos = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in hEOG_data_time.val.time]
    hEOG_text = text!(repeat(["h"], nrow(hEOG_data_time.val)), position=text_pos, space=:data, align=(:center, :center), fontsize=30, visible=false)
  end

  extreme_spans = []
  function update_extreme_spans!()
    if length(extreme_spans) > 0
      delete!(ax, extreme_spans[1])
      extreme_spans = []
    end
    tmp = findall(x -> x != 0, data[trial.val][!, :].is_extreme)
    if length(tmp) > 0
      extreme = splitgroups(tmp)
      if length(extreme) > 0
        # TODO: hard coded Toggle index!!!
        push!(extreme_spans, vspan!(ax, data[trial.val][extreme[1], :time], data[trial.val][extreme[2], :time], color="LightGrey", alpha=0.5, visible=toggles[5,1].active.val))
        println(toggles[5,2].visible)
      end
    end
  end

  #################### Extreme Values ###############################
  if ("is_extreme" in names(data[1]))
    update_extreme_spans!()
  end

  function draw(trial_number; plot_labels=true)
    alpha_orig = 1
    linewidth_orig = 2
    plot_filtered_data = !isnothing(data_filtered)
    if plot_filtered_data
      alpha_orig = 0.5
      linewidth_orig = 1
    end
    for col in channel_labels
      # original data
      channel_data_original[col] = lines!(ax, data[trial_number][!, :time], data[trial_number][!, col], color=@lift(abs.(dat.data[trial_number][!, col]) .>= $crit_val), colormap=[:darkgrey, :darkgrey, :red], linewidth=linewidth_orig, alpha=alpha_orig)
      #if plot_labels
      #  channel_data_labels[col] = text!(ax, @lift(data[trial][$xrange, :time][1]), @lift(data[trial][$xrange, col][$trial]), text=col, align=(:left, :center), fontsize=20)
      #end
      # # also show filtered data
      # if plot_filtered_data
      #   channel_data_filtered[col] = lines!(ax, data_filtered[trial.val][!, :time], data_filtered[trial.val][!, col], color=:blue, linewidth=2)
      # end
    end
  end

  # plot theme adjustments
  fontsize_theme = Theme(fontsize = 24)
  update_theme!(fontsize_theme)

  hideydecorations!(ax, label=true)
  yoffset!()
  draw(1, plot_labels=true)
  display(fig)
  # DataInspector(fig)

end




# # #################################################################
# # Data Browser: Epoched Data
# function plot_databrowser(dat::EpochData, channel_labels::Vector{<:AbstractString})
# 
#   fig = Figure()
#   ax = GLMakie.Axis(fig[1, 1])  # plot layout
# 
#   xrange = GLMakie.Observable(1:nrow(dat.data[1])) # default xrange
#   yrange = GLMakie.Observable(-1500:1500) # default yrange
#   trial = GLMakie.Observable(1) # first trial
# 
#   # toggle buttons
#   toggles = [Toggle(fig, active=active) for active in [false]]
#   toggle_labels = [Label(fig, lift(x -> x ? "$l (on)" : "$l (off)", t.active))
#                    for (t, l) in zip(toggles, ["Events"])]
#   fig[1, 2] = grid!(hcat(toggles, toggle_labels), tellheight=false)
# 
#   # keyboard events
#   on(events(fig).keyboardbutton) do event
#     if event.action in (Keyboard.press, Keyboard.repeat)
#       event.key == Keyboard.right && step_epoch_forward(trial)
#       event.key == Keyboard.left && step_epoch_backward(trial)
#       event.key == Keyboard.down && chans_less(ax, yrange)
#       event.key == Keyboard.up && chans_more(ax, yrange)
#     end
#     return Consume(false)
#   end
# 
#   xlims!(ax, dat.data[1].time[xrange.val[1]], dat.data[1].time[xrange.val[end]])
#   ylims!(ax, yrange.val[1], yrange.val[end])
#   ax.title = "Epoch $(trial.val)/$(length(dat.data))"
#   ax.xlabel = "Time (S)"
#   ax.ylabel = "Amplitude (mV)"
# 
#   # events
#   event_data = @lift(dat.data[$trial][$xrange, [:time, :triggers]].time[dat.data[$trial][$xrange, [:time, :triggers]].triggers.!=0])
#   event_lines = vlines!(event_data, color=:black, linewidth=2)
#   connect!(event_lines.visible, toggles[1].active)
# 
#   function step_epoch_forward(trial::Observable)
#     trial[] = min(length(dat.data), trial.val[1] + 1)
#     ax.title = "Epoch $(trial.val)/$(length(dat.data))"
#   end
# 
#   function step_epoch_backward(trial::Observable)
#     trial[] = max(1, trial.val[1] - 1)
#     ax.title = "Epoch $(trial.val)/$(length(dat.data))"
#   end
# 
#   function chans_less(ax::Axis, yrange::Observable)
#     (yrange.val[1] + 100 >= 0 || yrange.val[end] - 100 <= 0) && return
#     yrange.val = yrange[][1]+100:yrange[][end]-100
#     ylims!(ax, yrange.val[1], yrange.val[end])
#   end
# 
#   function chans_more(ax::Axis, yrange::Observable)
#     yrange.val = yrange[][1]-100:yrange[][end]+100
#     ylims!(ax, yrange.val[1], yrange.val[end])
#   end
# 
#   function draw(ax::Axis, xrange::Observable)
#     for col in names(dat.data[trial.val])
#       if col in channel_labels
#         lines!(ax, @lift(dat.data[$trial][$xrange, 1]), @lift(dat.data[$trial][$xrange, col]), color=:black)
#       end
#     end
#   end
# 
#   draw(ax, xrange)
#   display(fig)
# 
# end

function plot_databrowser(dat::EpochData)
  plot_databrowser(dat, dat.layout.label)
end

function plot_databrowser(dat::EpochData, channel_labels::Union{<:AbstractString,Vector})
  plot_databrowser(dat, [channel_labels])
end



# #################################################################
# plot epoch: Epoched Data
function plot_epochs(dat::EpochData, channels::Union{Vector{<:AbstractString},Vector{Symbol}}; xlim=nothing, ylim=nothing, xlabel="Time (S)", ylabel="mV")

  f = Figure()
  ax = GLMakie.Axis(f[1, 1])

  avg_data = zeros(nrow(dat.data[1]))
  for trial in eachindex(dat.data)
    trial_data = colmeans(dat.data[trial], channels)
    avg_data .+= trial_data
    GLMakie.lines!(dat.data[trial][!, :time], trial_data, color=:grey)
  end
  avg_data ./= length(dat.data)
  GLMakie.lines!(dat.data[1][!, :time], avg_data, color=:black)

  !isnothing(xlim) && xlims!(ax, xlim)
  !isnothing(ylim) && xlims!(ax, ylim)
  ax.xlabel = xlabel
  ax.ylabel = ylabel

  return f
end

function plot_epochs(dat::EpochData, channels::Union{AbstractString,Symbol}; xlim=nothing, ylim=nothing, xlabel="Time (S)", ylabel="mV")
  plot_epochs(dat::EpochData, [channels]; xlim=nothing, ylim=nothing, xlabel="Time (S)", ylabel="mV")
end



