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
  value::Bool
  label::String
end


function toggle_button_group(fig, labels)
  value_labels = []
  "triggers" in labels && push!(value_labels, ToggleButton(false, "Trigger"))
  "is_vEOG" in labels && push!(value_labels, ToggleButton(false, "vEOG"))
  "is_hEOG" in labels && push!(value_labels, ToggleButton(false, "hEOG"))
  "is_extreme" in labels && push!(value_labels, ToggleButton(false, "extreme"))

  toggle_buttons = [Toggle(fig, active=t.value) for t in value_labels]
  toggle_labels = [Label(fig, t.label, fontsize=30, halign=:left) for t in value_labels]

  return hcat(toggle_buttons, toggle_labels)

end

##################################################################
# Data Browser: Continuous Data
function plot_databrowser(dat::ContinuousData, channel_labels::Vector{<:AbstractString})

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

  # xrange/yrange
  xlimit = 8000
  xrange = GLMakie.Observable(1:xlimit)
  yrange = GLMakie.Observable(-1500:1500) # default yrange
  nchannels = length(channel_labels)
  channel_labels_original = channel_labels
  channels_to_plot = GLMakie.Observable(channel_labels)


  if nchannels > 1 # rough heuristic to space lines across y axis
    offset = GLMakie.Observable(collect(LinRange(1500 * ((nchannels / 100) + 0.2), -1500 * ((nchannels / 100) + 0.2), nchannels)))
  else # just centre
    offset = GLMakie.Observable(0.0)
  end

  @lift xlims!(ax, data.time[$xrange[1]], data.time[$xrange[end]])
  ylims!(ax, yrange.val[1], yrange.val[end])
  ax.xlabel = "Time (S)"
  ax.ylabel = "Amplitude (mV)"

  # toggle buttons for showing events (triggers, vEOG/hEOG, extreme values ...)
  toggles = toggle_button_group(fig, names(data))

  menu = hcat(Menu(fig, options=vcat(["All", "Left", "Right", "Central"], dat.layout.label), default="All"), Label(fig, "Labels"))
  on(menu[1].selection) do s
    channel_labels = [s]
    if s == "All"
      channel_labels = channel_labels_original
    else if s == "Left"
        channel_labels = 
    end
    nchannels = length(channel_labels)
    empty!(ax)
    data = copy(dat.data)
    if nchannels > 1 # rough heuristic to space lines across y axis
      offset = GLMakie.Observable(collect(LinRange(1500 * ((nchannels / 100) + 0.2), -1500 * ((nchannels / 100) + 0.2), nchannels)))
    else # just centre
      offset = GLMakie.Observable(0.0)
    end
    yoffset!()
    draw()
  end

  slider_extreme = Slider(fig[1, 2], range=0:5:200, startvalue=200, width=100)

  fig[1, 2] = grid!(vcat(toggles, hcat(slider_extreme, Label(fig, @lift("Extreme Values: $($(slider_extreme.value)) μV"), fontsize=30)), menu), tellheight=false)
  colsize!(fig.layout, 2, Relative(1 / 8))

  crit_val = lift(slider_extreme.value) do x
    x
  end

  slider_range = Slider(fig[3, 1], range=1000:10000, startvalue=xlimit, snap=false)
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

  #################### Triggers/Events ###############################
  trigger_data_time = @views data[findall(x -> x != 0, data[!, :].triggers), [:time, :triggers]]

  trigger_event_line = []
  trigger_event_label = []
  on(toggles[1].active) do x

    if length(trigger_event_line) >= 1
      delete!(ax, trigger_event_line[1])
      delete!(ax, trigger_event_label[1])
      trigger_event_line = []
      trigger_event_label = []
    else
      xpos = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in trigger_data_time.time]
      push!(trigger_event_line, vlines!(trigger_data_time.time, color=:grey, linewidth=1))
      push!(trigger_event_label, text!(string.(trigger_data_time.triggers), position=xpos, space=:data, align=(:center, :center), fontsize=30))
    end

  end


  if ("is_vEOG" in names(dat.data) && "is_hEOG" in names(dat.data))

    ################### vEOG ###############################
    vEOG_data_time = @views data[findall(x -> x != 0, data[!, :].is_vEOG), [:time, :is_vEOG]]

    vEOG_event_line = []
    vEOG_event_label = []
    on(toggles[2].active) do x
      if length(vEOG_event_line) >= 1
        delete!(ax, vEOG_event_line[1])
        delete!(ax, vEOG_event_label[1])
        vEOG_event_line = []
        vEOG_event_label = []
      else
        xpos = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in vEOG_data_time.time]
        push!(vEOG_event_line, vlines!(vEOG_data_time.time, color=:grey, linewidth=1))
        push!(vEOG_event_label, text!(repeat(["v"], nrow(vEOG_data_time)), position=xpos, space=:data, align=(:center, :center), fontsize=30))
      end
    end

    ################### hEOG ###############################
    hEOG_data_time = @views data[findall(x -> x != 0, data[!, :].is_hEOG), [:time, :is_hEOG]]

    hEOG_event_line = []
    hEOG_event_label = []
    on(toggles[3].active) do x
      if length(hEOG_event_line) >= 1
        delete!(ax, hEOG_event_line[1])
        delete!(ax, hEOG_event_label[1])
        hEOG_event_line = []
        hEOG_event_label = []
      else
        xpos = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in hEOG_data_time.time]
        push!(hEOG_event_line, vlines!(hEOG_data_time.time, color=:grey, linewidth=1))
        push!(hEOG_event_label, text!(repeat(["h"], nrow(hEOG_data_time)), position=xpos, space=:data, align=(:center, :center), fontsize=30))
      end
    end

  end

  if ("is_extreme" in names(data))
    ################### Extreme Values ###############################
    extreme = @views splitgroups(findall(x -> x != 0, data[!, :].is_extreme))

    extreme_event = []
    on(toggles[4].active) do x
      if length(extreme_event) >= 1
        delete!(ax, extreme_event[1])
        extreme_event = []
      else
        push!(extreme_event, vspan!(ax, data[extreme[1], :time], data[extreme[2], :time], color="LightGrey", alpha=0.5))
      end
    end

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

  function draw()
    for col in channel_labels
      lines!(ax, data[!, :time], data[!, col], color=@lift(abs.(dat.data[!, col]) .>= $crit_val), colormap=[:black, :black, :red], linewidth=2)
      text!(ax, @lift(data[$xrange, :time][1]), @lift(data[$xrange, col][1]), text=col, align=(:left, :center), fontsize=20)
    end
  end

  hideydecorations!(ax, label=true)
  yoffset!()
  draw()
  display(fig)
  # DataInspector(fig)

end

function plot_databrowser(dat::ContinuousData)
  plot_databrowser(dat, dat.layout.label)
end

function plot_databrowser(dat::ContinuousData, channel_labels::Union{<:AbstractString,Vector})
  plot_databrowser(dat, [channel_labels])
end



# #################################################################
# Data Browser: Epoched Data
function plot_databrowser(dat::EpochData, channel_labels::Union{Vector{<:AbstractString},Vector{Symbol}})

  fig = Figure()
  ax = GLMakie.Axis(fig[1, 1])  # plot layout

  xrange = GLMakie.Observable(1:nrow(dat.data[1])) # default xrange
  yrange = GLMakie.Observable(-1500:1500) # default yrange
  trial = GLMakie.Observable(1) # first trial

  # toggle buttons
  toggles = [Toggle(fig, active=active) for active in [false]]
  toggle_labels = [Label(fig, lift(x -> x ? "$l (on)" : "$l (off)", t.active))
                   for (t, l) in zip(toggles, ["Events"])]
  fig[1, 2] = grid!(hcat(toggles, toggle_labels), tellheight=false)

  # keyboard events
  on(events(fig).keyboardbutton) do event
    if event.action in (Keyboard.press, Keyboard.repeat)
      event.key == Keyboard.right && step_epoch_forward(trial)
      event.key == Keyboard.left && step_epoch_backward(trial)
      event.key == Keyboard.down && chans_less(ax, yrange)
      event.key == Keyboard.up && chans_more(ax, yrange)
    end
    return Consume(false)
  end

  xlims!(ax, dat.data[1].time[xrange.val[1]], dat.data[1].time[xrange.val[end]])
  ylims!(ax, yrange.val[1], yrange.val[end])
  ax.title = "Epoch $(trial.val)/$(length(dat.data))"
  ax.xlabel = "Time (S)"
  ax.ylabel = "Amplitude (mV)"

  # events
  event_data = @lift(dat.data[$trial][$xrange, [:time, :triggers]].time[dat.data[$trial][$xrange, [:time, :triggers]].triggers.!=0])
  event_lines = vlines!(event_data, color=:black, linewidth=2)
  connect!(event_lines.visible, toggles[1].active)

  function step_epoch_forward(trial::Observable)
    trial[] = min(length(dat.data), trial.val[1] + 1)
    ax.title = "Epoch $(trial.val)/$(length(dat.data))"
  end

  function step_epoch_backward(trial::Observable)
    trial[] = max(1, trial.val[1] - 1)
    ax.title = "Epoch $(trial.val)/$(length(dat.data))"
  end

  function chans_less(ax::Axis, yrange::Observable)
    (yrange.val[1] + 100 >= 0 || yrange.val[end] - 100 <= 0) && return
    yrange.val = yrange[][1]+100:yrange[][end]-100
    ylims!(ax, yrange.val[1], yrange.val[end])
  end

  function chans_more(ax::Axis, yrange::Observable)
    yrange.val = yrange[][1]-100:yrange[][end]+100
    ylims!(ax, yrange.val[1], yrange.val[end])
  end

  function draw(ax::Axis, xrange::Observable)
    for col in names(dat.data[trial.val])
      if col in channel_labels
        lines!(ax, @lift(dat.data[$trial][$xrange, 1]), @lift(dat.data[$trial][$xrange, col]), color=:black)
      end
    end
  end

  draw(ax, xrange)
  display(fig)

end

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



