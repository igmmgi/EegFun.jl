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

struct Marker
  xpos
  line
  label
end

function toggle_button_group(fig, labels)
  value_labels = []
  "triggers" in labels && push!(value_labels, ToggleButton(false, "Trigger"))
  "is_vEOG" in labels && push!(value_labels, ToggleButton(false, "vEOG"))
  "is_hEOG" in labels && push!(value_labels, ToggleButton(false, "hEOG"))
  "is_extreme" in labels && push!(value_labels, ToggleButton(false, "extreme"))

  toggle_buttons = [Toggle(fig, active=t.value) for t in value_labels]
  toggle_labels = [Label(fig, t.label) for t in value_labels]

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
  xrange = GLMakie.Observable(1:2000)
  yrange = GLMakie.Observable(-1500:1500) # default yrange
  nchannels = length(channel_labels)
  offset = GLMakie.Observable(collect(LinRange(1500 * ((nchannels / 100) + 0.2), -1500 * ((nchannels / 100) + 0.2), nchannels)))

  xlims!(ax, data.time[xrange.val[1]], data.time[xrange.val[end]])
  ylims!(ax, yrange.val[1], yrange.val[end])
  ax.xlabel = "Time (S)"
  ax.ylabel = "Amplitude (mV)"

  # toggle buttons for showing events (triggers, vEOG/hEOG, extreme values ...)
  toggles = toggle_button_group(fig, names(data))
  # menu = hcat(Menu(fig, options=vcat(["All", "Left", "Right", "Central"], dat.layout.label), default="Fp1"), Label(fig, "Labels"))
  fig[1, 2] = grid!(toggles, tellheight=false)
  # fig[1, 2] = grid!(vcat(toggles, menu), tellheight=false)

  # on(menu[1].selection) do s
  #   println(s)
  # end

  # keyboard events
  on(events(fig).keyboardbutton) do event
    if event.action in (Keyboard.press, Keyboard.repeat)
      event.key == Keyboard.left && step_back(ax, xrange)
      event.key == Keyboard.right && step_forward(ax, nrow(dat.data), xrange)
      event.key == Keyboard.down && chans_less(ax, yrange)
      event.key == Keyboard.up && chans_more(ax, yrange)
    end
    return Consume(false)
  end

  ################### Triggers/Events ###############################
  triggers = @lift(findall(x -> x != 0, data[$xrange, :].triggers))
  trigger_data_time = @lift(data[$xrange[$triggers], [:time, :triggers]])
  trigger_event = []

  function plot_event_markers!(event, events, visible, label)
    if length(event) > 0
      for i in eachindex(event)
        delete!(ax, event[i].line)
        delete!(ax, event[i].label)
      end
    end
    empty!(event)
    if visible
      xpos = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in events.val.time]
      if label == "val"
        label = string.(events.val[!, 2])
      else
        label = repeat([label], length(events.val.time))
      end
      push!(event,
        Marker(xpos,
          vlines!(events.val.time, color=:grey, linewidth=1),
          text!(label, position=xpos, space=:data, align=(:center, :center), fontsize=30)))
    end
  end

  on(toggles[1].active) do x
    plot_event_markers!(trigger_event, trigger_data_time, x, "val")
  end

  on(trigger_data_time) do _
    if toggles[1].active.val
      plot_event_markers!(trigger_event, trigger_data_time, toggles[1].active.val, "val")
    end
  end

  if ("is_vEOG" in names(dat.data) && "is_hEOG" in names(dat.data))

    ################### vEOG ###############################
    vEOG = @lift(findall(x -> x != 0, data[$xrange, :].is_vEOG))
    vEOG_data_time = @lift(data[$xrange[$vEOG], [:time, :is_vEOG]])
    vEOG_event = []

    on(toggles[2].active) do x
      plot_event_markers!(vEOG_event, vEOG_data_time, x, "v")
    end

    on(vEOG_data_time) do _
      if toggles[2].active.val
        plot_event_markers!(vEOG_event, vEOG_data_time, toggles[2].active.val, "v")
      end
    end

    ################### hEOG ###############################
    hEOG = @lift(findall(x -> x != 0, data[$xrange, :].is_hEOG))
    hEOG_data_time = @lift(data[$xrange[$hEOG], [:time, :is_hEOG]])
    hEOG_event = []

    on(toggles[3].active) do x
      plot_event_markers!(hEOG_event, hEOG_data_time, x, "h")
    end

    on(hEOG_data_time) do _
      if toggles[3].active.val
        plot_event_markers!(hEOG_event, hEOG_data_time, toggles[3].active.val, "h")
      end
    end

  end

  if ("is_extreme" in names(data))
    ################### Extreme Values ###############################
    extreme = @lift(findall(x -> x != 0, data[$xrange, :].is_extreme))
    extreme_data_time = @lift(data[$xrange[$extreme], [:time, :is_extreme]])
    extreme_event = []

    on(toggles[4].active) do x
      plot_event_markers!(extreme_event, extreme_data_time, x, "")
    end

    on(extreme_data_time) do _
      if toggles[4].active.val
        plot_event_markers!(extreme_event, extreme_data_time, toggles[4].active.val, "")
      end
    end

  end

  function step_back(ax::Axis, xrange::Observable)
    xrange.val[1] - 200 < 1 && return
    xrange[] = xrange.val .- 200
    xlims!(ax, data.time[xrange.val[1]], data.time[xrange.val[end]])
    ylims!(ax, yrange.val[1], yrange.val[end])
  end

  function step_forward(ax::Axis, xmax, xrange::Observable)
    xrange.val[1] + 200 > xmax && return
    xrange[] = xrange.val .+ 200
    xlims!(ax, data.time[xrange.val[1]], data.time[xrange.val[end]])
    ylims!(ax, yrange.val[1], yrange.val[end])
  end

  function chans_less(ax::Axis, yrange::Observable)
    (yrange.val[1] + 100 >= 0 || yrange.val[end] - 100 <= 0) && return
    yrange.val = yrange[][1]+100:yrange[][end]-100
    ylims!(ax, yrange.val[1], yrange.val[end])
    plot_event_markers!(event_data_time, toggles[1].active.val, "val")
    if ("is_vEOG" in names(data) && "is_hEOG" in names(dat.data))
      plot_event_markers!(vEOG_event, vEOG_data_time, toggles[2].active.val, "v")
      plot_event_markers!(hEOG_event, hEOG_data_time, toggles[3].active.val, "h")
    end
  end

  function chans_more(ax::Axis, yrange::Observable)
    yrange.val = yrange[][1]-100:yrange[][end]+100
    ylims!(ax, yrange.val[1], yrange.val[end])
    plot_event_markers!(event_data_time, toggles[1].active.val, "val")
    if ("is_vEOG" in names(data) && "is_hEOG" in names(dat.data))
      plot_event_markers!(vEOG_event, vEOG_data_time, toggles[2].active.val, "v")
      plot_event_markers!(hEOG_event, hEOG_data_time, toggles[3].active.val, "h")
    end

  end

  function yoffset!()
    for (idx, col) in enumerate(channel_labels)
      data[!, col] .+= offset.val[idx]
    end
  end

  function draw(ax, xrange::Observable)
    for col in names(data)
      if col in channel_labels
        # lines!(ax, @lift(data[$xrange, :time]), @lift(data[$xrange, col]), color= :black)
        # lines!(ax, @lift(data[$xrange, :time]), @lift(data[$xrange, col]), color= @lift(data[$xrange, :is_extreme]), colormap = [:black, :red])
        lines!(ax, @lift(data[$xrange, :time]), @lift(data[$xrange, col]), color=@lift(abs.(dat.data[$xrange, col]) .>= 100), colormap=[:black, :black, :red], linewidth=2)
        text!(ax, @lift(data[$xrange, :time][1]), @lift(data[$xrange, col][1]), text=col, align=(:left, :center), fontsize=20)
      end
    end
  end

  hideydecorations!(ax, label=true)
  yoffset!()
  draw(ax, xrange)
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



