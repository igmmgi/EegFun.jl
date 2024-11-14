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


##################################################################
# Data Browser: Continuous Data
function plot_databrowser(dat::ContinuousData, channel_labels::Vector{<:AbstractString})

  data = copy(dat.data)

  fig = Figure()
  ax = GLMakie.Axis(fig[1, 1])

  # controls
  # interactions(ax)
  # deregister_interaction!(ax, :rectanglezoom)
  # deregister_interaction!(ax, :dragpan)
  # deregister_interaction!(ax, :scrollzoom)
  # deregister_interaction!(ax, :limitreset)

  # xrange/yrange
  nchannels = length(channel_labels)
  bounds = 1500
  xrange = GLMakie.Observable(1:2000)
  yrange = GLMakie.Observable(-bounds:bounds) # default yrange
  offset = GLMakie.Observable(collect(LinRange(bounds * ((nchannels / 100) + 0.2), -bounds * ((nchannels / 100) + 0.2), nchannels)))

  xlims!(ax, data.time[xrange.val[1]], data.time[xrange.val[end]])
  ylims!(ax, yrange.val[1], yrange.val[end])
  ax.xlabel = "Time (S)"
  ax.ylabel = "Amplitude (mV)"

  scale = GLMakie.Observable(1.0)


  # toggle buttons for showing events (triggers, vEOG/hEOG, extreme values ...)
  toggle_values = []
  toggle_labels = []
  if "triggers" in names(data)
    push!(toggle_values, false)
    push!(toggle_labels, "Trigger")
  end
  if "is_vEOG" in names(data)
    push!(toggle_values, false)
    push!(toggle_labels, "vEOG")
  end
  if "is_hEOG" in names(data)
    push!(toggle_values, false)
    push!(toggle_labels, "hEOG")
  end
  if "is_extreme" in names(data)
    push!(toggle_values, false)
    push!(toggle_labels, "Extreme Value")
  end

  # is = IntervalSlider(fig[2, 1], range=@lift(data[$xrange, :time]), startvalues=(1, 2))
  # vl = lift(x -> x[1], is.interval)
  # vr = lift(x -> x[2], is.interval)
  # vspan!(fig[1, 1], vl, vr)

  toggles = [Toggle(fig, active=active) for active in toggle_values]
  toggle_labels = [Label(fig, lift(x -> x ? "$l (on)" : "$l (off)", t.active))
                   for (t, l) in zip(toggles, toggle_labels)]

  # menu = hcat(Menu(fig, options = dat.layout.label, default = "Fp1"), Label(fig, "menu"))
  # fig[1, 2] = grid!(vcat(hcat(toggles, toggle_labels), menu), tellheight=false)
  fig[1, 2] = grid!(hcat(toggles, toggle_labels), tellheight=false)

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
  event_data_time = @lift(data[$xrange[$triggers], [:time, :triggers]])
  event_lines = []
  event_text = []

  function plot_events(event_data_time, visible)
    if length(event_lines) > 0
      for i in eachindex(event_lines)
        delete!(ax, event_lines[i])
      end
    end
    if length(event_text) > 0
      for i in eachindex(event_text)
        delete!(ax, event_text[i])
      end
    end
    if visible
      push!(event_lines, vlines!(event_data_time.val.time, color=:grey, linewidth=1))
      event_label_pos = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in event_data_time.val.time]
      push!(event_text, text!(string.(event_data_time.val.triggers), position=event_label_pos, space=:data,
        align=(:center, :center), fontsize=30))
    end
  end

  on(toggles[1].active) do x
    plot_events(event_data_time, x)
  end

  on(event_data_time) do _
    if toggles[1].active.val
      plot_events(event_data_time, toggles[1].active.val)
    end
  end

  if ("is_vEOG" in names(dat.data) && "is_hEOG" in names(dat.data))

    ################### vEOG ###############################
    vEOG = @lift(findall(x -> x != 0, data[$xrange, :].is_vEOG))
    vEOG_data_time = @lift(data[$xrange[$vEOG], [:time, :is_vEOG]])
    vEOG_lines = []
    vEOG_text = []

    function plot_vEOG(vEOG_data_time, visible)
      if length(vEOG_lines) > 0
        for i in eachindex(vEOG_lines)
          delete!(ax, vEOG_lines[i])
        end
      end
      if length(vEOG_text) > 0
        for i in eachindex(vEOG_text)
          delete!(ax, vEOG_text[i])
        end
      end
      if visible
        push!(vEOG_lines, vlines!(vEOG_data_time.val.time, color=:grey, linewidth=1))
        vEOG_label_pos = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in vEOG_data_time.val.time]
        push!(vEOG_text, text!(repeat(["V"], length(vEOG_data_time.val.time)), position=vEOG_label_pos, space=:data,
          align=(:center, :center), fontsize=30))
      end
    end

    on(toggles[2].active) do x
      plot_vEOG(vEOG_data_time, x)
    end

    on(vEOG_data_time) do _
      if toggles[2].active.val
        plot_vEOG(vEOG_data_time, toggles[2].active.val)
      end
    end

    ################### hEOG ###############################
    hEOG = @lift(findall(x -> x != 0, data[$xrange, :].is_hEOG))
    hEOG_data_time = @lift(data[$xrange[$hEOG], [:time, :is_hEOG]])
    hEOG_lines = []
    hEOG_text = []

    function plot_hEOG(hEOG_data_time, visible)
      if length(hEOG_lines) > 0
        for i in eachindex(hEOG_lines)
          delete!(ax, hEOG_lines[i])
        end
      end
      if length(hEOG_text) > 0
        for i in eachindex(hEOG_text)
          delete!(ax, hEOG_text[i])
        end
      end
      if visible
        push!(hEOG_lines, vlines!(hEOG_data_time.val.time, color=:grey, linewidth=1))
        hEOG_label_pos = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in hEOG_data_time.val.time]
        push!(hEOG_text, text!(repeat(["H"], length(hEOG_data_time.val.time)), position=hEOG_label_pos, space=:data,
          align=(:center, :center), fontsize=30))
      end
    end

    on(toggles[3].active) do x
      plot_hEOG(hEOG_data_time, x)
    end

    on(hEOG_data_time) do _
      if toggles[3].active.val
        plot_hEOG(hEOG_data_time, toggles[3].active.val)
      end
    end

  end

  if ("is_extreme" in names(data))
    ################### Extreme Values ###############################
    extreme = @lift(findall(x -> x != 0, data[$xrange, :].is_extreme))
    extreme_data_time = @lift(data[$xrange[$extreme], [:time, :is_extreme]])
    extreme_lines = []

    function plot_extreme(extreme_data_time, visible)
      if length(extreme_lines) > 0
        for i in eachindex(extreme_lines)
          delete!(ax, extreme_lines[i])
        end
      end
      if visible
        push!(extreme_lines, vlines!(extreme_data_time.val.time, color=:grey, linewidth=1))
      end
    end

    on(toggles[4].active) do x
      plot_extreme(extreme_data_time, x)
    end

    on(extreme_data_time) do _
      if toggles[4].active.val
        plot_extreme(extreme_data_time, toggles[4].active.val)
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
    #(yrange.val[1] + 100 >= 0 || yrange.val[end] - 100 <= 0) && return
    #yrange.val = yrange[][1]+100:yrange[][end]-100
    #ylims!(ax, yrange.val[1], yrange.val[end])
    scale[] = scale.val * 0.9
    # plot_events(event_data_time, toggles[1].active.val)
    # if ("is_vEOG" in names(data) && "is_hEOG" in names(dat.data))
    #   plot_vEOG(vEOG_data_time, toggles[2].active.val)
    #   plot_hEOG(hEOG_data_time, toggles[3].active.val)
    # end
  end

  function chans_more(ax::Axis, yrange::Observable)
    #yrange.val = yrange[][1]-100:yrange[][end]+100
    #ylims!(ax, yrange.val[1], yrange.val[end])
    scale[] = scale.val * 1.1
    # plot_events(event_data_time, toggles[1].active.val)
    # if ("is_vEOG" in names(data) && "is_hEOG" in names(dat.data))
    #   plot_vEOG(vEOG_data_time, toggles[2].active.val)
    #   plot_hEOG(hEOG_data_time, toggles[3].active.val)
    # end
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
        lines!(ax, @lift(data[$xrange, :time]), @lift(data[$xrange, col] .* $scale), color=@lift(abs.(dat.data[$xrange, col]) .>= 100), colormap=[:black, :black, :red], linewidth=2)
        # text!(ax, @lift(data[$xrange, :time][1]), @lift(data[$xrange, col][1] .* scale.val[1]), text=col, align=(:left, :center), fontsize=20)
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



