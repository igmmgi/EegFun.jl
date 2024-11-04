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
function plot_databrowser(dat::ContinuousData, channel_labels::Union{Vector{<:AbstractString},Vector{Symbol}})

  fig = Figure()
  ax = GLMakie.Axis(fig[1, 1])  # plot layout

  xrange = GLMakie.Observable(1:2000) # default xrange
  yrange = GLMakie.Observable(-1500:1500) # default yrange

  # toggle buttons
  toggles = [Toggle(fig, active=active) for active in [false]]
  toggle_labels = [Label(fig, lift(x -> x ? "$l (on)" : "$l (off)", t.active))
                   for (t, l) in zip(toggles, ["Events"])]
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

  xlims!(ax, dat.data.time[xrange.val[1]], dat.data.time[xrange.val[end]])
  ylims!(ax, yrange.val[1], yrange.val[end])
  ax.xlabel = "Time (S)"
  ax.ylabel = "Amplitude (mV)"

  # events
  event_data_time = @lift(dat.data[$xrange, [:time, :events]].time[dat.data[$xrange, [:time, :events]].events.!=0])
  event_data_trigger = @lift(dat.data[$xrange, [:time, :events]].events[dat.data[$xrange, [:time, :events]].events.!=0])

  event_lines = vlines!(event_data_time, color=:grey, linewidth=2)
  text_lines = text!(string.(event_data_trigger), position=event_data_time)

  connect!(text_lines.visible, toggles[1].active)
  connect!(event_lines.visible, toggles[1].active)


  function step_back(ax::Axis, xrange::Observable)
    xrange.val[1] - 100 < 1 && return
    xrange[] = xrange.val .- 100
    xlims!(ax, dat.data.time[xrange.val[1]], dat.data.time[xrange.val[end]])
    ylims!(ax, yrange.val[1], yrange.val[end])
  end

  function step_forward(ax::Axis, xmax, xrange::Observable)
    println(event_data_time)
    println(event_data_trigger)
    xrange.val[1] + 100 > xmax && return
    xrange[] = xrange.val .+ 100
    xlims!(ax, dat.data.time[xrange.val[1]], dat.data.time[xrange.val[end]])
    ylims!(ax, yrange.val[1], yrange.val[end])
  end

  function chans_less(ax::Axis, yrange::Observable)
    (yrange.val[1] + 100 >= 0 || yrange.val[end] - 100 <= 0) && return
    yrange.val = yrange[][1]+100:yrange[][end]-100
    xlims!(ax, dat.data.time[xrange.val[1]], dat.data.time[xrange.val[end]])
    ylims!(ax, yrange.val[1], yrange.val[end])
  end

  function chans_more(ax::Axis, yrange::Observable)
    yrange.val = yrange[][1]-100:yrange[][end]+100
    xlims!(ax, dat.data.time[xrange.val[1]], dat.data.time[xrange.val[end]])
    ylims!(ax, yrange.val[1], yrange.val[end])
  end

  function draw(ax::Axis, xrange::Observable)
    for col in names(dat.data)
      if col in channel_labels
        lines!(ax, @lift(dat.data[$xrange, 1]), @lift(dat.data[$xrange, col]), color=:black)
      end
    end
  end

  draw(ax, xrange)
  display(fig)

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
  event_data = @lift(dat.data[$trial][$xrange, [:time, :events]].time[dat.data[$trial][$xrange, [:time, :events]].events.!=0])
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



