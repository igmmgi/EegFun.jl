function test_plot_eog_detection(dat, xlim, channel, detected)
    fig = Figure()
    ax = Axis(fig[1, 1])  # plot layout
    lines!(ax, dat.data.time[xlim], dat.data[!, channel][xlim])
    vlines!(ax, dat.data.time[xlim][dat.data[!, detected][xlim]], color = :black)
    display(fig)
    return fig, ax
end
