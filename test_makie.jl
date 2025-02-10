using GLMakie





fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, [-10, -2, 0,1,2,3, 6], [0, -5, 0,1,4,9, 16])

# Hacky function to draw cartesian axes by copying the default Makie axes settings
function cartesian_axes!(ax)

    # x-axis major
    lines!(ax, [ax.xaxis.tickvalues.val[1], ax.xaxis.tickvalues.val[end]], [0, 0],
        color=ax.xaxis.attributes.spinecolor)
    for tick in ax.xaxis.tickvalues.val
        lines!(ax, [tick, tick], [-ax.xtickwidth.val / 2, ax.xtickwidth.val / 2],
            color=ax.xtickcolor,
            linewidth=ax.xtickwidth)
        if tick != 0
            text!(ax, tick, -ax.xtickwidth.val, text=string(tick), align=(:center, :top), 
            fontsize = ax.xticklabelsize)
        end
    end

    # y-axis major
    lines!(ax, [0, 0], [ax.yaxis.tickvalues.val[1], ax.yaxis.tickvalues.val[end]], 
        color=ax.yaxis.attributes.spinecolor)
    for tick in ax.yaxis.tickvalues.val
        lines!(ax, [-ax.ytickwidth.val / 2, ax.ytickwidth.val / 2], [tick, tick],
            color=ax.ytickcolor,
            linewidth=ax.xtickwidth)
        if tick != 0
            text!(ax, -ax.ytickwidth.val, tick, text=string(tick), align=(:right, :center), 
            fontsize = ax.yticklabelsize)
        end
    end

    hidedecorations!(ax)
    hidespines!(ax)

end

cartesian_axes!(ax)






