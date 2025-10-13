

# TODO

# ERP Plot
fig, ax = eegfun.plot_erp(erps, average_channels = true, title = "Custom Title")
display(fig)


# ERP Plot
fig, ax = eegfun.plot_erp(erps, channel_selection = eegfun.channels([:Cz, :PO7, :PO8]), average_channels = true)
display(fig)


# ERP Plot
fig, ax = eegfun.plot_erp(erps[1], layout = [2, 2])
display(fig)


fig, ax = eegfun.plot_erp([erps[1], copy(erps[1])])
display(fig)

fig, ax = eegfun.plot_erp(erps[1], layout = :grid)
display(fig)

fig, ax = eegfun.plot_erp(erps, layout = :grid, title = "Custom Title")
display(fig)

fig, ax = eegfun.plot_erp(erps[1], layout = :topo)
display(fig)

fig, ax = eegfun.plot_erp(erps, layout = :topo, channel_selection = eegfun.channels([:Fp1, :Fp2, :PO8]))
display(fig)

fig, ax = eegfun.plot_erp(erps, layout = :grid, channel_selection = x -> startswith.(string.(x), "F"))
display(fig)




fig, ax = eegfun.plot_erp([erps[1], copy(erps[1])], channel_selection = eegfun.channels([:Fp1, :Fp2]), layout = :grid)
display(fig)

fig, ax = eegfun.plot_erp(erps[1], linewidth = 2)
fig, ax = eegfun.plot_erp(erps[1], layout = :grid)
fig, ax = eegfun.plot_erp(erps[1], layout = :topo)





# ERP Plot
fig, ax = eegfun.plot_erp(erps[1], channel_selection = eegfun.channels([:Fp1, :Fp2]))
display(fig)


fig, ax = eegfun.plot_erp(erps[1], channel_selection = eegfun.channels([:Fp1, :Fp2]), layout = :grid)
display(fig)

fig, ax = eegfun.plot_erp(erps[1], channel_selection = eegfun.channels([:Fp1, :Fp2, :O1]), layout = :topo)
display(fig)








# ERP Plot
fig, ax = eegfun.plot_erp(erps)
display(fig)

fig, ax = eegfun.plot_erp(erps[1]; layout = :grid)
display(fig)

fig, ax = eegfun.plot_erp(erps[1]; layout = :topo)
display(fig)




# ERP Plot
fig, ax = eegfun.plot_erp(erps[1]; channel_selection = eegfun.channels([:Fp1, :Fp2]), sample_selection = x -> -1 .< x.time .< 1.5)
display(fig)

# ERP Plot
@time fig, ax = eegfun.plot_erp(erps; channel_selection = eegfun.channels([:Fp1, :P08]), sample_selection = x -> -1 .< x.time .< 1.5); display(fig);

fig, ax = eegfun.plot_erp(erps); display(fig);


# ERP Plot
plot_erp(erps[1])
plot_erp(erps[1], :Fp1)
plot_erp(erps[1], [:Fp1, :Fp2])
plot_erp(erps[2], [:Fp1, :Fp2, :Cz])
plot_erp(erps[1], [:Fp1, :Fp2], kwargs = Dict(:average_channels => true))
plot_erp(erps[1], [:Fp1, :Fp2], kwargs = Dict(:average_channels => true))
plot_erp(erps[1], [:Fp1, :Fp2], kwargs = Dict(:average_channels => false, :add_topoplot => true))
plot_erp(erps[1], erps[3], [:PO7, :PO8])
plot_erp([erps[1], erps[2], erps[3]], [:PO7, :PO8], kwargs = Dict(:average_channels => true, :add_topoplot => true))
