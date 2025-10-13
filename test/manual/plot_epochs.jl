eegfun.plot_epochs(epochs[1], channel_selection = eegfun.channels([:Fp1, :Fp2]))
eegfun.plot_epochs(epochs[1])


eegfun.plot_epochs(epochs[1], layout = :grid)
eegfun.plot_epochs(epochs[1])

eegfun.plot_epochs(epochs[1], layout = :topo)





# These work exactly as before
eegfun.plot_epochs(epochs[1], channel_selection = eegfun.channels([:Fp1, :Fp2, :Cz]), layout_type = :topo, layout_kwargs = Dict(:plot_width => 0.15))

# But now you can easily control layout
eegfun.plot_epochs(epochs[1], channel_selection = eegfun.channels([:Fp1, :Fp2, :Cz]); 
            layout_type = :topo, layout_kwargs = Dict(:plot_width => 0.15))


fig, ax = eegfun.plot_epochs(epochs[1])
fig, ax = eegfun.plot_epochs(epochs[1], layout = :grid)
fig, ax = eegfun.plot_epochs(epochs[1], layout = :topo)


# And the system automatically chooses the best option
plot_epochs(epochs, channel_selection = channels([:Fp1]))  # Single plot
plot_epochs(epochs, channel_selection = channels([:Fp1, :Fp2, :Cz]))  # Grid

plot_epochs(epochs[1], :Cz)

