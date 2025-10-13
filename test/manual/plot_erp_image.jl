eegfun.plot_erp_image(epochs[1], channel_selection = eegfun.channels([:Fp1]))
eegfun.plot_erp_image(epochs[1], channel_selection = eegfun.channels([:Fp1]), plot_erp = false)

# TODO: x limit consistency between erp image and erp waveform
eegfun.plot_erp_image(epochs[1], layout = :single) 
fig, axes = eegfun.plot_erp_image(epochs[1], channel_selection = eegfun.channels([:Fp1, :Fp2]), layout = :single, boxcar_average = 20, colorrange = (-20, 20)) 
# TODO: no electrode labels in title
eegfun.plot_erp_image(epochs[1], layout = :topo) 
eegfun.plot_erp_image(epochs[1], layout = :grid) 

