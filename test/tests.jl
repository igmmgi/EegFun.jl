###########################################################################
# PLOT EPOCHS
# Basic Tests
# TODO: Implement proper tests
layout = read_layout("./layouts/biosemi72.csv");
dat = read_bdf("../Flank_C_3.bdf");
dat = create_eeg_dataframe(dat, layout);
filter_data!(dat, "hp", "iir", 1, order=1)

# Epoch Data
epoch = extract_epochs(dat, 1, 1, -2, 4)
plot_epochs(epoch, :Fp1) 
plot_epochs(epoch, [:Fp1, :Fp2, :Fpz, :C1, :Cz, :Pz, :O1, :O2, :Oz], kwargs = Dict(:average_channels => true, :ylim => [-100, 100])) 
plot_epochs(epoch, epoch.layout.label, kwargs = Dict(:average_channels => false, :ylim => [-100, 100])) 
plot_epochs(epoch, epoch.layout.label[1:20], kwargs = Dict(:average_channels => false, :ylim => [-100, 100])) 

###########################################################################
# PLOT ERP
# # Basic Tests
# # TODO: Implement proper tests
# layout = read_layout("./layouts/biosemi72.csv");
# dat = read_bdf("../Flank_C_3.bdf");
# dat = create_eeg_dataframe(dat, layout);
# filter_data!(dat, "hp", "iir", 1, order = 1)

# # Epoch Data
# epoch = extract_epochs(dat, 1, 1, -2, 4)

# # ERP Data
# erp = average_epochs(epoch)
# fig, ax = plot_erp(erp, [:Fp1, :Fp2])
# plot_erp(erp, [:Fp1, :Fp2], kwargs = Dict(:average_channels => true))


# fig, ax = plot_erp(erp, [:Fp1, :Fp2], kwargs = Dict(:add_topoplot => false))
# topo_ax = Axis(fig[1,1], 
#                width=Relative(0.2),
#                height=Relative(0.2),
#                halign=0.5,
#                valign=0.5)
# layout = filter(row -> row.label in [:Fp1, :Fp2], erp.layout)
# plot_layout_2d!(fig, topo_ax, layout, 
#                 point_kwargs=Dict(:colormap => :jet, 
#                                  :color => 1:2, 
#                                  :markersize => 18))


plot_erp(erp, erp)


###########################################################################
# PLOT ERP GRID
# # Basic Tests
# # TODO: Implement proper tests
# layout = read_layout("./layouts/biosemi72.csv");
# dat = read_bdf("../Flank_C_3.bdf");
# dat = create_eeg_dataframe(dat, layout);
# filter_data!(dat, "hp", "iir", 1, order = 1)

# # Epoch Data
# epoch = extract_epochs(dat, 1, 1, -2, 4)

# # ERP Data
# erp = average_epochs(epoch)

plot_grid_rect(erp, channels = [:F8, :Fp1, :Fp2, :Fpz, :Fz, :F3, :F4, :F7, :F8])


###########################################################################
# PLOT ERP IMAGE
# Basic Tests
# TODO: Implement proper tests
layout = read_layout("./layouts/biosemi72.csv");
dat = read_bdf("../Flank_C_3.bdf");
dat = create_eeg_dataframe(dat, layout);
filter_data!(dat, "hp", "iir", 1, order=1)

# Epoch Data
epoch = extract_epochs(dat, 1, 1, -2, 4)

plot_erp_image(epoch, :Fp1)
plot_erp_image(epoch, [:Fp1, :IO1] )


###########################################################################
# PLOT TOPO
# Basic Tests
# TODO: Implement proper tests
layout = read_layout("./layouts/biosemi72.csv");
dat = read_bdf("../Flank_C_3.bdf");
dat = create_eeg_dataframe(dat, layout);
filter_data!(dat, "hp", "iir", 1, order = 1)
# Epoch Data
epoch = extract_epochs(dat, 1, 1, -2, 4)
# ERP Data
erp = average_epochs(epoch)
plot_grid_topo(erp)

# Example usage with callback function
function handle_plot_click(channel)
    println("Selected channel: ", channel)
    fig, ax = plot_erp(erp, channel, kwargs = Dict(:average_channels => true))
end

# Create the grid plot with click handling
fig, ax = plot_grid_topo(erp, on_click = handle_plot_click)
# Display the grid plot in a new window
# display(GLMakie.Screen(), fig)

# Alternative example with multiple channel selection
selected_channels = Symbol[]
function handle_multiple_plot_click(channel)
    if channel ∉ selected_channels
        push!(selected_channels, channel)
        println("Selected channels: ", selected_channels)
        # Create a new figure for the ERP plot
        fig, ax = plot_erp(erp, selected_channels)
        # Display in a new window without closing the original
        # display(GLMakie.Screen(), fig)
    end
end

# Create the grid plot with multiple channel selection
fig, pos_map = plot_grid_topo(erp, on_click=handle_multiple_plot_click)
# Display the grid plot in a new window
display(GLMakie.Screen(), fig)

