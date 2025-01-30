using DimensionalData

# Create a custom sampled dimension with both time and sample indices
time_coords = range(-1, step=1/300, length=size(signal_2d, 1))
sample_coords = 1:size(signal_2d, 1)
time_dim = Sampled(1:size(signal_2d, 1), 
                   (time=time_coords, samples=sample_coords),
                   :time)
trial_dim = Dim(1:n_trials, :trial)

# Create DimArray with the multi-coordinate dimension
ds = DimArray(signal_2d, (time_dim, trial_dim); name="signal")

# Now you can access either coordinate
ds[Ti(time=0.5)]  # Access by time
ds[Ti(samples=150)]  # Access by sample number

# Or get the coordinates
ds[Ti].val.time  # Get time values
ds[Ti].val.samples  # Get sample indices