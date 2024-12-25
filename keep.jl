###########################################################
mapcols(x -> filtfilt(filter, x), dat.data[:, 3:end])
mapcols!(x -> filtfilt(filter, x), dat.data[:, 3:end])
