include("tf.jl")


#####################################################################################
# Test tf_morlet with simulated signal
# Generate 100 trials with multiple frequency components
times, signal = generate_signal(
    100,                             # n_trials
    [-2.0, 3.0],                     # time window
    256.0,                           # sampling rate
    [10.0, 20.0, 40.0],        # frequencies
    [1.0, 1.0, 1.0],            # amplitudes
    [
        [0.25, 0.5],      # time windows
        [0.5, 1.0],
        [1.0, 1.5],
    ],
    0,                              # noise level
)
#lines(times, vec(mean(signal, dims = 2)))
sample_rate = 256;
times = (-2:(1/sample_rate):3-(1/sample_rate));
tf_trials, times_out, freqs_out = tf_morlet(signal, times, 256, 1:1:50, [10]; log_freqs = false)
# tf_trials, times_out, freqs_out = tf_morlet(signal, times, 256, 1:1:50, [4]; log_freqs = false);
tf = average_over_trials(tf_trials);

plot_tf(tf, times_out, freqs_out)


#####################################################################################
# Test tf_morlet with dataFIC from FieldTrip tutorial
signal = DataFrame(CSV.File("dataFIC.csv")).data
signal = reshape(signal, 900, 76)
sample_rate = 300
times = (-1:(1/sample_rate):2-(1/sample_rate))

# TODO: frequency offset?
tf_trials, times_out, freqs_out = tf_morlet(signal, times, 256, 1:1:50, [7], log_freqs = false, tois = -0.5:0.01:1.5)
tf = average_over_trials(tf_trials);
tf = apply_tf_baseline_db(tf, times_out, (-0.5, -0.1));


df = DataFrame(
    time = repeat(times_out, outer = length(freqs_out)),
    frequency = repeat(freqs_out, inner = length(times_out)),
    power = reshape(permutedims(tf), length(times_out) * length(freqs_out)),
)



plot_tf(tf, times_out, freqs_out, log_yscale = false, colorrange = [-2, 1.5], xlim = [-0.5, 1.5], interpolate = false)






#####################################################################################
# Test tf_morlet with dataFIC from book (Figure 13.11)
signal = DataFrame(CSV.File("data2.csv")).data
signal = reshape(signal, 640, 99)
sample_rate = 256
times = (-1:(1/sample_rate):1.5-(1/sample_rate))

# Figure A
tf_trials, times_out, freqs_out = tf_morlet(signal, times, 256, 2:1:80, [3 10]; log_freqs = true);
tf = average_over_trials(tf_trials);
tf = apply_tf_baseline_db(tf, times_out, (-0.5, -0.1));
plot_tf(tf, times_out, freqs_out, log_yscale = true, colorrange = [-3, 3], xlim = [-0.2, 1.0])

# Figure B
tf_trials, times_out, freqs_out = tf_morlet(signal, times, 256, 2:1:80, [3 10]; log_freqs = true);
tf = average_over_trials(tf_trials);
tf = apply_tf_baseline_db(tf, times_out, (-0.5, -0.1));
plot_tf(tf, times_out, freqs_out, log_yscale = false, colorrange = [-3, 3], xlim = [-0.2, 1.0], interpolate = true)



#####################################################################################
# Test tf_hanning with simulated signal
# Generate 100 trials with multiple frequency components
times, signal = generate_signal(
    1000,                             # n_trials
    [-2.0, 3.0],                     # time window
    256.0,                           # sampling rate
    [10.0, 10.0, 20.0, 10.0],        # frequencies
    [4.0, 4.0, 4.0, 4.0],            # amplitudes
    [
        [0.0, 0.25],
        [0.0, 0.25],      # time windows
        [0.0, 0.25],
        [0.0, 0.25],
    ],
    0,                              # noise level
)
#lines(times, vec(mean(signal, dims = 2)))
sample_rate = 256;
times = (-2:(1/sample_rate):3-(1/sample_rate));
tf_trials, times_out, freqs_out = tf_hanning(signal, times, sample_rate, 1:1:50, -0.5:0.2:1.5, window_length = 0.5)
#tf_trials, times_out, freqs_out = tf_hanning(signal, times, sample_rate, 1:1:50, -0.5:0.2:1.5, cycles = 5)
tf = average_over_trials(tf_trials);
plot_tf(tf, times_out, freqs_out)


#####################################################################################
# Test tf_morlet with dataFIC from FieldTrip tutorial
signal = DataFrame(CSV.File("dataFIC.csv")).data
signal = reshape(signal, 900, 76)
sample_rate = 300
times = (-1:(1/sample_rate):2-(1/sample_rate))
#tf_trials, times_out, freqs_out = tf_hanning(signal, times, sample_rate, 1:1:80, -0.5:0.2:1.5, window_length = 0.2)
tf_trials, times_out, freqs_out = tf_hanning(signal, times, sample_rate, 1:1:50, -0.5:0.01:1.5, cycles = 7)
tf = average_over_trials(tf_trials)
tf = apply_tf_baseline_db(tf, times_out, (-0.5, -0.1))

plot_tf(tf, times_out, freqs_out)



####################################################################
signal = DataFrame(CSV.File("data_baseline.csv")).data
sample_rate = 256
times = -1:(1/sample_rate):(1.5-(1/sample_rate))
tf_trials, times_out, freqs_out = tf_morlet(signal, times, sample_rate, 2:1:128, [3 7]; log_freqs = true);
tf = average_over_trials(tf_trials)

# apply different types of baseline normalisations
tf_db_base = apply_tf_baseline_db(tf, times_out, (-0.5, -0.2))
tf_perchange_base = apply_tf_baseline_perchange(tf, times_out, (-0.5, -0.2))
tf_relchange_base = apply_tf_baseline_relchange(tf, times_out, (-0.5, -0.2))
tf_ztransform_base = apply_tf_baseline_ztransforms(tf, times_out, (-0.5, -0.2))

plot_tf(
    tf_db_base,
    times_out,
    freqs_out,
    xlim = [-0.5, 1.5],
    colorrange = [-10, 10],
    log_yscale = true,
    interpolate = true,
)
plot_tf(tf_perchange_base, times_out, freqs_out, xlim = [-0.5, 1.5], colorrange = [-500, 500], log_yscale = true)
plot_tf(tf_relchange_base, times_out, freqs_out, xlim = [-0.5, 1.5], colorrange = [-7.5, 7.5], log_yscale = true)
plot_tf(tf_ztransform_base, times_out, freqs_out, xlim = [-0.5, 1.5], colorrange = [-3.5, 3.5], log_yscale = true)

