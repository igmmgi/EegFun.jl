# Demo: Response-Locked Realignment
# Shows how to realign stimulus-locked epochs to a different time point
# (e.g., response time) for response-locked ERP analysis.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
using GLMakie

const DEMO_OUTPUT = "./tutorials/output/"
mkpath(DEMO_OUTPUT)

#######################################################################
# LOAD EPOCHED DATA
#######################################################################

# Load stimulus-locked epochs for a single participant
epochs = EegFun.read_data(EegFun.example_path("data/julia/epochs/example1_epochs.jld2"))

#######################################################################
# REALIGN EPOCHS (IN-MEMORY)
#######################################################################

# Realign epochs to response triggers (e.g., 201, 202). 
# This shifts the time axis of each epoch so that the response trigger becomes t=0,
# and drops epochs that do not contain any of these triggers.
EegFun.realign!(epochs, [201, 202])

# The time axis is now response-locked! Let's compute the ERP.
# We will use condition 1.
# Note that we use grand_average since we are taking the mean over epochs
response_locked_erp = EegFun.grand_average([epochs[1]])

# Plot the response-locked ERP
EegFun.plot_erp(response_locked_erp)

#######################################################################
# TRIGGER INTERVALS (REACTION TIMES)
#######################################################################

# We can also calculate the time interval between the stimulus and the response 
# (i.e. the reaction time) and append it as a column to the epoch data.
epochs_rt = EegFun.read_data(EegFun.example_path("data/julia/epochs/example1_epochs.jld2"))

# Calculate the interval between stimulus (e.g., 101) and response (e.g., 201)
EegFun.calculate_trigger_interval!(epochs_rt, [101], [201], column_name = :reaction_time)


#######################################################################
# BATCH REALIGNMENT (FROM DISK)
#######################################################################

# Realign all epoch files matching the pattern to response triggers
EegFun.realign(
    "epochs",
    [201, 202],
    input_dir = EegFun.example_path("data/julia/epochs/"),
    output_dir = DEMO_OUTPUT
)
