using EegFun

# Suppose we have ERPs from three participants
# results = [erp_p1, erp_p2, erp_p3]
# For this demo, let's assume 'results' is already defined

# Create grand averages for all conditions found in the list
# grand_avgs = grand_average(results)

# --- Batch Averaging (from Disk) ---

# Averages all JLD2 files in the folder that match the "erps_cleaned" pattern
# This creates a new directory "grand_average_erps_cleaned" containing the results
# grand_average("erps_cleaned", input_dir = "derivatives/erp_analysis/")

# Filtering Participants and Conditions
# grand_average("erps_cleaned",
#     participant_selection = participants(1:10),  # Only first 10 participants
#     condition_selection = conditions([1, 2]),    # Only "Standard" and "Target"
#     output_dir = "results/grand_averages"
# )

# --- Visualizing the Result ---

# Load the batch results
# ga_results = read_data("results/grand_averages/grand_average_erps_cleaned.jld2")

# Plot with standard error shading (if n > 1)
# plot_erp(ga_results, labels = ["Standard", "Target"])
