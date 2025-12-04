#!/usr/bin/env julia
"""
Test script for prepare_permutation_data function.

This script loads ERP data from the specified directory and tests
the prepare_permutation_data function to ensure it works correctly.
"""

using eegfun

# Data directory
input_dir = "/home/ian/Documents/Julia/output_data/filtered_erps_good_lp_30hz"

# File pattern to match
file_pattern = "erps_good"

prepared = eegfun.prepare_permutation_data(
    file_pattern,  # file_pattern (positional)
    1,             # condition_A (positional)
    2;             # condition_B (positional)
    design = :paired,
    input_dir = input_dir,
    participant_selection = eegfun.participants(3:18),  # MATLAB: 3:18
    channel_selection = eegfun.channels(1:72),  # MATLAB: 1:72
    sample_selection = eegfun.samples((0.0, 2.0)),  # MATLAB: [0 1] analysis window
    baseline_window = eegfun.samples((-0.2, 0.0)),  # MATLAB: [-0.2 0] baseline
)
    
# ============================================
# OPTION 1: Cluster-based Monte Carlo (permutation test)
# ============================================
# Equivalent to MATLAB: 'montecarlo', 'cluster'
# Example: analyseMyERPstats(3:18, 1:2, [-0.2 0], [0 1], false, 1:72, 'depsamplesT', 'montecarlo', 'cluster')
result_cluster = eegfun.cluster_permutation_test(
    prepared, 
    n_permutations=1000, 
    min_cluster_size=0, 
    min_num_neighbors=3
)

# ============================================
# OPTION 2: Analytic t-test with NO correction
# ============================================
# Equivalent to MATLAB: 'analytic', 'no'
# Example: analyseMyERPstats(3:18, 1:2, [-0.2 0], [0 1], false, 1:72, 'depsamplesT', 'analytic', 'no')
result_analytic_no = eegfun.analytic_ttest( prepared, alpha=0.05, tail=:both, correction_method=:no)

# ============================================
# OPTION 3: Analytic t-test with BONFERRONI correction
# ============================================
# Equivalent to MATLAB: 'analytic', 'bonferroni'
# Example: analyseMyERPstats(3:18, 1:2, [-0.2 0], [0 1], false, 1:72, 'depsamplesT', 'analytic', 'bonferroni')
result_analytic_bonf = eegfun.analytic_ttest(
    prepared,
    alpha=0.05,
    tail=:both,
    correction_method=:bonferroni
)

# ============================================
# VISUALIZATION: Plot results
# ============================================
# Plot analytic t-test results for a specific channel
# Note: Restart Julia session if function is not found
using GLMakie  # or CairoMakie for static plots

# Try a channel that exists in your data (adjust as needed)
# Check available channels: println(prepared.electrodes)
# Example: Plot t-values with critical t-values
fig = eegfun.plot_analytic_ttest(result_analytic_no, prepared, channel=:PO7, 
                     plot_erp=false, plot_tvalues=true, show_critical_t=true)
display(fig) 

# 1. ERP waveforms only
fig = eegfun.plot_analytic_ttest(result_analytic_no, prepared, channel=:PO7, 
                     plot_erp=true, plot_difference=false, show_significance=false)
display(fig) 

# 2. ERP waveforms + sig points
fig = eegfun.plot_analytic_ttest(result_analytic_no, prepared, channel=:PO7, 
                     plot_erp=true, plot_difference=true, show_significance=true)
display(fig) 

# 3. ERP waveforms + diff waveforms only
fig = eegfun.plot_analytic_ttest(result_analytic_no, prepared, channel=:PO7, 
                     plot_erp=true, plot_difference=true, show_significance=false, show_critical_t=true)
display(fig) 

# 4. ERP waveforms + diff waveforms + sig points
fig = eegfun.plot_analytic_ttest(result_analytic_no, prepared, channel=:PO7, 
                     plot_erp=true, plot_difference=true, show_significance=true)
display(fig) 

# 5. diff waveforms only
fig = eegfun.plot_analytic_ttest(result_analytic_no, prepared, channel=:PO7, 
                     plot_erp=false, plot_difference=true, show_significance=false)
display(fig) 

# 6. diff waveforms + sig points
fig = eegfun.plot_analytic_ttest(result_analytic_no, prepared, channel=:PO7, 
                     plot_erp=false, plot_difference=false, show_significance=true, show_critical_t=true)
display(fig) 

# 6. diff waveforms + sig points
fig = eegfun.plot_analytic_ttest(result_analytic_bonf, prepared, channel=:PO7, 
                     plot_erp=false, plot_difference=false, show_significance=true, show_critical_t=true)
display(fig) 

# 7. Cluster permutation test results (now works with the same function!)
fig = eegfun.plot_analytic_ttest(result_cluster, prepared, channel=:PO7, 
                     plot_erp=false, plot_difference=false, show_significance=true, show_critical_t=true, 
                     sig_bar_position=4.0, sig_bar_color=:red)
display(fig) 