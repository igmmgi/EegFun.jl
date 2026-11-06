# Demo: Interactive ERP Simulation & MVPA Decoding GUIs
# Short tutorial script for launching and testing EegFun's interactive
# ERP trial-averaging simulator and MVPA linear decoding simulator.
# Step through each section in VS Code (Alt+Enter / Shift+Enter) or the REPL.

using EegFun
using GLMakie

#######################################################################
# 1. ERP TRIAL-AVERAGING & JITTER SIMULATOR (simulate_erp)
#######################################################################

# Controls to test:
# - Number of trials (1-500) & single-trial noise amplitude sliders (watch 1/sqrt(N) SNR gain)
# - "Regenerate" button & confidence interval toggle
# - Up to 3 ERP peak components: Toggle active, frequency, amplitude, peak latency,
#   amplitude jitter, and latency jitter (watch latency jitter smear & attenuate peaks!)
EegFun.simulate_erp()

#######################################################################
# 2. MVPA LINEAR CLASSIFIER DECODING DEMO (signal_example_decoding)
#######################################################################

# Controls to test:
# - Condition effect size, trial count, noise level, and cross-validation controls
# - "Decode!" button to inspect time-resolved classification accuracy vs chance
EegFun.signal_example_decoding()
