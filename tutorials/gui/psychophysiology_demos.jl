# Demo: Interactive Psychophysiology GUIs (hEOG, EMG, ECG & HRV)
# Short tutorial script for launching and testing EegFun's 4 interactive
# peripheral physiology simulators.
# Step through each section in VS Code (Alt+Enter / Shift+Enter) or the REPL.

using EegFun
using GLMakie

#######################################################################
# 1. HORIZONTAL EOG DIPOLE & SACCADE SIMULATOR (signal_example_heog)
#######################################################################

# Controls to test:
# - Gaze angle (-45° to +45°), saccade duration, and common-mode noise sliders
# - "Play Saccade" animation button
# - Electric field lines checkbox & monopolar vs bipolar (CMR) derivation checkbox
EegFun.signal_example_heog()

#######################################################################
# 2. SURFACE EMG RECRUITMENT & ENVELOPE CASCADE (signal_example_emg)
#######################################################################

# Controls to test:
# - Contraction force, motor-unit recruitment, and filter/smoothing sliders
# - "▶ Play Contraction" button showing raw sEMG -> full-wave rectified -> low-pass/RMS -> integrated EMG
EegFun.signal_example_emg()

#######################################################################
# 3. CARDIAC DIPOLE & EINTHOVEN'S TRIANGLE (signal_example_ecg_einthoven)
#######################################################################

# Controls to test:
# - Heart rate, electrical heart axis angle, and rhythm/waveform controls
# - Live rotating cardiac vector projected onto Limb Leads I, II, and III
EegFun.signal_example_ecg_einthoven()

#######################################################################
# 4. HEART RATE VARIABILITY & AUTONOMIC BALANCE (signal_example_ecg_hrv)
#######################################################################

# Controls to test:
# - Sympathetic vs parasympathetic (vagal/RSA) tone, respiration rate, and heart rate sliders
# - Inspect beat-to-beat R-R tachogram, time-domain RMSSD/SDNN, and LF/HF spectrum
EegFun.signal_example_ecg_hrv()
