# Demo: Interactive Psychophysiology GUIs (hEOG, EMG, ECG, HRV & EDA)
# Short tutorial script for launching and testing EegFun's 6 interactive
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

#######################################################################
# 5. EDA SWEAT GLAND BIOPHYSICS & SCR MEASURES (signal_example_eda_scr)
#######################################################################

# Controls to test:
# - Baseline SCL, SCR amplitude, latency, rise time, and 50% half-recovery (Rec t/2) sliders
# - "▶ Trigger Stimulus" button (animates sweat duct poral valve filling & SCR waveform)
# - Low vs. High Tonic SCL presets
EegFun.signal_example_eda_scr()

#######################################################################
# 6. EDA HABITUATION, TONIC/PHASIC & OVERLAP (signal_example_eda_habituation)
#######################################################################

# Controls to test:
# - Initial SCR (S1), habituation rate, stimulus ISI, tonic drift, and spontaneous fluctuations (SFs/min)
# - Shorten ISI (< 7 s) or click "Fast ISI Overlap (Superposition)" to see response pile-up
# - Toggle "Novel Stimulus at S6 (Dishabituation)" or "Non-Responder (~10%)"
EegFun.signal_example_eda_habituation()
