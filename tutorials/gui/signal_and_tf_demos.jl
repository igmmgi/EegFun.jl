# Demo: Interactive Signal Processing & Time-Frequency GUIs
# Short tutorial script for launching and testing EegFun's 5 interactive
# signal processing, Fourier, convolution, and time-frequency simulators.
# Step through each section in VS Code (Alt+Enter / Shift+Enter) or the REPL.

using EegFun
using GLMakie

#######################################################################
# 1. SAMPLING & NYQUIST ALIASING (signal_example_sampling)
#######################################################################

# Controls to test:
# - True signal frequency vs sampling rate sliders
# - Observe Nyquist limit and reconstructed aliased sine wave when fs < 2*f
EegFun.signal_example_sampling()

#######################################################################
# 2. SIGNAL COMPOSITION, FFT & FILTERING (signal_example_composition)
#######################################################################

# Controls to test:
# - Enable/adjust multiple sine components (frequency, amplitude, phase) + white noise
# - Inspect composite time-domain waveform, FFT amplitude spectrum, and live filter effects
EegFun.signal_example_composition()

#######################################################################
# 3. DOT PRODUCT & COMPLEX SINE WAVES (signal_example_dotproduct)
#######################################################################

# Controls to test:
# - Probe frequency & phase sliders vs target signal
# - Toggle real sine vs Euler complex wavelet (sin + i*cos) to see phase-invariant magnitude
EegFun.signal_example_dotproduct()

#######################################################################
# 4. SLIDING CONVOLUTION & WAVELETS (signal_example_convolution)
#######################################################################

# Controls to test:
# - Kernel menu (Gaussian, Boxcar, Morlet Wavelet), kernel width/frequency sliders
# - Manual position slider or "▶ Play Animation", plus preset buttons
EegFun.signal_example_convolution()

#######################################################################
# 5. TIME-FREQUENCY METHODS & BASELINES (signal_example_tf)
#######################################################################

# Controls to test:
# - Method menu: Morlet Wavelet vs STFT (Hanning) vs Multitaper (DPSS)
# - Cycles/window parameter slider, noise slider, linear vs log10 frequency toggle
# - Baseline normalization menu (None, dB, % Change, Z-score, Relative, etc.) & window IntervalSlider
# - Up to 3 independent oscillatory components with onset-offset IntervalSliders
EegFun.signal_example_tf()
