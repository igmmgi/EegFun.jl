# Demo: Interactive ICA Teaching & Geometry GUIs
# Short tutorial script for launching and testing EegFun's 6 interactive
# Independent Component Analysis (ICA) educational simulators.
# Step through each section in VS Code (Alt+Enter / Shift+Enter) or the REPL.

using EegFun
using GLMakie

#######################################################################
# 1. CENTRAL LIMIT THEOREM (signal_example_ica_clt)
#######################################################################

# Controls to test:
# - Number of mixed sources slider & "Regenerate Sources" button
# - Watch non-Gaussian source distributions converge toward a Gaussian bell curve
EegFun.signal_example_ica_clt()

#######################################################################
# 2. LINEAR SOURCE MIXING (signal_example_ica_mixture)
#######################################################################

# Controls to test:
# - Mixing weight slider blending two independent sources into observed mixtures
EegFun.signal_example_ica_mixture()

#######################################################################
# 3. INTERACTIVE SIGNAL UNMIXING (signal_example_ica_separation)
#######################################################################

# Controls to test:
# - Manual unmixing matrix sliders, "Unmix!" button, and "Reset" button
EegFun.signal_example_ica_separation()

#######################################################################
# 4. CENTERING, DECORRELATION & SPHERING (signal_example_ica_sphering)
#######################################################################

# Controls to test:
# - Continuous morph slider (0 -> 3) & stage buttons: Raw -> Centered -> Decorrelated -> Whitened (Sphered)
EegFun.signal_example_ica_sphering()

#######################################################################
# 5. ICA GEOMETRY: INFOMAX VS EXTENDED INFOMAX (signal_example_ica_geometry)
#######################################################################

# Controls to test:
# - Source distribution types (super-Gaussian vs sub-Gaussian), mixing & manual unmixing angles
# - "Infomax" vs "Extended" buttons to compare recovery on sub-Gaussian sources
EegFun.signal_example_ica_geometry()

#######################################################################
# 6. UNMIXING ANGLE OPTIMIZATION LANDSCAPE (signal_example_ica_optimization)
#######################################################################

# Controls to test:
# - Unmixing angle slider (0-180°), "Step Algorithm", and "Run to Convergence"
# - Watch kurtosis / non-Gaussianity objective climb to an independent axis
EegFun.signal_example_ica_optimization()
