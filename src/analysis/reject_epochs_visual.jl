"""
Visual epoch rejection - analysis module.

This module re-exports the interactive rejection functions from the plotting module
for convenient access. The main interactive GUI is implemented in plot_epoch_rejection.jl.

Use reject_epochs_interactive() to launch the GUI, then extract results using:
- get_rejected_epochs(state) - get indices of rejected epochs
- get_clean_epochs(state) - get cleaned EpochData
- save_rejection_decisions(state, filename) - save report
"""

# The interactive rejection functions are defined in plot_epoch_rejection.jl
# and automatically available when the module is loaded.
# This file exists for organization and potential future extensions.

