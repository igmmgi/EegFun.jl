# === ARTIFACT DETECTION ===
#
# Detection, reporting and repair of artifacts in continuous and epoched EEG data.
# Implementation lives in artifact_detection/:
# - eog_onsets.jl:             EOG onset detection
# - extreme.jl:                Extreme value (absolute threshold) detection
# - step.jl:                   Step/jump detection
# - flatline.jl:               Flatline (low variance) detection
# - peak_to_peak.jl:           Peak-to-peak amplitude detection
# - bridged.jl:                Bridged channel detection
# - epoch_rejection_types.jl:  Rejection info types and their reporting
# - epoch_detection.jl:        Automatic bad epoch detection
# - repair.jl:                 Repairing detected artifacts
# - subset_bad_data.jl:        File-based exclusion of participants with low data retention

include("artifact_detection/eog_onsets.jl")
include("artifact_detection/extreme.jl")
include("artifact_detection/step.jl")
include("artifact_detection/flatline.jl")
include("artifact_detection/peak_to_peak.jl")
include("artifact_detection/bridged.jl")
include("artifact_detection/epoch_rejection_types.jl")
include("artifact_detection/epoch_detection.jl")
include("artifact_detection/repair.jl")
include("artifact_detection/subset_bad_data.jl")
include("artifact_detection/autoreject.jl")
