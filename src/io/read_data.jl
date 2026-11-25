# Unified raw EEG and saved `.jld2` dataset loading interface.

######################################
########### User Functions ###########
######################################

"""
    read_raw_data(filepath::AbstractString; header_only::Bool = false, kwargs...)

Read raw EEG data from `filepath` (dispatching on `.bdf`, `.vhdr`/`.eeg`/`.vmrk`, `.edf`, `.fif`/`.fiff`, or `.xdf`/`.xdfz`/`.xdf.gz`).
"""
function read_raw_data(filepath::AbstractString; header_only::Bool = false, kwargs...)
    path = String(filepath)
    ext = _get_file_extension(path)
    if ext == ".bdf"
        return read_bdf(path; header_only = header_only, kwargs...)
    elseif ext == ".vhdr" || ext == ".eeg" || ext == ".vmrk"
        return read_brainvision(path; kwargs...)
    elseif ext == ".edf"
        return read_edf(path; kwargs...)
    elseif ext == ".fif" || ext == ".fiff"
        return read_fif(path; kwargs...)
    elseif ext == ".xdf" || ext == ".xdfz" || ext == ".xdf.gz"
        return read_xdf(path; kwargs...)
    end
    @minimal_error "Unsupported file extension: $ext"
end
