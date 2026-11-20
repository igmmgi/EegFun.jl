"""
    read_raw_data(filepath::String; header_only::Bool = false, kwargs...)

Read raw EEG data from `filepath` (dispatching on `.bdf`, `.vhdr`/`.eeg`/`.vmrk`, `.edf`, `.fif`/`.fiff`, or `.xdf`/`.xdfz`/`.xdf.gz`).
"""
function read_raw_data(filepath::String; header_only::Bool = false, kwargs...)
    ext = get_file_extension(filepath)
    if ext == ".bdf"
        return read_bdf(filepath; header_only = header_only, kwargs...)
    elseif ext == ".vhdr" || ext == ".eeg" || ext == ".vmrk"
        return read_brainvision(filepath; kwargs...)
    elseif ext == ".edf"
        return read_edf(filepath; kwargs...)
    elseif ext == ".fif" || ext == ".fiff"
        return read_fif(filepath; kwargs...)
    elseif ext == ".xdf" || ext == ".xdfz" || ext == ".xdf.gz"
        return read_xdf(filepath; kwargs...)
    end
    @minimal_error "Unsupported file extension: $ext"
end
