
"""
    read_raw_data(filepath::String; kwargs...) -> Union{BiosemiDataFormat.BiosemiData, BrainVisionDataFormat.BrainVisionData, EuropeanDataFormat.EdfData, FunctionalImageFormat.FifData, ExtensibleDataFormat.XdfData}

Read raw EEG data from various file formats (BDF, BrainVision, EDF, FIF, XDF).

# Arguments
- `filepath::String`: Path to the raw data file.

# Returns
- Raw data object from the underlying reader library.
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
    else
        @minimal_error "Unsupported file extension: $ext"
    end
end
