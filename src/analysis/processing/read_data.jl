
"""
    read_raw_data(filepath::String; kwargs...) -> Union{BiosemiDataFormat.BiosemiData, BrainVisionDataFormat.BrainVisionData}

Read raw EEG data from various file formats (BDF, BrainVision).

# Arguments
- `filepath::String`: Path to the raw data file.

# Returns
- Raw data object from the underlying reader library.
"""
function read_raw_data(filepath::String; kwargs...)
    ext = get_file_extension(filepath)
    if ext == ".bdf"
        return read_bdf(filepath; kwargs...)
    elseif ext == ".vhdr" || ext == ".eeg" || ext == ".vmrk"
        return read_brainvision(filepath; kwargs...)
    else
        @minimal_error "Unsupported file extension: $ext"
    end
end
