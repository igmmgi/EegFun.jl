
"""
    read_raw_data(file_path::String; kwargs...) -> Union{BiosemiDataFormat.BiosemiData, BrainVisionDataFormat.BrainVisionData}

Read raw EEG data from various file formats (BDF, BrainVision).

# Arguments
- `file_path::String`: Path to the raw data file.

# Returns
- Raw data object from the underlying reader library.
"""
function read_raw_data(file_path::String; kwargs...)
    ext = get_file_extension(file_path)
    if ext == ".bdf"
        return read_bdf(file_path; kwargs...)
    elseif ext == ".vhdr" || ext == ".eeg" || ext == ".vmrk"
        return read_brainvision(file_path; kwargs...)
    else
        @minimal_error "Unsupported file extension: $ext"
    end
end
