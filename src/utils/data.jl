#  === EEG DATA UTILITIES ===
#
# This file provides utilities for accessing EegFun data structures.
# It includes functions for column identification, data access, convenience functions,
# data viewing, mathematical utilities, and subsetting operations.
#
# Organization (implementation lives in utils/data/):
# - columns.jl:        Column identification (metadata, channel, extra columns)
# - access.jl:         Data access, convenience functions, data viewing
# - math.jl:           Mathematical utilities (_datarange, _colmeans, data_limits, y-ranges)
# - subset.jl:         Subsetting EEG data objects, selection helpers, epoch_to_continuous
# - subset_files.jl:   File-based (batch) subsetting
# - rename.jl:         Renaming channels
# - from_raw.jl:       create_eegfun_data constructors from raw readers
# - base_interface.jl: Base interface extensions

include("data/columns.jl")
include("data/access.jl")
include("data/math.jl")
include("data/subset.jl")
include("data/subset_files.jl")
include("data/rename.jl")
include("data/from_raw.jl")
include("data/base_interface.jl")
