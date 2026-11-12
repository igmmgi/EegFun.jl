"""
    rename_channel!(dat::EegData, rename_dict::Dict{Symbol,Symbol})

Rename channels in `dat` and its `Layout` in place according to `rename_dict` (`old => new`),
supporting simultaneous swaps (e.g., `Dict(:A => :B, :B => :A)`).
"""
function rename_channel!(dat::EegData, rename_dict::Dict{Symbol,Symbol})
    original_channels = Set{Symbol}(dat.layout.data.label)
    n_found = count(ch -> ch in original_channels, keys(rename_dict))
    if n_found == 0
        rename_channel!(dat.layout, rename_dict)
        @info "rename_channel!: No channels found to rename in data"
        return dat
    end

    _rename_data_columns!(dat, rename_dict, original_channels)
    rename_channel!(dat.layout, rename_dict)
    @info "rename_channel!: Renamed $n_found channels in data and layout"
    return dat
end

function _rename_data_columns!(df::DataFrame, rename_dict::Dict{Symbol,Symbol}, existing_channels::Set{Symbol})
    pairs = Pair{Symbol,Symbol}[]
    targets = Symbol[]
    renamed_cols = Set{Symbol}()
    for (old_name, new_name) in rename_dict
        if old_name in existing_channels && hasproperty(df, old_name)
            push!(targets, new_name)
            push!(renamed_cols, old_name)
            old_name !== new_name && push!(pairs, old_name => new_name)
        end
    end

    if length(targets) != length(unique(targets))
        dup = filter(x -> count(==(x), targets) > 1, unique(targets))
        @minimal_error "Cannot rename channels to duplicate names: $(join(dup, ", "))"
    end

    unrenamed = Set{Symbol}(Symbol(c) for c in names(df) if Symbol(c) ∉ renamed_cols)
    collisions = filter(in(unrenamed), targets)
    if !isempty(collisions)
        @minimal_error "Cannot rename channels to existing column names: $(join(unique(collisions), ", "))"
    end

    isempty(pairs) || rename!(df, pairs)
    return nothing
end

function _rename_data_columns!(dat::EegData, rename_dict::Dict{Symbol,Symbol}, existing_channels::Set{Symbol})
    for container in _data_frames(dat)
        _foreach_df(df -> _rename_data_columns!(df, rename_dict, existing_channels), container)
    end
    return nothing
end

"""
    rename_channel(dat::EegData, rename_dict::Dict{Symbol,Symbol})

Return a copy of `dat` with channels renamed according to `rename_dict`.
"""
rename_channel(dat::EegData, rename_dict::Dict{Symbol,Symbol}) = rename_channel!(copy(dat), rename_dict)
