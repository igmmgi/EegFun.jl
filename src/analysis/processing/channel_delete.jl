"""
    channel_delete!(dat::EegData, channels_to_delete) -> EegData

Delete specified channel(s) (`Symbol`, `String`, or vector thereof) from `dat` in place, updating
both data columns and layout rows (clearing cached layout neighbors if layout channels change).
"""
function channel_delete!(dat::EegData, channels_to_delete::AbstractVector{Symbol})
    to_del = Set(channels_to_delete)
    channels_in_data = intersect(Set(all_labels(dat)), to_del)
    if isempty(channels_in_data)
        @info "channel_delete!: No channels found to delete in data"
        return dat
    end

    channels_in_layout = intersect(Set(dat.layout.data.label), to_del)
    if !isempty(channels_in_layout)
        filter!(:label => ∉(channels_in_layout), dat.layout.data)
        if has_neighbors(dat.layout)
            @info "channel_delete!: Clearing neighbors since channels have changed"
            clear_neighbors!(dat.layout)
        end
    end

    for container in _data_frames(dat)
        _delete_data_columns!(container, channels_in_data)
    end

    @info "channel_delete!: Deleted $(length(channels_in_data)) channel(s): $(join(string.(channels_in_data), ", "))"
    return dat
end
channel_delete!(dat::EegData, channel::Symbol) = channel_delete!(dat, [channel])
channel_delete!(dat::EegData, channel::AbstractString) = channel_delete!(dat, [Symbol(channel)])
channel_delete!(dat::EegData, channels::AbstractVector{<:AbstractString}) = channel_delete!(dat, Symbol.(channels))

function channel_delete!(dat::AbstractVector{<:EegData}, channels_to_delete)
    for d in dat
        channel_delete!(d, channels_to_delete)
    end
    return dat
end

function _delete_data_columns!(df::DataFrame, channels_to_delete::Set{Symbol})
    cols_to_remove = [c for c in propertynames(df) if c in channels_to_delete]
    !isempty(cols_to_remove) && select!(df, Not(cols_to_remove))
    return nothing
end

function _delete_data_columns!(dfs::AbstractVector{<:DataFrame}, channels_to_delete::Set{Symbol})
    isempty(dfs) && return nothing
    cols_to_remove = [c for c in propertynames(dfs[1]) if c in channels_to_delete]
    if !isempty(cols_to_remove)
        for df in dfs
            select!(df, Not(cols_to_remove))
        end
    end
    return nothing
end

"""
    channel_delete(dat, channels_to_delete)

Non-mutating version of [`channel_delete!`](@ref).
"""
channel_delete(dat::EegData, channels_to_delete::AbstractVector{Symbol}) = channel_delete!(copy(dat), channels_to_delete)
channel_delete(dat::EegData, channel::Symbol) = channel_delete(dat, [channel])
channel_delete(dat::EegData, channel::AbstractString) = channel_delete(dat, [Symbol(channel)])
channel_delete(dat::EegData, channels::AbstractVector{<:AbstractString}) = channel_delete(dat, Symbol.(channels))
channel_delete(dat::AbstractVector{<:EegData}, channels_to_delete) = [channel_delete(d, channels_to_delete) for d in dat]
