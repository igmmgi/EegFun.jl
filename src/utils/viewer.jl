# Interactive tabular `DataFrame` and `EegData` inspection helpers.

######################################
########### User Functions ###########
######################################

"""
    viewer(dat)

Display `dat` (or `all_data(dat)` when defined) using `Main.vscodedisplay` if available, falling back to `display`.
"""
function viewer(dat)
    target = applicable(all_data, dat) ? all_data(dat) : dat
    if isdefined(Main, :vscodedisplay)
        try
            Main.vscodedisplay(target)
        catch e
            @debug "vscodedisplay failed, falling back to display()" exception = e
            display(target)
        end
    else
        display(target)
    end
end

"""
    head(dat::EegData; n = 5)

Display and return a view of the first `n` rows of `all_data(dat)`.
"""
function head(dat::EegData; n = nothing)
    n_rows_show = isnothing(n) ? 5 : n
    data = _head_source_df(dat, n_rows_show)
    k = clamp(n_rows_show, 0, nrow(data))
    result = view(data, 1:k, :)
    viewer(result)
    return result
end

"""
    tail(dat::EegData; n = 5)

Display and return a view of the last `n` rows of `all_data(dat)`.
"""
function tail(dat::EegData; n = nothing)
    n_rows_show = isnothing(n) ? 5 : n
    data = _tail_source_df(dat, n_rows_show)
    nrows = nrow(data)
    k = clamp(n_rows_show, 0, nrows)
    result = view(data, (nrows-k+1):nrows, :)
    viewer(result)
    return result
end

######################################
######### Internal Functions #########
######################################

function _head_source_df(dat::EegData, n_rows_show::Int)::DataFrame
    return all_data(dat)
end

function _head_source_df(dat::MultiDataFrameEeg, n_rows_show::Int)::DataFrame
    isempty(dat.data) && return DataFrame()
    (n_rows_show <= 0 || nrow(dat.data[1]) >= n_rows_show || length(dat.data) == 1) && return dat.data[1]
    acc = 0
    last_ep = length(dat.data)
    for i in eachindex(dat.data)
        acc += nrow(dat.data[i])
        if acc >= n_rows_show
            last_ep = i
            break
        end
    end
    return reduce(vcat, @view(dat.data[1:last_ep]))
end

function _tail_source_df(dat::EegData, n_rows_show::Int)::DataFrame
    return all_data(dat)
end

function _tail_source_df(dat::MultiDataFrameEeg, n_rows_show::Int)::DataFrame
    isempty(dat.data) && return DataFrame()
    (n_rows_show <= 0 || nrow(dat.data[end]) >= n_rows_show || length(dat.data) == 1) && return dat.data[end]
    acc = 0
    first_ep = 1
    for i = length(dat.data):-1:1
        acc += nrow(dat.data[i])
        if acc >= n_rows_show
            first_ep = i
            break
        end
    end
    return reduce(vcat, @view(dat.data[first_ep:end]))
end
