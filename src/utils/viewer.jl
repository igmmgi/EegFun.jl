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
    data = all_data(dat)
    nrows = nrow(data)
    k = min(n_rows_show, nrows)
    result = k > 0 ? view(data, 1:k, :) : DataFrame()
    viewer(result)
    return result
end

"""
    tail(dat::EegData; n = 5)

Display and return a view of the last `n` rows of `all_data(dat)`.
"""
function tail(dat::EegData; n = nothing)
    n_rows_show = isnothing(n) ? 5 : n
    data = all_data(dat)
    nrows = nrow(data)
    k = min(n_rows_show, nrows)
    result = k > 0 ? view(data, max(1, nrows-k+1):nrows, :) : DataFrame()
    viewer(result)
    return result
end
