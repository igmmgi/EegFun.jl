# === DATA VIEWING FUNCTIONS ===

"""
    viewer(dat)

Display data using VS Code viewer if available, otherwise use standard display.
If `all_data(dat)` is defined for the type, it is called first to extract a
DataFrame — so adding viewer support for a new type only requires adding an
`all_data` method in data.jl.

# Arguments
- `dat`: Any data object to display
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
    head(dat::EegData; n=5)

Display and return the first `n` rows of the EEG data.

# Examples
```julia
head(erps[1])        # First 5 rows
head(erps[1], n=10)  # First 10 rows
```
"""
function head(dat::EegData; n = nothing)
    isnothing(n) && (n = 5)
    data = all_data(dat)
    nrows = nrow(data)
    n = min(n, nrows)
    result = n > 0 ? data[1:n, :] : DataFrame()
    viewer(result)
    return result
end

"""
    tail(dat::EegData; n=5)

Display and return the last `n` rows of the EEG data.

# Examples
```julia
tail(erps[1])        # Last 5 rows
tail(erps[1], n=10)  # Last 10 rows
```
"""
function tail(dat::EegData; n = nothing)
    isnothing(n) && (n = 5)
    data = all_data(dat)
    nrows = nrow(data)
    n = min(n, nrows)
    result = n > 0 ? data[max(1, nrows-n+1):nrows, :] : DataFrame()
    viewer(result)
    return result
end
