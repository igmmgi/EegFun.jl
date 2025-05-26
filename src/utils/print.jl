
function print_vector_(v::Vector; max_length::Int = 10, n_ends::Int = 5)
    if length(v) > max_length
        v = vcat(first(v, n_ends), "...", last(v, n_ends))
    end
    return join(v, ", ")
end

function print_vector(v::UnitRange; max_length::Int = 10, n_ends::Int = 5)
    print_vector_(collect(v), max_length = max_length, n_ends = n_ends)
end

function print_vector(v::Vector; max_length::Int = 10, n_ends::Int = 5)
    print_vector_(collect(v), max_length = max_length, n_ends = n_ends)
end


"""
    print_config(config, [io=stdout]; [indent=0])
    print_config(config, filename::String; [indent=0])

Print configuration in a readable format.

# Arguments
- `config`: Configuration dictionary (typically loaded from TOML)
- `io=stdout`: Optional IO object to print to (default: standard output)
- `filename`: Optional filename to write output to
- `indent=0`: Indentation level (used internally for recursive calls)

# Example
```julia
# Print to console
print_config(config)

# Print to file
print_config(config, "config_dump.txt")
```
"""
function print_config(config, io::IO=stdout; indent::Int=0)
    # Sort just the keys, not the key-value pairs
    for key in sort(collect(keys(config)))
        value = config[key]
        if value isa Dict
            println(io, " "^indent, "$(key):")
            print_config(value, io; indent=indent+2)
        else
            # Format based on value type
            val_str = if value isa AbstractArray
                "[" * join(string.(value), ", ") * "]"
            else
                string(value)
            end
            println(io, " "^indent, "$(key): $(val_str)")
        end
    end
end

# Convenience method for printing to a file
function print_config(config, filename::String; indent::Int=0)
    open(filename, "w") do file
        print_config(config, file; indent=indent)
    end
    println("Configuration written to: $filename")
end

