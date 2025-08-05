# Simple script to fix the TOML file
using Base.Filesystem

# Read the file
content = read("src/config/default.toml", String)

# Replace function references with string values
content = replace(content, "filter_func = filtfilt" => "filter_func = \"filtfilt\"")
content = replace(content, "# Type: Function" => "# Type: String")

# Write back
write("src/config/default.toml", content)

println("Fixed TOML file!")
