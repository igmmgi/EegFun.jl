using Coverage

# Find all .cov files
src_dir = joinpath(@__DIR__, "..", "src")
test_dir = @__DIR__

println("Looking for .cov files in:")
println("  - src: $src_dir")
println("  - test: $test_dir")

# Get all .cov files from src directory and its subdirectories
src_cov_files = String[]
for (root, dirs, files) in walkdir(src_dir)
    for file in files
        if endswith(file, ".cov")
            push!(src_cov_files, joinpath(root, file))
        end
    end
end

println("\nFound .cov files in src:")
for file in src_cov_files
    println("  - $file")
end

# Extract the original source files from .cov files, preserving subdirectory structure
src_files = unique(map(f -> begin
    rel_path = relpath(f, src_dir)
    replace(rel_path, r"\.\d+\.cov$" => "")
end, src_cov_files))

println("\nFound coverage data for source files:")
for file in src_files
    println("  - $file")
end

# Process coverage information
println("\nCoverage by file:")

# Process source files
if !isempty(src_files)
    coverage = Coverage.CoverageTools.process_file.([joinpath(src_dir, file) for file in src_files])
    for file in coverage
        covered, total = get_summary(file)
        if total > 0
            println("$(relpath(file.filename, src_dir)): $(covered/total * 100)%")
        end
    end
end
