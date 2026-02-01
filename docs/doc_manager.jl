# Now load the packages
using Documenter
using DocumenterTools
using DocumenterVitepress
using JuliaFormatter
using LiveServer
using Logging
using Printf

# Add the parent directory to the load path so we can load the local package
push!(LOAD_PATH, dirname(@__DIR__))
using EegFun

"""
Documentation Manager for EegFun.jl

This is a comprehensive documentation management tool that provides:
- Building documentation with Documenter.jl
- Documentation coverage analysis
- Code formatting with JuliaFormatter
- Cleanup and maintenance tasks
- Interactive menu system

Usage:
    julia --project=. docs/build_docs.jl [command] [options]

Commands:
    build                    - Build documentation
    coverage                 - Check documentation coverage
    clean                    - Clean build artifacts
    interactive              - Show interactive menu
    all                      - Run complete documentation workflow

If no command is provided, shows interactive menu.
"""


# Colors for output
const RED = "\033[0;31m"
const GREEN = "\033[0;32m"
const YELLOW = "\033[1;33m"
const BLUE = "\033[0;34m"
const CYAN = "\033[0;36m"
const NC = "\033[0m" # No Color

function print_colored(color::String, message::String)
    println("$color$message$NC")
end

function print_header()
    print_colored(BLUE, "=== EegFun.jl Documentation Manager ===")
    println()
end

function get_project_root()
    # The script is in docs/, so project root is one level up
    return dirname(@__DIR__)
end


function build_documentation(project_root::String)

    try
        # Build documentation
        build_jl_path = joinpath(project_root, "docs", "make.jl")
        print_colored(GREEN, " Building documentation with Documenter.jl...")

        # Suppress warnings during build
        old_logger = global_logger()
        try
            # Use a logger that only shows errors
            logger = ConsoleLogger(stderr, Logging.Error)
            global_logger(logger)
            include(build_jl_path)
        finally
            # Restore original logger
            global_logger(old_logger)
        end
        print_colored(GREEN, " Documentation built successfully")

    catch e
        print_colored(RED, " Error building documentation: $e")
        return false
    end

    println()
    return true
end




function check_doc_coverage(project_root::String; skip_build_check::Bool = false)
    print_colored(YELLOW, "Checking documentation coverage...")

    try
        # Check for basic documentation files
        doc_files = [
            joinpath(project_root, "docs", "src", "index.md"),
            joinpath(project_root, "docs", "src", "api.md"),
            joinpath(project_root, "docs", "make.jl"),
        ]
        missing_files = []

        for file in doc_files
            if !isfile(file)
                push!(missing_files, file)
            end
        end

        if !isempty(missing_files)
            print_colored(RED, " Missing documentation files:")
            for file in missing_files
                println("  - $file")
            end
            return false
        end

        # Check if documentation has been built (only if not skipping)
        if !skip_build_check
            build_dir = joinpath(project_root, "docs", "build")
            if !isdir(build_dir)
                print_colored(YELLOW, " Documentation not built yet. Run 'build' first.")
                return true
            end
        end

        # Analyze source code documentation (docstrings)
        println("\nSource Code Documentation Analysis:")
        println("="^50)

        # Find all Julia source files
        src_dir = joinpath(project_root, "src")
        source_files = String[]
        for (root, dirs, files) in walkdir(src_dir)
            for file in files
                if endswith(file, ".jl")
                    push!(source_files, joinpath(root, file))
                end
            end
        end

        total_functions = 0
        documented_functions = 0
        total_docstring_chars = 0
        files_with_docs = 0

        println(" Analyzing $(length(source_files)) source files...")

        for file in source_files
            try
                content = read(file, String)
                lines = split(content, '\n')

                # Count functions and docstrings
                file_functions = 0
                file_documented = 0
                file_doc_chars = 0

                i = 1
                while i <= length(lines)
                    line = strip(lines[i])

                    # Look for function definitions
                    if occursin(r"^function\s+\w+", line) || occursin(r"^\w+\(.*\)\s*=", line)
                        file_functions += 1
                        total_functions += 1

                        # Check if there's a docstring above (look back up to 3 lines)
                        docstring_found = false
                        for j = max(1, i - 3):i-1
                            if j <= length(lines)
                                prev_line = strip(lines[j])
                                if startswith(prev_line, "\"\"\"") || (startswith(prev_line, "\"") && !endswith(prev_line, "\""))
                                    file_documented += 1
                                    documented_functions += 1
                                    docstring_found = true

                                    # Count docstring characters
                                    if startswith(prev_line, "\"\"\"")
                                        # Multi-line docstring
                                        doc_start = j
                                        doc_end = j
                                        for k = j+1:length(lines)
                                            if endswith(strip(lines[k]), "\"\"\"")
                                                doc_end = k
                                                break
                                            end
                                        end
                                        for k = doc_start:doc_end
                                            file_doc_chars += length(strip(lines[k]))
                                        end
                                    else
                                        # Single line docstring
                                        file_doc_chars += length(prev_line)
                                    end
                                    break
                                end
                            end
                        end
                    end
                    i += 1
                end

                if file_documented > 0
                    files_with_docs += 1
                    total_docstring_chars += file_doc_chars
                    println("  $(basename(file)): $file_documented/$file_functions functions documented")
                end

            catch e
                println("  Error reading $file: $e")
            end
        end

        # Calculate coverage percentage
        coverage_percent = total_functions > 0 ? round((documented_functions / total_functions) * 100, digits = 1) : 0

        println("\n Documentation Coverage Summary:")
        println("   Total functions: $total_functions")
        println("   Documented functions: $documented_functions")
        println("   Coverage: $coverage_percent%")
        println("   Files with documentation: $files_with_docs/$(length(source_files))")
        println("   Total docstring characters: $total_docstring_chars")

        print_colored(GREEN, "\n Documentation coverage analysis completed")

    catch e
        print_colored(RED, " Error checking documentation coverage: $e")
        return false
    end

    println()
    return true
end

function clean_docs(project_root::String)
    print_colored(YELLOW, "Cleaning documentation build artifacts...")

    # Clean build directory
    build_dir = joinpath(project_root, "docs", "build")
    if isdir(build_dir)
        rm(build_dir, recursive = true)
        print_colored(GREEN, "✓ Removed docs/build directory")
    else
        print_colored(YELLOW, "No build directory found")
    end

    # Clean other common build artifacts
    artifacts = ["site", ".documenter", "Manifest.toml"]
    for artifact in artifacts
        artifact_path = joinpath(project_root, "docs", artifact)
        if isdir(artifact_path) || isfile(artifact_path)
            rm(artifact_path, recursive = true)
            print_colored(GREEN, "✓ Removed docs/$artifact")
        end
    end

    print_colored(GREEN, "✓ Documentation cleanup completed")
    println()
end



function format_source_files(project_root::String)
    print_colored(YELLOW, "Formatting Julia source files...")

    try
        # Find all Julia files in src directory
        src_dir = joinpath(project_root, "src")
        julia_files = String[]
        for (root, dirs, files) in walkdir(src_dir)
            for file in files
                if endswith(file, ".jl")
                    push!(julia_files, joinpath(root, file))
                end
            end
        end

        if isempty(julia_files)
            print_colored(YELLOW, "No Julia files found in src directory")
            return true
        end

        println("Found $(length(julia_files)) Julia file(s) to format:")
        for file in julia_files
            println("  - $file")
        end

        # Format files
        formatted_count = 0
        for file in julia_files
            try
                JuliaFormatter.format_file(file)
                formatted_count += 1
            catch e
                print_colored(RED, " Error formatting $file: $e")
            end
        end

        print_colored(GREEN, " Successfully formatted $formatted_count/$(length(julia_files)) file(s)")

        # Also format test files
        test_dir = joinpath(project_root, "test")
        test_files = String[]
        for (root, dirs, files) in walkdir(test_dir)
            for file in files
                if endswith(file, ".jl") && !endswith(file, ".cov")
                    push!(test_files, joinpath(root, file))
                end
            end
        end

        if !isempty(test_files)
            println("\nFormatting $(length(test_files)) test file(s)...")
            for file in test_files
                try
                    JuliaFormatter.format_file(file)
                catch e
                    print_colored(RED, " Error formatting $file: $e")
                end
            end
        end

    catch e
        print_colored(RED, " Error during formatting: $e")
        return false
    end

    println()
    return true
end

function format_and_check(project_root::String)
    print_colored(YELLOW, "Formatting and checking Julia files...")

    if !format_source_files(project_root)
        return false
    end

    # Run a quick syntax check
    print_colored(YELLOW, "Running syntax check...")
    try
        Pkg.precompile()
        print_colored(GREEN, " Syntax check passed")
    catch e
        print_colored(RED, " Syntax check failed: $e")
        return false
    end

    print_colored(GREEN, " Formatting and syntax check completed")
    println()
    return true
end

function view_documentation(project_root::String)
    print_colored(YELLOW, "Starting VitePress dev server...")

    build_dir = joinpath(project_root, "docs", "build")
    if !isdir(build_dir)
        print_colored(RED, " Documentation not built yet. Run 'Build documentation' first.")
        return false
    end

    try
        LiveServer.serve(dir = "docs/build/1")
        return true
    catch e
        if isa(e, InterruptException)
            print_colored(YELLOW, "\n Server stopped")
        else
            print_colored(RED, " Error starting server: $e")
            return false
        end
    end

    println()
    return true
end




function run_all_docs(project_root::String)
    print_colored(GREEN, "Running complete documentation workflow...")
    println()

    # Step 1: Format source files
    print_colored(YELLOW, "Step 1: Formatting source files...")
    if !format_source_files(project_root)
        print_colored(RED, " Formatting failed")
        return false
    end

    # Step 2: Build documentation
    print_colored(YELLOW, "Step 2: Building documentation...")
    if !build_documentation(project_root)
        print_colored(RED, " Documentation build failed")
        return false
    end

    # Step 3: Check documentation coverage
    print_colored(YELLOW, "Step 3: Checking documentation coverage...")
    check_doc_coverage(project_root)

    print_colored(GREEN, "=== Documentation Workflow Complete ===")
    println("Next steps:")
    println("1. Review the built documentation in docs/build/")
    println("2. Open docs/build/index.html in your browser to view documentation")
    println("3. Deploy to GitHub Pages when ready")

    return true
end

function show_interactive_menu(project_root::String)
    print_header()

    while true
        println("\nChoose an option:")
        println("1. Build documentation")
        println("2. Check documentation coverage")
        println("3. Format source files")
        println("4. Format and check syntax")
        println("5. Clean build artifacts")
        println("6. Run complete workflow")
        println("7. View documentation (via LiveServer)")
        println("8. Exit")

        print("\nEnter your choice (1-8): ")
        choice = readline()

        if choice == "1"
            build_documentation(project_root)
            if isfile(joinpath(project_root, "docs", "build", "index.html"))
                print_colored(GREEN, " Documentation built successfully!")
                print_colored(CYAN, "Open docs/build/index.html in your browser to view the documentation")
            end
        elseif choice == "2"
            check_doc_coverage(project_root)
        elseif choice == "3"
            format_source_files(project_root)
        elseif choice == "4"
            format_and_check(project_root)
        elseif choice == "5"
            clean_docs(project_root)
        elseif choice == "6"
            run_all_docs(project_root)
        elseif choice == "7"
            view_documentation(project_root)
        elseif choice == "8"
            break
        else
            print_colored(RED, "Invalid choice. Please enter 1-8.")
        end
    end
end

function main()
    project_root = get_project_root()

    # Simple command line argument parsing
    if length(ARGS) == 0
        command = "interactive"
    elseif length(ARGS) == 1
        command = ARGS[1]
    else
        command = ARGS[1]
        # Additional arguments could be handled here
    end

    if command == "build"
        build_documentation(project_root)
    elseif command == "coverage"
        check_doc_coverage(project_root)
    elseif command == "format"
        format_source_files(project_root)
    elseif command == "format-check"
        format_and_check(project_root)
    elseif command == "clean"
        clean_docs(project_root)
    elseif command == "all"
        run_all_docs(project_root)
    elseif command == "interactive"
        show_interactive_menu(project_root)
    else
        print_header()
        println("Usage: julia --project=. docs/build_docs.jl [command]")
        println()
        println("Commands:")
        println("  build                    - Build documentation")
        println("")
        println("")
        println("  coverage                 - Check documentation coverage")
        println("  format                   - Format Julia source files")
        println("  format-check             - Format files and run syntax check")
        println("  clean                    - Clean build artifacts")
        println("  interactive              - Show interactive menu")
        println("  all                      - Run complete documentation workflow")
        println()
        println("Examples:")
        println("  julia --project=. docs/build_docs.jl build")
        println("  julia --project=. docs/build_docs.jl format")
        println("  julia --project=. docs/build_docs.jl clean")
        println("  julia --project=. docs/build_docs.jl interactive")
    end
end

# Run the main function
main()
