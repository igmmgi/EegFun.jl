#!/usr/bin/env julia

"""
Test Runner and Coverage Analysis Tool for eegfun

This is a pure Julia equivalent of test.sh that provides:
- Running tests with coverage
- Coverage analysis and reporting
- HTML report generation
- Cleanup of .cov files
- Interactive menu system

Usage:
    julia --project=. test/run_tests.jl [command] [options]

Commands:
    test                    - Run tests with coverage
    summary                 - Show coverage summary
    detailed                - Show detailed analysis
    file FILENAME           - Analyze specific file
    missed FILENAME         - Show missed code branches
    html                    - Generate HTML report
    clean                   - Remove all .cov files
    interactive             - Show interactive menu
    all                     - Run complete workflow

If no command is provided, shows interactive menu.
"""

using Printf

# Load coverage packages from extras
using Pkg
Pkg.activate(; temp = true)
Pkg.add(["Coverage", "CoverageTools"])
using Coverage
using CoverageTools


# Colors for output
const RED = "\033[0;31m"
const GREEN = "\033[0;32m"
const YELLOW = "\033[1;33m"
const BLUE = "\033[0;34m"
const NC = "\033[0m" # No Color

function print_colored(color::String, message::String)
    println("$color$message$NC")
end

function print_header()
    print_colored(BLUE, "=== eegfun Test Runner and Coverage Analysis ===")
    println()
end



function show_coverage_summary()
    print_colored(YELLOW, "Step 2: Coverage Summary")

    try
        coverage = process_folder("src")

        println("Coverage Summary:")
        println("=================")

        total_covered = 0
        total_uncovered = 0

        for c in coverage
            if c.coverage !== nothing
                covered = count(x -> x !== nothing && x > 0, c.coverage)
                uncovered = count(x -> x !== nothing && x == 0, c.coverage)
                total = covered + uncovered

                if total > 0
                    percentage = round(covered / total * 100, digits = 2)
                    filename = replace(c.filename, "src/" => "")
                    println("$filename: $percentage% ($covered/$total lines)")

                    total_covered += covered
                    total_uncovered += uncovered
                end
            end
        end

        if total_covered + total_uncovered > 0
            overall_percentage = round(total_covered / (total_covered + total_uncovered) * 100, digits = 2)
            println(
                "\nOverall Coverage: $overall_percentage% ($total_covered/$(total_covered + total_uncovered) lines)",
            )
        end

    catch e
        print_colored(RED, "Error: $e")
        println("Make sure you've run tests with coverage first!")
    end
    println()
end

function show_detailed_analysis()
    print_colored(YELLOW, "Step 3: Detailed Analysis")

    try
        coverage = process_folder("src")

        for c in coverage
            if c.coverage !== nothing
                covered_lines = count(x -> x !== nothing && x > 0, c.coverage)
                uncovered_lines = count(x -> x !== nothing && x == 0, c.coverage)
                not_executable = count(x -> x === nothing, c.coverage)
                total_lines = length(c.coverage)

                if covered_lines + uncovered_lines > 0
                    percentage = round(covered_lines / (covered_lines + uncovered_lines) * 100, digits = 2)
                    filename = replace(c.filename, "src/" => "")

                    println("\n--- $filename ---")
                    println("Total lines: $total_lines")
                    println("Covered lines: $covered_lines")
                    println("Uncovered lines: $uncovered_lines")
                    println("Not executable lines: $not_executable")
                    println("Coverage percentage: $percentage%")

                    # Show uncovered line numbers (first 20)
                    uncovered_lines_list = Int[]
                    for (i, cov) in enumerate(c.coverage)
                        if cov !== nothing && cov == 0
                            push!(uncovered_lines_list, i)
                        end
                    end

                    if !isempty(uncovered_lines_list)
                        if length(uncovered_lines_list) <= 20
                            println("Uncovered lines: $(join(uncovered_lines_list, ", "))")
                        else
                            println(
                                "Uncovered lines: $(join(uncovered_lines_list[1:20], ", ")) ... (and $(length(uncovered_lines_list) - 20) more)",
                            )
                        end
                    end
                end
            end
        end

    catch e
        print_colored(RED, "Error: $e")
        println("Make sure you've run tests with coverage first!")
    end
    println()
end

function analyze_specific_file(target_file::String)
    print_colored(YELLOW, "Analyzing: $target_file")

    try
        coverage = process_folder("src")
        data_coverage = filter(c -> occursin(target_file, c.filename), coverage)

        if isempty(data_coverage)
            println("No coverage data found for $target_file")
            return
        end

        c = data_coverage[1]
        covered_lines = count(x -> x !== nothing && x > 0, c.coverage)
        uncovered_lines = count(x -> x !== nothing && x == 0, c.coverage)
        not_executable = count(x -> x === nothing, c.coverage)
        total_lines = length(c.coverage)

        println("File: $(c.filename)")
        println("Total lines: $total_lines")
        println("Covered lines: $covered_lines")
        println("Uncovered lines: $uncovered_lines")
        println("Not executable lines: $not_executable")

        if covered_lines + uncovered_lines > 0
            percentage = round(covered_lines / (covered_lines + uncovered_lines) * 100, digits = 2)
            println("Coverage percentage: $percentage%")
        end

        # Show uncovered line numbers
        uncovered_lines_list = Int[]
        for (i, cov) in enumerate(c.coverage)
            if cov !== nothing && cov == 0
                push!(uncovered_lines_list, i)
            end
        end

        if !isempty(uncovered_lines_list)
            println("\nUncovered line numbers:")
            if length(uncovered_lines_list) <= 50
                println(join(uncovered_lines_list, ", "))
            else
                println(join(uncovered_lines_list[1:50], ", "), " ... (and $(length(uncovered_lines_list) - 50) more)")
            end
        end

    catch e
        print_colored(RED, "Error: $e")
        println("Make sure you've run tests with coverage first!")
    end
    println()
end

function show_missed_branches(target_file::String)
    print_colored(YELLOW, "Missed Code Branches: $target_file")

    try
        coverage = process_folder("src")
        data_coverage = filter(c -> occursin(target_file, c.filename), coverage)

        if isempty(data_coverage)
            println("No coverage data found for $target_file")
            return
        end

        c = data_coverage[1]

        # Read source file
        if isfile(c.filename)
            lines = readlines(c.filename)

            println("File: $(c.filename)")
            println("Total lines: $(length(c.coverage))")

            # Show uncovered lines with context
            uncovered_count = 0
            for (i, cov) in enumerate(c.coverage)
                if cov !== nothing && cov == 0 && uncovered_count < 30
                    uncovered_count += 1
                    line_num = i
                    if line_num <= length(lines)
                        start_line = max(1, line_num - 2)
                        end_line = min(length(lines), line_num + 2)

                        println("\n--- Around line $line_num ---")
                        for j = start_line:end_line
                            marker = j == line_num ? ">>> " : "    "
                            println("$marker$j: $(lines[j])")
                        end
                    end
                end
            end

            if uncovered_count >= 30
                println("\n... (showing first 30 uncovered lines)")
            end
        else
            println("Source file not found: $(c.filename)")
        end

    catch e
        print_colored(RED, "Error: $e")
        println("Make sure you've run tests with coverage first!")
    end
    println()
end

function generate_html_report()
    print_colored(YELLOW, "Step 4: Generating HTML Coverage Report")

    try
        coverage = process_folder("src")

        # Generate LCOV file in test directory
        lcov_file = "test/coverage.lcov"
        CoverageTools.LCOV.writefile(lcov_file, coverage)
        print_colored(GREEN, "✓ LCOV file generated: $lcov_file")

        # Check if genhtml is available
        try
            run(`which genhtml`)
            print_colored(GREEN, "✓ genhtml found, generating HTML report...")
            run(`genhtml test/coverage.lcov -o test/coverage_html`)
            print_colored(GREEN, "✓ HTML report generated: test/coverage_html/index.html")
            println("Open with: open test/coverage_html/index.html")
        catch
            print_colored(YELLOW, "⚠ genhtml not found. Install with: brew install lcov")
            println("Then run: genhtml test/coverage.lcov -o test/coverage_html")
        end

    catch e
        print_colored(RED, "Error: $e")
        println("Make sure you've run tests with coverage first!")
    end
    println()
end

function clean_coverage_files()
    print_colored(YELLOW, "Cleaning up .cov files...")

    # Find .cov files in current directory and test subdirectory
    cov_files = String[]

    # Search current directory
    for (root, dirs, files) in walkdir(".")
        for file in files
            if endswith(file, ".cov")
                push!(cov_files, joinpath(root, file))
            end
        end
    end

    # Search test directory if it exists
    if isdir("test")
        for (root, dirs, files) in walkdir("test")
            for file in files
                if endswith(file, ".cov")
                    push!(cov_files, joinpath(root, file))
                end
            end
        end
    end

    if isempty(cov_files)
        print_colored(GREEN, "No .cov files found to clean")
        return
    end

    println("Found $(length(cov_files)) .cov file(s) to remove:")
    for file in cov_files
        println("  - $file")
    end

    # Remove files
    for file in cov_files
        rm(file, force = true)
    end

    print_colored(GREEN, "✓ Successfully removed $(length(cov_files)) .cov file(s)")
    println()
end






function show_interactive_menu()
    print_header()

    while true
        println("\nChoose an option:")
        println("1. Run tests with coverage")
        println("2. Show coverage summary")
        println("3. Show detailed analysis")
        println("4. Analyze specific file")
        println("5. Show missed code branches")
        println("6. Generate HTML report")
        println("7. Clean .cov files")
        println("8. Exit")

        print("\nEnter your choice (1-8): ")
        choice = readline()

        if choice == "1"
            print_colored(YELLOW, "Running tests with coverage...")
            try
                run(`julia --project=. -e "using Pkg; Pkg.test(coverage=true)"`)
                print_colored(GREEN, "✓ Tests completed successfully")
            catch e
                print_colored(RED, "✗ Error running tests: $e")
            end
            println()
        elseif choice == "2"
            show_coverage_summary()
        elseif choice == "3"
            show_detailed_analysis()
        elseif choice == "4"
            print("Enter filename to analyze (e.g., utils/data.jl): ")
            filename = readline()
            if !isempty(filename)
                analyze_specific_file(filename)
            end
        elseif choice == "5"
            print("Enter filename to show missed branches (e.g., utils/data.jl): ")
            filename = readline()
            if !isempty(filename)
                show_missed_branches(filename)
            end
        elseif choice == "6"
            generate_html_report()
        elseif choice == "7"
            clean_coverage_files()
        elseif choice == "8"
            break
        else
            print_colored(RED, "Invalid choice. Please enter 1-8.")
        end
    end
end


function main()
    # Simple command line argument parsing
    if length(ARGS) == 0
        command = "interactive"
        filename = ""
    elseif length(ARGS) == 1
        command = ARGS[1]
        filename = ""
    else
        command = ARGS[1]
        filename = ARGS[2]
    end

    if command == "test"
        run_tests_with_coverage()
    elseif command == "summary"
        show_coverage_summary()
    elseif command == "detailed"
        show_detailed_analysis()
    elseif command == "file"
        if isempty(filename)
            print_colored(RED, "Error: Please specify a filename")
            println("Usage: julia --project=. test/run_tests.jl file utils/data.jl")
        else
            analyze_specific_file(filename)
        end
    elseif command == "missed"
        if isempty(filename)
            print_colored(RED, "Error: Please specify a filename")
            println("Usage: julia --project=. test/run_tests.jl missed utils/data.jl")
        else
            show_missed_branches(filename)
        end
    elseif command == "html"
        generate_html_report()
    elseif command == "clean"
        clean_coverage_files()
    elseif command == "all"
        run_all_analyses()
    elseif command == "interactive"
        show_interactive_menu()
    else
        print_header()
        println("Usage: julia --project=. test/run_tests.jl [command] [options]")
        println()
        println("Commands:")
        println("  test                    - Run tests with coverage")
        println("  summary                 - Show coverage summary")
        println("  detailed                - Show detailed analysis")
        println("  file FILENAME           - Analyze specific file")
        println("  missed FILENAME         - Show missed code branches")
        println("  html                    - Generate HTML report")
        println("  clean                   - Remove all .cov files")
        println("  interactive             - Show interactive menu")
        println("  all                     - Run complete workflow")
        println()
        println("Output files are saved in the test/ directory:")
        println("  - test/coverage.lcov     - LCOV coverage data")
        println("  - test/coverage_html/    - HTML coverage report")
        println()
        println("Examples:")
        println("  julia --project=. test/test_manager.jl all")
        println("  julia --project=. test/test_manager.jl file utils/data.jl")
        println("  julia --project=. test/test_manager.jl missed utils/data.jl")
        println("  julia --project=. test/test_manager.jl clean")
    end
end

# Run the main function
main()
