"""
Test Runner and Coverage Analysis Tool for EyeFun

- Running tests with coverage
- Coverage analysis and reporting
- HTML report generation
- Cleanup of .cov files
- Interactive menu system

Usage:
    julia --project=. test/test_manager.jl [command] [options]

Commands:
    test                    - Run tests with coverage
    summary                 - Show coverage summary
    detailed                - Show detailed analysis
    file FILENAME           - Analyze specific file
    missed FILENAME         - Show missed code branches
    html                    - Generate HTML report
    clean                   - Remove all .cov files
    all                     - Run complete workflow (test + summary + html)
    interactive             - Show interactive menu
    help                    - Show this help message

If no command is provided, shows interactive menu.

Run from root directory with: julia --project=. test/test_manager.jl
"""

using Printf
using Pkg
# activate temp env and add packages needed for coverage
Pkg.activate(; temp=true)
Pkg.add(["Coverage", "CoverageTools"])
using Coverage
using CoverageTools

# Colours for output
const RED = :red
const GREEN = :green
const YELLOW = :yellow
const BLUE = :blue

print_colored(color::Symbol, message::String) = printstyled(message * "\n"; color=color)

function print_header()
    print_colored(BLUE, "=== EyeFun Test Runner and Coverage Analysis ===")
    println()
end

function run_tests_with_coverage()
    clean_coverage_files() # ensure fresh coverage data
    print_colored(YELLOW, "Step 1: Running tests with coverage...")
    try # Run tests with coverage=true
        run(`julia --project=. -e "using Pkg; Pkg.test(coverage=true)"`)
        print_colored(GREEN, "✓ Tests completed successfully")
    catch e
        print_colored(RED, "✗ Error running tests: $e")
    end
    println()
end

function get_coverage_stats(c)
    if isnothing(c.coverage)
        return (; covered=0, uncovered=0, not_executable=0, total_lines=0, percentage=0.0)
    end
    covered = count(x -> !isnothing(x) && x > 0, c.coverage)
    uncovered = count(x -> !isnothing(x) && x == 0, c.coverage)
    not_executable = count(x -> isnothing(x), c.coverage)
    total_lines = length(c.coverage)
    executable = covered + uncovered
    percentage = executable > 0 ? round(covered / executable * 100, digits=2) : 0.0
    return (; covered, uncovered, not_executable, total_lines, percentage)
end

function print_file_analysis(c, stats, max_uncovered=50)
    filename = replace(c.filename, "src/" => "")
    println("\n--- $filename ---")
    println("Total lines: $(stats.total_lines)")
    println("Covered lines: $(stats.covered)")
    println("Uncovered lines: $(stats.uncovered)")
    println("Not executable lines: $(stats.not_executable)")
    if stats.covered + stats.uncovered > 0
        println("Coverage percentage: $(stats.percentage)%")
    end

    uncovered_lines_list = findall(x -> !isnothing(x) && x == 0, c.coverage)
    if !isempty(uncovered_lines_list)
        if length(uncovered_lines_list) <= max_uncovered
            println("Uncovered lines: $(join(uncovered_lines_list, ", "))")
        else
            println(
                "Uncovered lines: $(join(uncovered_lines_list[1:max_uncovered], ", ")) ... (and $(length(uncovered_lines_list) - max_uncovered) more)",
            )
        end
    end
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
            stats = get_coverage_stats(c)
            if stats.covered + stats.uncovered > 0
                filename = replace(c.filename, "src/" => "")
                println("$filename: $(stats.percentage)% ($(stats.covered)/$(stats.covered + stats.uncovered) lines)")
                total_covered += stats.covered
                total_uncovered += stats.uncovered
            end
        end

        if total_covered + total_uncovered > 0
            overall_percentage =
                round(total_covered / (total_covered + total_uncovered) * 100, digits=2)
            println(
                "\nOverall Coverage: $overall_percentage% ($total_covered/$(total_covered + total_uncovered) lines)",
            )
        end

    catch e
        print_colored(RED, "Error: $e")
    end
    println()
end

function show_detailed_analysis()
    print_colored(YELLOW, "Step 3: Detailed Analysis")
    try
        coverage = process_folder("src")

        for c in coverage
            stats = get_coverage_stats(c)
            if stats.covered + stats.uncovered > 0
                print_file_analysis(c, stats, 20)
            end
        end

    catch e
        print_colored(RED, "Error: $e")
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
        stats = get_coverage_stats(c)
        print_file_analysis(c, stats, 50)

    catch e
        print_colored(RED, "Error: $e")
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
        if isfile(c.filename)
            lines = readlines(c.filename)

            println("File: $(c.filename)")
            println("Total lines: $(length(c.coverage))")

            # Show uncovered lines with context
            uncovered_count = 0
            for (i, cov) in enumerate(c.coverage)
                if !isnothing(cov) && cov == 0 && uncovered_count < 30
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
            open_cmd = Sys.isapple() ? "open" : (Sys.islinux() ? "xdg-open" : "start")
            println("Open with: $open_cmd test/coverage_html/index.html")
        catch
            print_colored(YELLOW, "⚠ genhtml not found. Install lcov")
            println("Then run: genhtml test/coverage.lcov -o test/coverage_html")
        end

    catch e
        print_colored(RED, "Error: $e")
    end
    println()
end

function clean_coverage_files()
    print_colored(YELLOW, "Cleaning up .cov files...")

    CoverageTools.clean_folder(".")
    print_colored(GREEN, "✓ Cleaned .cov files")

    # Also remove secondary artifacts
    lcov_file = "test/coverage.lcov"
    if isfile(lcov_file)
        println("Removing $lcov_file...")
        rm(lcov_file, force=true)
    end

    html_dir = "test/coverage_html"
    if isdir(html_dir)
        println("Removing $html_dir/ directory...")
        rm(html_dir, recursive=true, force=true)
    end

    println()
end

function run_all_analyses()
    print_colored(BLUE, "=== Running Complete Workflow ===")
    run_tests_with_coverage()
    show_coverage_summary()
    generate_html_report()
    print_colored(GREEN, "✓ Complete workflow finished")
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
            run_tests_with_coverage()
        elseif choice == "2"
            show_coverage_summary()
        elseif choice == "3"
            show_detailed_analysis()
        elseif choice == "4"
            print("Enter filename to analyze (e.g., eyelink_edf/ascii_exporter.jl): ")
            filename = readline()
            if !isempty(filename)
                analyze_specific_file(filename)
            end
        elseif choice == "5"
            print(
                "Enter filename to show missed branches (e.g., eyelink_edf/ascii_exporter.jl): ",
            )
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
    command = isempty(ARGS) ? "interactive" : ARGS[1]
    filename = length(ARGS) >= 2 ? ARGS[2] : ""

    if command == "test"
        run_tests_with_coverage()
    elseif command == "summary"
        show_coverage_summary()
    elseif command == "detailed"
        show_detailed_analysis()
    elseif command == "file"
        if isempty(filename)
            print_colored(RED, "Error: Please specify a filename")
            println(
                "Usage: julia --project=. test/test_manager.jl file eyelink_edf/ascii_exporter.jl",
            )
        else
            analyze_specific_file(filename)
        end
    elseif command == "missed"
        if isempty(filename)
            print_colored(RED, "Error: Please specify a filename")
            println(
                "Usage: julia --project=. test/test_manager.jl missed eyelink_edf/ascii_exporter.jl",
            )
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
    elseif command == "help" || command == "-h" || command == "--help"
        print_header()
        println("Usage: julia --project=. test/test_manager.jl [command] [options]")
        println()
        println("Commands:")
        println("  test                    - Run tests with coverage")
        println("  summary                 - Show coverage summary")
        println("  detailed                - Show detailed analysis")
        println("  file FILENAME           - Analyze specific file")
        println("  missed FILENAME         - Show missed code branches")
        println("  html                    - Generate HTML report")
        println("  clean                   - Remove all .cov files")
        println("  all                     - Run complete workflow (test + summary + html)")
        println("  interactive             - Show interactive menu")
        println("  help                    - Show this help message")
        println()
        println("Output files are saved in the test/ directory:")
        println("  - test/coverage.lcov     - LCOV coverage data")
        println("  - test/coverage_html/    - HTML coverage report")
        println()
        println("Examples:")
        println("  julia --project=. test/test_manager.jl")
        println("  julia --project=. test/test_manager.jl summary")
        println("  julia --project=. test/test_manager.jl file utils/data.jl")
        println("  julia --project=. test/test_manager.jl clean")
    else
        print_colored(RED, "Unknown command: $command")
        println("Use 'help' to see available commands.")
    end
end

main()
