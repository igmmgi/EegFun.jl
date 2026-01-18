#!/usr/bin/env julia

"""
Coverage Analysis Script for eegfun

This script provides various coverage analysis functions that can be run
individually or as a complete workflow.

Usage:
    julia --project=. analyze_coverage.jl [command]

Commands:
    summary     - Show coverage summary for all files
    detailed    - Show detailed analysis for all files
    file FILE   - Analyze specific file (e.g., "file utils/data.jl")
    missed FILE - Show missed code branches for specific file
    html        - Generate HTML coverage report
    all         - Run all analyses

If no command is provided, runs 'summary' by default.
"""

using Coverage
using CoverageTools

function run_tests_with_coverage()
    println("=== Step 1: Running tests with coverage ===")
    println("Running: julia --project=. -e \"using Pkg; Pkg.test(coverage=true)\"")
    println("(Please run this command manually first to generate coverage data)")
    println()
end

function show_coverage_summary()
    println("=== Step 2: Coverage Summary ===")

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
        println("Error: $e")
        println("Make sure you've run tests with coverage first!")
    end
    println()
end

function show_detailed_analysis()
    println("=== Step 3: Detailed Analysis ===")

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
        println("Error: $e")
        println("Make sure you've run tests with coverage first!")
    end
    println()
end

function analyze_specific_file(target_file::String)
    println("=== Analyzing: $target_file ===")

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
        println("Error: $e")
        println("Make sure you've run tests with coverage first!")
    end
    println()
end

function show_missed_branches(target_file::String)
    println("=== Missed Code Branches: $target_file ===")

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
        println("Error: $e")
        println("Make sure you've run tests with coverage first!")
    end
    println()
end

function generate_html_report()
    println("=== Step 4: Generating HTML Coverage Report ===")

    try
        coverage = process_folder("src")

        # Generate LCOV file in test directory
        lcov_file = "test/coverage.lcov"
        CoverageTools.LCOV.writefile(lcov_file, coverage)
        println("✓ LCOV file generated: $lcov_file")

        # Check if genhtml is available
        try
            run(`which genhtml`)
            println("✓ genhtml found, generating HTML report...")
            run(`genhtml test/coverage.lcov -o test/coverage_html`)
            println("✓ HTML report generated: test/coverage_html/index.html")
            println("Open with: open test/coverage_html/index.html")
        catch
            println("⚠ genhtml not found. Install with: brew install lcov")
            println("Then run: genhtml test/coverage.lcov -o test/coverage_html")
        end

    catch e
        println("Error: $e")
        println("Make sure you've run tests with coverage first!")
    end
    println()
end

function run_all_analyses()
    println("=== Complete Coverage Analysis Workflow ===")
    println()

    run_tests_with_coverage()
    show_coverage_summary()
    show_detailed_analysis()
    generate_html_report()

    println("=== Analysis Complete ===")
    println("Next steps:")
    println("1. Run tests with coverage: julia --project=. -e \"using Pkg; Pkg.test(coverage=true)\"")
    println("2. Re-run this script to see results")
    println("3. Use 'file FILENAME' to analyze specific files")
    println("4. Use 'missed FILENAME' to see missed code branches")
end

# Main execution
function main()
    if length(ARGS) == 0
        command = "summary"
    else
        command = ARGS[1]
    end

    if command == "summary"
        show_coverage_summary()
    elseif command == "detailed"
        show_detailed_analysis()
    elseif command == "file" && length(ARGS) >= 2
        analyze_specific_file(ARGS[2])
    elseif command == "missed" && length(ARGS) >= 2
        show_missed_branches(ARGS[2])
    elseif command == "html"
        generate_html_report()
    elseif command == "all"
        run_all_analyses()
    else
        println("Usage: julia --project=. analyze_coverage.jl [command]")
        println()
        println("Commands:")
        println("  summary     - Show coverage summary for all files")
        println("  detailed    - Show detailed analysis for all files")
        println("  file FILE   - Analyze specific file (e.g., \"file utils/data.jl\")")
        println("  missed FILE - Show missed code branches for specific file")
        println("  html        - Generate HTML coverage report")
        println("  all         - Run all analyses")
        println()
        println("Examples:")
        println("  julia --project=. analyze_coverage.jl")
        println("  julia --project=. analyze_coverage.jl file utils/data.jl")
        println("  julia --project=. analyze_coverage.jl missed utils/data.jl")
    end
end

# Run if called directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
