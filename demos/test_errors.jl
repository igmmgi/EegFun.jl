"""
Test script for the new error handling system.

Run interactively in the REPL to compare error output:
    julia> include("demos/test_errors.jl")
"""

using EegFun

println("="^60)
println("Testing EegFun Error Handling")
println("="^60)

# ──────────────────────────────────────────────────────────
# 1. EegFunError — custom exception with clean output
# ──────────────────────────────────────────────────────────
println("\n--- 1. Direct EegFunError throw ---")
try
    throw(EegFun.EegFunError("This is a clean user-facing error message"))
catch e
    println("Caught: ", e)
    println("Type:   ", typeof(e))
    println("Sprint: ", sprint(showerror, e))
end

# ──────────────────────────────────────────────────────────
# 2. @minimal_error macro — logs + throws EegFunError
# ──────────────────────────────────────────────────────────
println("\n--- 2. @minimal_error macro ---")
try
    EegFun.@minimal_error "Something went wrong with your data"
catch e
    println("Caught: ", typeof(e))
    println("Message: ", e.msg)
end

# ──────────────────────────────────────────────────────────
# 3. Compare: Julia default error (with stacktrace spew)
# ──────────────────────────────────────────────────────────
println("\n--- 3. Standard Julia error (for comparison) ---")
try
    error("Standard Julia error with full stacktrace")
catch e
    println("Caught: ", typeof(e))
    println("Message: ", e.msg)
end

# ──────────────────────────────────────────────────────────
# 4. @minimal_warning — stays the same
# ──────────────────────────────────────────────────────────
println("\n--- 4. @minimal_warning ---")
EegFun.@minimal_warning "This is a clean warning"

# ──────────────────────────────────────────────────────────
# 5. Real-world: trigger an actual EegFun error
# ──────────────────────────────────────────────────────────
println("\n--- 5. Real-world: plot_topography with empty data ---")
try
    EegFun.plot_topography(EegFun.TimeFreqData[], freq_range = (4.0, 12.0))
catch e
    if e isa EegFun.EegFunError
        println("✓ Got clean EegFunError: ", e.msg)
    else
        println("✗ Got $(typeof(e)) instead of EegFunError")
        println("  Message: ", sprint(showerror, e))
    end
end

# ──────────────────────────────────────────────────────────
# 6. Uncaught error — see what the REPL shows
# ──────────────────────────────────────────────────────────
println("\n--- 6. Uncomment below to see uncaught EegFunError in REPL ---")
println("# throw(EegFun.EegFunError(\"Clean uncaught error\"))")
println("# EegFun.@minimal_error \"Clean uncaught error via macro\"")

println("\n✓ All error tests completed!")
