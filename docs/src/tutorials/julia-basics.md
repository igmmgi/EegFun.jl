# Julia Basics

EegFun.jl is designed so that you do not need to be a Julia programmer to use
it. In practice, most EEG analysis with EegFun involves calling functions and
interacting with plots — there is very little traditional "programming"
required, and in that sense the choice of language is almost secondary.
That said, a basic understanding of Julia always helps, especially when you
want to write small scripts, customise a pipeline, or understand an error
message. This page covers the essentials.

## The REPL as a Calculator

The simplest way to start with Julia is to type expressions directly into the
REPL. Julia evaluates them immediately and prints the result:

```julia-repl
julia> 2 + 3
5

julia> 7 * 8.5
59.5

julia> 2^10
1024

julia> sqrt(144)
12.0
```

Standard mathematical operators work as expected: `+`, `-`, `*`, `/`, `^`
(exponentiation), and `%` (remainder). Parentheses control evaluation order:

```julia-repl
julia> (3 + 4) * 2
14
```

## Variables

Assign values to names with `=`. Julia is dynamically typed, so you do not
need to declare a type:

```julia-repl
julia> sample_rate = 512
512

julia> duration = 2.5
2.5

julia> n_samples = sample_rate * duration
1280.0
```

Variable names can include Unicode characters, which is useful for writing
code that reads like textbook notation:

```julia-repl
julia> μ = 0.0
0.0

julia> σ = 1.5
1.5

julia> Δt = 1 / 512
0.001953125
```

> [!TIP]
> In the REPL or VS Code, type the LaTeX name and press `Tab` to insert
> Unicode symbols: `\mu` → `μ`, `\sigma` → `σ`, `\Delta` → `Δ`.

To see all variables defined in your current session (similar to MATLAB's
`whos`), use `varinfo()`:

```julia-repl
julia> varinfo()
  name                    size summary
  –––––––––––––––– ––––––––––– –––––––
  Δt                    8 bytes Float64
  μ                     8 bytes Float64
  σ                     8 bytes Float64
  duration              8 bytes Float64
  n_samples             8 bytes Float64
  sample_rate           8 bytes Int64
```

## Types

Every value in Julia has a type. You can inspect it with `typeof`:

```julia-repl
julia> typeof(42)
Int64

julia> typeof(3.14)
Float64

julia> typeof("hello")
String

julia> typeof(true)
Bool
```

## Vectors and Arrays

Square brackets create vectors (1D arrays):

```julia-repl
julia> voltages = [1.2, -0.5, 3.1, 0.8]
4-element Vector{Float64}:
 1.2
 -0.5
  3.1
  0.8

julia> voltages[1]       # Julia uses 1-based indexing
1.2

julia> voltages[end]     # 'end' refers to the last element
0.8

julia> voltages[2:3]     # slicing
2-element Vector{Float64}:
 -0.5
  3.1
```

Create ranges and regular sequences:

```julia-repl
julia> 1:5               # a range from 1 to 5
1:5

julia> collect(1:5)      # materialise into a vector
5-element Vector{Int64}:
 1
 2
 3
 4
 5

julia> 0:0.5:2           # start:step:stop
0.0:0.5:2.0
```

Matrices (2D arrays, channels × time points in EEG) use semicolons or spaces:

```julia-repl
julia> data = [1.0 2.0 3.0; 4.0 5.0 6.0]
2×3 Matrix{Float64}:
 1.0  2.0  3.0
 4.0  5.0  6.0

julia> size(data)
(2, 3)

julia> data[1, :]        # first row (e.g. first channel)
3-element Vector{Float64}:
 1.0
 2.0
 3.0
```

## Control Flow

### if-else

```julia
if p_value < 0.05
    println("Significant")
else
    println("Not significant")
end
```

### for loops

```julia
channels = ["Fz", "Cz", "Pz"]
for ch in channels
    println("Processing channel: $ch")
end
```

Loops in Julia are fast — unlike MATLAB or Python, there is no performance
penalty for writing explicit loops instead of vectorised code.

### Comprehensions

Concise syntax for building arrays from loops:

```julia-repl
julia> [i^2 for i in 1:5]
5-element Vector{Int64}:
  1
  4
  9
 16
 25
```

## Strings

Strings use double quotes. String interpolation uses `$`:

```julia-repl
julia> participant = "P01"
"P01"

julia> condition = "congruent"
"congruent"

julia> filename = "data/$(participant)_$(condition).bdf"
"data/P01_congruent.bdf"
```

Single characters use single quotes (`'a'`). Use `string()` or `*` to
concatenate strings:

```julia-repl
julia> "hello" * " " * "world"
"hello world"
```

## Dictionaries

Key–value storage, useful for mapping trigger codes to condition labels:

```julia-repl
julia> triggers = Dict(1 => "standard", 2 => "deviant", 3 => "novel")
Dict{Int64, String} with 3 entries:
  2 => "deviant"
  3 => "novel"
  1 => "standard"

julia> triggers[2]
"deviant"
```

## Tuples

Tuples are fixed-size, immutable collections. You have already seen them —
`size()` returns a tuple:

```julia-repl
julia> dims = (64, 512)
(64, 512)

julia> dims[1]
64

julia> typeof(dims)
Tuple{Int64, Int64}
```

Unlike vectors, tuples cannot be modified after creation. They are useful for
grouping a small number of related values, such as a time window:

```julia-repl
julia> baseline_window = (-0.2, 0.0)
(-0.2, 0.0)
```

## Named Tuples

Named tuples add field names to each element, giving you lightweight,
self-documenting containers without defining a custom type:

```julia-repl
julia> participant = (id = "P01", age = 25, group = "control")
(id = "P01", age = 25, group = "control")

julia> participant.age
25

julia> participant[:group]
"control"
```

Named tuples are immutable like regular tuples. They are commonly used for
passing groups of options or returning multiple values from a function.

## Functions

### Standard Form

Functions are defined with the `function ... end` block:

```julia-repl
julia> function add(a, b)
           return a + b
       end
add (generic function with 1 method)

julia> add(3, 5)
8
```

### Short Form

For simple one-liners, the same function can be written more concisely:

```julia-repl
julia> add(a, b) = a + b
add (generic function with 1 method)

julia> add(3, 5)
8
```

### Anonymous Functions

Short throwaway functions use the `->` syntax. These are common with `map` and
`filter`:

```julia-repl
julia> map(x -> x^2, [1, 2, 3])
3-element Vector{Int64}:
 1
 4
 9
```

### Mutating Functions

By convention, functions that modify their input end with `!`:

```julia
lowpass_filter!(data, 30.0)   # modifies data in-place
lowpass_filter(data, 30.0)    # returns a new copy, data is unchanged
```

This is a naming convention, not enforced by the language, but Julia packages
(including EegFun.jl) follow it consistently.

### Broadcasting (The Dot Syntax)

Adding a dot (`.`) before an operator or after a function name applies it
element-wise to each value in an array:

```julia-repl
julia> a = [1, 2, 3, 4]

julia> sin.(a)            # apply sin to each element
4-element Vector{Float64}:
 0.8414709848078965
 0.9092974268256817
 0.1411200080598672
 -0.7568024953079282

```

## Using Packages

Julia packages are loaded with `using`. For example, the `Statistics` standard
library provides `mean`:

```julia-repl
julia> using Statistics

julia> data = [1.2, -0.5, 3.1, 0.8, 2.4]

julia> mean(data)
1.4

julia> std(data)
1.3266499161421599
```

The first `using` in a session triggers compilation (see
[Why Julia? — Trade-offs](why-julia.md#trade-offs)). Subsequent calls in the
same session are instant.

## Multiple Dispatch

In Julia, the same function name can have different **methods** depending on
the types of its arguments. Julia automatically picks the right method:

```julia-repl
julia> describe(x::Int) = println("$x is an integer")
julia> describe(x::Float64) = println("$x is a floating-point number")
julia> describe(x::String) = println("$x is a string")

julia> describe(42)
42 is an integer

julia> describe(3.14)
3.14 is a floating-point number

julia> describe("hello")
hello is a string
```

You can check how many methods a function has with `methods`:

```julia-repl
julia> methods(describe)
# 3 methods for generic function "describe"
```

Types in Julia are organised in a hierarchy that you can explore with
`supertype` and `subtypes`. Both `Int64` and `Float64` are subtypes of
`Number`, so you can write a single method that accepts any number:

```julia-repl
julia> double(x::Number) = x * 2

julia> double(3)
6

julia> double(3.14)
6.28
```

This is how EegFun.jl works internally — many functions have different methods
for different data types (continuous, epoched, averaged), so you use the same
function name regardless of what you pass in. See
[**Data Structures**](../explanations/data-structures.md) for EegFun.jl's own
type hierarchy.

## Getting Help

The REPL help mode (press `?`) gives you inline documentation:

```julia-repl
help?> sqrt
search: sqrt isqrt

  sqrt(x)

  Return √x.
```

For EegFun functions:

```julia-repl
help?> EegFun.lowpass_filter!
```

## Next Steps

| Resource | Link |
| --- | --- |
| Julia Manual | [docs.julialang.org](https://docs.julialang.org/en/v1/) |
| Julia learning resources | [julialang.org/learning](https://julialang.org/learning/) |
| MATLAB–Python–Julia cheat sheet | [cheatsheets.quantecon.org](https://cheatsheets.quantecon.org/) |
| Think Julia (free book) | [benlauwens.github.io/ThinkJulia.jl](https://benlauwens.github.io/ThinkJulia.jl/latest/book.html) |
