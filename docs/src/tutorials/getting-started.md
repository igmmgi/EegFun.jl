# Getting Started with EegFun.jl

## Requirements

- [Julia](https://julialang.org/) 1.10 (LTS) or later

### Installing Julia

If you need to install Julia, visit the [Julia Downloads page](https://julialang.org/downloads/).

### Using the Julia REPL

Julia is an interactive language that can be used from a Read-Eval-Print Loop (REPL). The REPL provides different modes accessed by special keys:

- **Julia mode** (default) - Execute Julia code
- **Package mode** - Press `]` to manage packages
- **Help mode** - Press `?` to access documentation
- **Shell mode** - Press `;` to run shell commands

Press `Backspace` to return to Julia mode from any other mode.

### Installing Julia Packages

Julia packages are installed using the built-in package manager. In the Julia REPL, use `Pkg.add("PackageName")` to install packages, or enter package mode by pressing `]` and typing `add PackageName`.

> [!NOTE]
> Unregistered packages can be installed directly from GitHub using the repository URL: `Pkg.add(url="https://github.com/username/package.jl")`

### Installation of EegFun

> [!IMPORTANT]
> EegFun.jl is not currently registered in the Julia General registry and must be installed directly from GitHub.

```julia
using Pkg
Pkg.add(url="https://github.com/igmmgi/EegFun.jl")
```

Load the package:

```julia
using EegFun
```

Using the package:

```julia
output = EegFun.XXX # Replace XXX with the actual function name
```

#### Julia Function Naming Conventions

Julia uses specific naming conventions used within EegFun that are important to understand:

- Functions ending with `!` (e.g., `filter!`) mutate their input arguments
- Functions without `!` (e.g., `filter`) return new data without modifying inputs
- Functions starting with `_` (e.g., `_internal_helper`) are considered internal and not part of the public API

## EegFun Philosophy

EegFun.jl is designed with ease-of-use as a core principle, making it accessible even for those without extensive programming experience. While EegFun provides many interactive GUIs for data visualization and exploration, it is not a full GUI application. Instead, the package emphasizes a code-based workflow that remains simple and intuitive.

The package offers a mix of high-level and lower-level functions, providing complete analysis pipelines that take you from raw data through to ERP analyses, while still allowing fine-grained control when needed. In practice, a complete EEG analysis pipeline can be accomplished with little to zero traditional "coding" — simply typing commands in the Julia REPL and/or combining them into small, readable scripts. This approach provides the flexibility and reproducibility of a programming environment while maintaining the accessibility of interactive tools needed for data visualization and exploration.

### Useful Julia Resources

- [Julia Learning Resources](https://julialang.org/learning/)

#### Julia Cheat Sheets

- [Julia Cheat Sheets](https://cheatsheet.juliadocs.org/)
- [Matlab-Python-Julia Cheat Sheet](https://cheatsheets.quantecon.org/)

### Key Julia Packages

All interactive plotting functionality in EegFun.jl is provided by the excellent [Makie.jl](https://docs.makie.org/) package for high-performance data visualization. EegFun uses:

- [GLMakie](https://docs.makie.org/stable/explanations/backends/glmakie/) - GPU-accelerated backend for interactive plots
- [CairoMakie](https://docs.makie.org/stable/explanations/backends/cairomakie/) - Vector graphics backend for publication-quality figures

Other packages used within EegFun:

- [DataFrames.jl](https://dataframes.juliadata.org/) - Tabular data manipulation and analysis
- [DSP.jl](https://docs.juliadsp.org/) - Digital signal processing (filtering, spectral analysis)
- [JLD2.jl](https://juliaio.github.io/JLD2.jl/) - Efficient Julia data serialization

## Supported Data Formats

EegFun.jl currently supports reading EEG data from the following formats:

- **BioSemi** (.bdf) - Supported via [BiosemiDataFormat.jl](https://github.com/igmmgi/BiosemiDataFormat.jl)
- **BrainVision** (.vhdr/.eeg/.vmrk) - Supported via [BrainVisionDataFormat.jl](https://github.com/igmmgi/BrainVisionDataFormat.jl)

Basic support is also available for EEGLAB (.set/.fdt) and FieldTrip (.mat) file formats.

> [!NOTE]
> Additional file format support is planned for future releases.

The format-specific packages (BiosemiDataFormat.jl and BrainVisionDataFormat.jl) are automatically installed as dependencies when you install EegFun.jl.
