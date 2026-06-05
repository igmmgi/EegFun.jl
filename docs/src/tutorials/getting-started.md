# Getting Started with EegFun.jl

## Installing Julia

EegFun.jl requires [Julia](https://julialang.org/) 1.12.

The recommended way to install and manage Julia versions is with **[juliaup](https://github.com/JuliaLang/juliaup)**.
Alternatively, download an installer directly from the [Julia Downloads page](https://julialang.org/downloads/).

## The Julia REPL

Julia is an interactive language built around a Read-Eval-Print Loop (REPL). The REPL provides different modes accessed by special keys:

| Key | Mode | Purpose |
| --- | --- | --- |
| (default) | Julia mode | Execute Julia code |
| `]` | Package mode | Install and manage packages |
| `?` | Help mode | Access inline documentation |
| `;` | Shell mode | Run shell commands |

Press `Backspace` to return to Julia mode from any other mode.

## IDE Workflows

Most users pair the REPL with an editor that adds syntax highlighting and — most importantly — lets you send code directly into the live Julia REPL session.
See [**IDE Workflows**](ide-workflows.md) for a comparison of VS Code/VSCodium, Positron, JetBrains, and the Neovim + Iron.nvim terminal workflow.

## Installing EegFun

You can install EegFun.jl through the standard Julia package manager.

### Standard Installation

From the Julia REPL, enter Pkg mode by pressing `]` and run:

```julia
add EegFun
```

Or using `Pkg` in the code:

```julia
using Pkg
Pkg.add("EegFun")
```

### Development Version (vía GitHub)

To install the latest development version directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/igmmgi/EegFun.jl")
```

Then load the package in any Julia REPL session:

```julia
using EegFun
```

## First Steps

All EegFun functions are called with the `EegFun.` prefix. A minimal session could look something like this:

```julia
using EegFun

# Load a raw BioSemi recording
dat = EegFun.read_raw_data("participant1.bdf")

# Attach an electrode layout
layout = EegFun.read_layout("biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

# Browse the raw data interactively
EegFun.plot_databrowser(dat)
```

> [!TIP]
> Functions ending with `!` (e.g. `filter!`) mutate their input in-place. Functions without `!` return a new copy. Functions starting with `_` are internal helpers and not part of the public API.

## EegFun Philosophy

EegFun.jl is designed with ease-of-use as a core principle, making it accessible even for those without extensive programming experience. However, whilst EegFun provides many interactive GUIs for data visualization and exploration, it is not a full GUI application — the package emphasises a code-based workflow that is intended to be simple and intuitive.

The package offers a mix of high-level and lower-level functions, including complete analysis pipelines that take you from raw data through to ERP analyses, while still allowing fine-grained control when needed. In practice, a complete EEG analysis can be accomplished with little to zero traditional "coding" — simply typing commands in the Julia REPL and/or combining them into small, readable/runnable scripts.

## Supported Data Formats

| Format | Extension(s) | Notes |
| --- | --- | --- |
| BioSemi | `.bdf` | via [BiosemiDataFormat.jl](https://github.com/igmmgi/BiosemiDataFormat.jl) |
| BrainVision | `.vhdr` / `.eeg` / `.vmrk` | via [BrainVisionDataFormat.jl](https://github.com/igmmgi/BrainVisionDataFormat.jl) |
| EEGLAB | `.set` / `.fdt` | basic support |
| FieldTrip | `.mat` | basic support |

> [!NOTE]
> Additional file format support is planned for future releases. The format-specific packages are automatically installed as dependencies of EegFun.jl.

## Next Steps and Resources

| Resource | Link |
| --- | --- |
| EegFun.jl GitHub | [github.com/igmmgi/EegFun.jl](https://github.com/igmmgi/EegFun.jl) |
| Manual Preprocessing tutorial | [Manual Preprocessing](manual-preprocessing.md) |
| All how-to guides | [How-to Guides](../tutorials/workflows/preprocessing_workflow.md) |
| Julia learning resources | [julialang.org/learning](https://julialang.org/learning/) |
| Julia cheat sheet | [cheatsheet.juliadocs.org](https://cheatsheet.juliadocs.org/) |
| MATLAB–Python–Julia cheat sheet | [cheatsheets.quantecon.org](https://cheatsheets.quantecon.org/) |
| Makie.jl (plotting) | [docs.makie.org](https://docs.makie.org/) |
| DataFrames.jl | [dataframes.juliadata.org](https://dataframes.juliadata.org/) |
