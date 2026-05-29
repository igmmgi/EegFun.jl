# Benchmark

`EegFun.jl` includes a benchmark suite to compare its execution time against other standard EEG toolboxes like **MNE-Python** and **EEGLAB** (MATLAB).

> [!NOTE]
> The goal of this pipeline is **not** to demonstrate the ideal EEG/ERP processing methodology (it is certainly not!). Rather, it provides a standardized, basic workflow—read data, re-reference, filter, run ICA, extract epochs, and average—to measure realistic performance characteristics across different software platforms.

## Downloading the Dataset

The benchmark uses the same raw Visual Attention (Posner Cueing) dataset featured in our other tutorials. You can download the complete dataset from Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19045958.svg)](https://doi.org/10.5281/zenodo.19045958)

Once downloaded, you can extract the files into any folder on your computer (e.g., `~/Downloads/AttentionExp`). The benchmark scripts will accept this directory as a command line argument and save all output files directly into that same folder.

## Prerequisites & Installation

To run the full suite, you need to ensure the respective toolboxes are set up on your machine.

**Julia (EegFun.jl)**
All dependencies (Makie, BenchmarkTools, DelimitedFiles) are automatically resolved when you run the scripts using the `--project=.` flag from the repository root.

**Python (MNE-Python)**
We recommend using the lightning-fast [`uv`](https://github.com/astral-sh/uv) package manager to create an isolated Python virtual environment for the benchmark:

```bash
# 1. Create a virtual environment
uv venv

# 2. Activate it (Linux/macOS)
source .venv/bin/activate
# (On Windows: .venv\Scripts\activate)

# 3. Install MNE-Python and Matplotlib
uv pip install mne matplotlib
```

**MATLAB (EEGLAB)**
Ensure you have downloaded EEGLAB and added it to your MATLAB path. The benchmark script relies on standard EEGLAB core functions and requires the `biosig` plugin to load the raw `.bdf` files.

## Running the Benchmark Suite

The official benchmark scripts are located in the repository under the `benchmark/` folder.

This script demonstrates a completely realistic, full-scale pipeline from raw data to a Grand Average plot:

1. Loads raw BDF data and custom electrode layouts.
2. Re-references to the average and applies high-pass filtering (automatically distributed across CPU threads).
3. Identifies and flags extreme artifacts.
4. Executes Extended Infomax Independent Component Analysis (ICA).
    - *Note on ICA:* In this automated script, we simply remove the single largest variance component (IC1) without manual review. While IC1 is highly likely to reflect an eye blink, no care is taken to verify this. This pipeline is strictly designed as a **timing benchmark** for computational throughput, and to provide a rough comparison of the resulting visual attention effect across parietal-occipital electrodes.
5. Extracts specific epoch conditions and applies baseline correction.
6. Averages the trials into Event-Related Potentials (ERPs).
7. Computes a Grand Average across multiple subjects.
8. Generates final comparison figures.

### Running the Benchmarks

All benchmark scripts (Julia, Python, and MATLAB) require you to pass the absolute path to the downloaded dataset. They also accept optional arguments to control the execution length and whether to run Independent Component Analysis (ICA).

**Arguments:**

1. `data_dir` *(Required)*: Absolute path to the folder containing the downloaded Zenodo `.bdf` files.
2. `n_files` *(Optional, Default: `0`)*: Number of files to process. `0` means process all files in the directory.
3. `run_ica` *(Optional, Default: `true`)*: Boolean flag to enable or disable the ICA step. Disabling ICA is highly recommended for isolating raw I/O and mathematical filtering throughput.

**Julia Execution:**
We highly recommend launching Julia with multi-threading enabled (`--threads=auto`) to accurately measure the pipeline's concurrent performance characteristics.

```bash
# Example 1: Process 3 files and disable ICA
julia --threads=auto --project=. benchmark/julia/eegfun_benchmark.jl /path/to/downloaded/AttentionExp 3 false

# Example 2: Process ALL files (0) and disable ICA
julia --threads=auto --project=. benchmark/julia/eegfun_benchmark.jl /path/to/downloaded/AttentionExp 0 false
```

**Python (MNE) Execution:**
```bash
python benchmark/python/mne_benchmark.py /path/to/downloaded/AttentionExp 3 false
```

**MATLAB (EEGLAB) Execution:**
Because MATLAB handles command-line arguments poorly, you must run the benchmark function directly from your MATLAB command window. Make sure you first add EEGLAB and the benchmark script folder to your MATLAB path!

```matlab
% 1. Add EEGLAB to path and initialize it
addpath('/path/to/eeglab');
eeglab;

% 2. Add the benchmark script folder to path
addpath('/path/to/EegFun.jl/benchmark/matlab');

% 3. Run the benchmark function
eeglab_benchmark('/path/to/downloaded/AttentionExp', 3, false)
```

### Plotting the Results

Once you have run the benchmarks for the different toolboxes, they will automatically save their execution times (`.txt`) and exported ERP arrays (`.csv`) directly into your dataset directory.

You can then run the Julia comparison script to generate a single, unified plot that overlays the resulting ERP waveforms and automatically embeds the total execution times directly into the legend!

```bash
# Generate the cross-pipeline comparison plot
julia --project=. benchmark/plot_benchmark_comparison.jl /path/to/downloaded/AttentionExp
```
