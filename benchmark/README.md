# EegFun.jl Benchmark Suite

This directory contains the cross-platform benchmarking suite for comparing EEG processing pipelines in Julia (EegFun.jl), Python (MNE-Python), and MATLAB (EEGLAB).

## Running the Benchmarks

All benchmark scripts expect the following arguments:
1. `data_dir` - The path to the directory containing your raw EEG data files (in `.bdf` format).
2. `n_files_to_process` - (Optional) The number of files to process. For a quick test, you can set this to `2`.
3. `run_ica` - (Optional) A boolean flag to determine if ICA should be run (defaults to `true` or equivalent in each script).

Replace `/path/to/your/data/directory` in the examples below with the actual folder path containing your `.bdf` files.

### 1. Julia (EegFun.jl)

Run from your terminal:

```bash
cd julia
julia eegfun_benchmark.jl /path/to/your/data/directory 2
```

### 2. Python (MNE-Python)

Run from your terminal:

```bash
cd python
python mne_benchmark.py /path/to/your/data/directory 2
```

### 3. MATLAB (EEGLAB)

You can run the MATLAB script directly from your terminal using the `-batch` flag:

```bash
cd matlab
matlab -batch "eeglab_benchmark('/path/to/your/data/directory', 2)"
```

*(Alternatively, you can run it directly inside the MATLAB command window):*

```matlab
cd matlab
eeglab_benchmark('/path/to/your/data/directory', 2)
```

## Generating the Comparison Plots

Once all three benchmarks have finished successfully, they will output their timing results and ERP data as `.csv` and `.txt` files into your data directory. You can then generate the comparison plots using the included visualization script:

```bash
# Run from the root of the benchmark directory
julia plot_benchmark_comparison.jl /path/to/your/data/directory
```

This will generate two PDF files in your data directory:
- `cross_pipeline_comparison.pdf` - A combined plot overlaying the ERPs from all three platforms.
- `individual_pipelines.pdf` - Individual subplots for each platform.
