#!/bin/bash

# Check if data directory is provided
if [ -z "$1" ]; then
    echo "Error: Data directory not provided."
    echo "Usage: ./run_all_benchmarks.sh /path/to/data/directory [n_files] [run_ica]"
    echo "Example: ./run_all_benchmarks.sh /path/to/data 2 true"
    exit 1
fi

# Convert to absolute path
DATA_DIR=$(realpath "$1")
N_FILES=${2:-2}
RUN_ICA=${3:-true}

# Get the absolute path of the benchmark directory
BENCHMARK_DIR=$(dirname "$(realpath "$0")")

echo "=================================================="
echo "Starting EegFun.jl Cross-Platform Benchmark Suite"
echo "Data Directory: $DATA_DIR"
echo "Files to Process: $N_FILES"
echo "Run ICA: $RUN_ICA"
echo "=================================================="

# 1. Julia Benchmark
echo ""
echo ">>> Running Julia (EegFun.jl) Benchmark..."
cd "$BENCHMARK_DIR/julia" || exit
julia -t auto --project="$BENCHMARK_DIR/.." eegfun_benchmark.jl "$DATA_DIR" "$N_FILES" "$RUN_ICA"

# 2. Python Benchmark
echo ""
echo ">>> Running Python (MNE-Python) Benchmark..."
cd "$BENCHMARK_DIR/python" || exit
python mne_benchmark.py "$DATA_DIR" "$N_FILES" "$RUN_ICA"

# 3. MATLAB Benchmark
echo ""
echo ">>> Running MATLAB (EEGLAB) Benchmark..."
cd "$BENCHMARK_DIR/matlab" || exit
matlab -batch "eeglab_benchmark('$DATA_DIR', $N_FILES, '$RUN_ICA')"

# 4. Generate Comparison Plots
echo ""
echo ">>> Generating Comparison Plots..."
cd "$BENCHMARK_DIR" || exit
julia --project="$BENCHMARK_DIR/.." plot_benchmark_comparison.jl "$DATA_DIR"

echo ""
echo "=================================================="
echo "All benchmarks completed!"
echo "Data and plots have been generated in:"
echo "$DATA_DIR/benchmarks"
echo "=================================================="
