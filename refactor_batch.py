import re
import os

files = [
    "src/analysis/erps/lrp.jl",
    "src/analysis/erps/realign.jl",
    "src/analysis/processing/baseline.jl",
    "src/analysis/processing/channel_average.jl",
    "src/analysis/processing/condition_average.jl",
    "src/analysis/processing/condition_combine.jl",
    "src/analysis/processing/condition_difference.jl",
    "src/analysis/processing/epochs.jl",
    "src/analysis/processing/filter.jl",
    "src/analysis/processing/rereference.jl",
    "src/analysis/processing/resample.jl",
    "src/analysis/time_frequency/baseline.jl",
    "src/analysis/time_frequency/tf_morlet.jl",
    "src/analysis/time_frequency/tf_multitaper.jl",
    "src/analysis/time_frequency/tf_stft.jl",
    "src/utils/data.jl" # subset
]

for file in files:
    if "lrp.jl" in file: continue # already did this one manually
    with open(file, 'r') as f:
        content = f.read()

    # Find the pattern
    # log_file = "..."
    # setup_global_logging(log_file)
    # 
    # [optional result = ...]
    # try
    #    ...
    #    @info "Found ... files"
    #    [optional extra @infos]
    #    process_fn = ...
    #    results = _run_batch_operation(..., operation_name="...")
    #    _log_batch_summary(...)
    # finally
    #    _cleanup_logging(...)
    # end
    
    # We can do this file by file manually if there are only 15 left, or script it if it's too tedious.
    # Actually, a simpler regex just looks for `process_fn = ...` and `batch_process(process_fn, ...)`
    # Oh wait! The files in my workspace currently HAVE the try...finally block.
    # To refactor them, I need to do what I did manually.
