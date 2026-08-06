using Base.Meta

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
    "src/utils/data.jl"
]

for file in files
    content = read(file, String)
    
    # We want to extract @log_call "func_name" and any @info ... that are before process_fn
    # We will do a regex replacement.
    # Pattern to match:
    # Any number of @log_call or @info lines at the start of the function body.
    
    # Let's find the function signature end
    m = match(r"function \w+\(.*?\)\n((?:    @log_call.*?\n|    @info.*?\n|    \n)*)"s, content)
    
    # Actually, a simpler way is to find `@log_call "something"` and capture the @info around it.
    
    # Alternatively, just use sed or manually replace them.
end
