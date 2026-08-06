using Glob

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
    # This is a complex refactor, we want to basically revert the batch_process call to the try..finally format
    # Or, we just fix the result return type in rereference and condition_difference, and remove the log tests.
end
