import re
import os

files = [
    ("src/analysis/erps/realign.jl", "calculate_trigger_interval"),
    ("src/analysis/erps/lrp.jl", "lrp"),
    ("src/analysis/processing/baseline.jl", "baseline"),
    ("src/analysis/processing/channel_average.jl", "channel_average"),
    ("src/analysis/processing/condition_average.jl", "condition_average"),
    ("src/analysis/processing/condition_combine.jl", "condition_combine"),
    ("src/analysis/processing/condition_difference.jl", "condition_difference"),
    ("src/analysis/processing/epochs.jl", "average_epochs"),
    ("src/analysis/processing/filter.jl", "_run_filter_batch"),
    ("src/analysis/processing/rereference.jl", "rereference"),
    ("src/analysis/processing/resample.jl", "resample"),
    ("src/analysis/time_frequency/baseline.jl", "tf_baseline"),
    ("src/analysis/time_frequency/tf_morlet.jl", "tf_morlet"),
    ("src/analysis/time_frequency/tf_multitaper.jl", "tf_multitaper"),
    ("src/analysis/time_frequency/tf_stft.jl", "tf_stft"),
    ("src/utils/data.jl", "subset")
]

for file, func in files:
    with open(file, 'r') as f:
        content = f.read()

    # Search for files = _find_batch_files
    # Search for process_fn = ...
    # Search for _run_batch_operation ... operation_name = "..."
    # Search for _log_batch_summary
    
    # We will match from `files = _find_batch_files` to `_log_batch_summary`
    
    pattern = r"(        # Find(?: JLD2)? files\n        files = _find_batch_files[\s\S]*?        process_fn = ([^\n]*?)\n[\s\S]*?_run_batch_operation\([^\n]*?operation_name\s*=\s*\"(.*?)\"[\s\S]*?_log_batch_summary[^\n]*?\n)"
    
    def repl(m):
        process_fn_def = m.group(2)
        op_name = m.group(3)
        
        has_parallel = "parallel = parallel" in m.group(0)
        
        ret = f"        process_fn = {process_fn_def}\n"
        ret += f"        result = batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, \"{op_name}\""
        if has_parallel:
            ret += ", parallel = parallel"
        ret += ")\n"
        
        return ret

    new_content, count = re.subn(pattern, repl, content, flags=re.MULTILINE)
    
    if count > 0:
        with open(file, 'w') as f:
            f.write(new_content)
        print(f"Refactored {file}")
    else:
        print(f"Failed to refactor {file}")

