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
        
    # We want to replace the whole block starting from `files = _find_batch_files` 
    # up to `_log_batch_summary` or `result = _log_batch_summary`
    # and replace it with process_fn definition and batch_process.
    
    # regex to find it:
    pattern = r"(        # Find files\n        files = _find_batch_files[\s\S]*?        process_fn = ([^\n]*?)\n\n        (?:results = )?_run_batch_operation\([^\n]*?; operation_name = \"(.*?)\"(?:, parallel = parallel)?\)\n\n        (?:result = )?_log_batch_summary\([^\n]*?\)\n)"
    
    def repl(m):
        process_fn_def = m.group(2)
        op_name = m.group(3)
        
        # some use parallel, some don't. We can pass it if the original had it.
        # But for now let's just pass `parallel = parallel` if the function signature had it.
        # It's safer to just look at the original call.
        has_parallel = ", parallel = parallel" in m.group(0)
        
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

