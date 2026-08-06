import re
import os

files = [
    "src/analysis/erps/realign.jl",
    "src/analysis/processing/channel_average.jl",
    "src/analysis/processing/filter.jl",
    "src/analysis/processing/rereference.jl",
    "src/analysis/time_frequency/baseline.jl",
    "src/analysis/time_frequency/tf_morlet.jl",
    "src/analysis/time_frequency/tf_multitaper.jl",
    "src/analysis/time_frequency/tf_stft.jl",
    "src/utils/data.jl"
]

for file in files:
    with open(file, 'r') as f:
        content = f.read()
    
    # We will match from `files = _find_batch_files` to `_log_batch_summary`
    # Some functions have `# Transform output filenames` or `process_fn =` inside `_run_batch_operation`.
    # Wait, some might just do `process_fn = ...` before `files = _find_batch_files`.
    # No, all of them have `process_fn = ...` AFTER `files = _find_batch_files` EXCEPT maybe a few?
    # Let's just find `process_fn = ...` and `_run_batch_operation` anywhere in the try block!
    
    # We will search for:
    # 1. `files = _find_batch_files`
    # 2. `process_fn = ...`
    # 3. `_run_batch_operation`
    # 4. `_log_batch_summary`
    
    pattern = r"(        (?:# Find[^\n]*\n)?        files = _find_batch_files[\s\S]*?        process_fn = (.*?)\n[\s\S]*?_run_batch_operation\([^\n]*?operation_name\s*=\s*\"(.*?)\"[\s\S]*?_log_batch_summary[^\n]*?\n)"
    
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
    
    if count == 0:
        # maybe process_fn = ... is multiple lines
        pattern2 = r"(        (?:# Find[^\n]*\n)?        files = _find_batch_files[\s\S]*?        process_fn = ([\s\S]*?\n        end)\n[\s\S]*?_run_batch_operation\([^\n]*?operation_name\s*=\s*\"(.*?)\"[\s\S]*?_log_batch_summary[^\n]*?\n)"
        new_content, count = re.subn(pattern2, repl, content, flags=re.MULTILINE)
        
    if count == 0:
        # maybe no process_fn? wait, they all have process_fn.
        print(f"Failed to refactor {file}")
    else:
        with open(file, 'w') as f:
            f.write(new_content)
        print(f"Refactored {file}")

