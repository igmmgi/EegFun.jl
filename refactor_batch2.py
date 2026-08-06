import re
import os

files = [
    ("src/analysis/erps/realign.jl", "calculate_trigger_interval"),
    ("src/analysis/processing/baseline.jl", "baseline"),
    ("src/analysis/processing/channel_average.jl", "channel_average"),
    ("src/analysis/processing/condition_average.jl", "condition_average"),
    ("src/analysis/processing/condition_combine.jl", "condition_combine"),
    ("src/analysis/processing/condition_difference.jl", "condition_difference"),
    ("src/analysis/processing/epochs.jl", "average_epochs"),
    ("src/analysis/processing/filter.jl", "lowpass_filter"),
    ("src/analysis/processing/rereference.jl", "rereference"),
    ("src/analysis/processing/resample.jl", "resample"),
    ("src/analysis/time_frequency/baseline.jl", "tf_baseline"),
    ("src/analysis/time_frequency/tf_morlet.jl", "tf_morlet"),
    ("src/analysis/time_frequency/tf_multitaper.jl", "tf_multitaper"),
    ("src/analysis/time_frequency/tf_stft.jl", "tf_stft"),
    ("src/utils/data.jl", "subset")
]

for file, func_name in files:
    with open(file, 'r') as f:
        content = f.read()

    # The block we want to replace starts after `setup_global_logging` and ends at `finally\n        _cleanup_logging`
    
    # We will do a generic regex: 
    # find `log_file = ...\n    setup_global_logging(log_file)\n\n    try\n(.*?)\n        process_fn = (.*?)\n        results = _run_batch_operation\(.*?; operation_name = "(.*?)"\)\n\n        _log_batch_summary\(.*?\)\n\n    finally\n        _cleanup_logging\(.*?\)
    # 
    # Some functions like condition_difference have `result = (success = 0, errors = 0)` before `try`.

    pattern = r"(?:    result = \(success = 0, errors = 0\)\n)?    log_file = [^\n]*\n    setup_global_logging\([^\n]*\)\n\n(?:    result = [^\n]*\n)?    try\n(.*?)\n        process_fn = (.*?)\n        results = _run_batch_operation\([^\n]*; operation_name = \"(.*?)\"\)\n\n        (?:result = )?_log_batch_summary\([^\n]*\)\n\n    finally\n        _cleanup_logging\([^\n]*\)\n    end\n(?:    return result\n)?"
    
    def repl(m):
        body = m.group(1)
        process_fn = m.group(2)
        op_name = m.group(3)
        
        # We parse the body to extract validation vs logging vs directory setup vs file finding.
        
        # Extract validation (anything with _validate_ or @minimal_error or general checks)
        # Extract logging (@info, @log_call)
        
        lines = body.split('\n')
        validation_lines = []
        logging_lines = []
        default_dir_fn = "nothing"
        
        for line in lines:
            line_str = line.strip()
            if line_str.startswith("@info \"Batch") or line_str.startswith("@info \"Found") or "JLD2 files matching pattern" in line_str:
                continue
            if line_str.startswith("output_dir = something(output_dir,"):
                # extract default dir fn
                m2 = re.search(r"something\(output_dir, (.*?)\)", line_str)
                if m2:
                    call = m2.group(1)
                    # if it uses input_dir and file_pattern, we can make it a lambda.
                    # but it might use other arguments.
                    default_dir_fn = f"(dir, pat) -> {call.replace('input_dir', 'dir').replace('file_pattern', 'pat')}"
                continue
            if line_str.startswith("mkpath(output_dir)") or line_str.startswith("files = _find_batch_files") or line_str.startswith("if isempty(files)") or line_str.startswith("@minimal_warning \"No JLD2 files") or line_str.startswith("return nothing") or line_str.startswith("return result") or line_str == "end":
                continue
                
            if line_str.startswith("@info") or line_str.startswith("@log_call"):
                logging_lines.append(line)
            elif line_str:
                validation_lines.append(line)

        # Assemble new block
        new_block = ""
        if validation_lines:
            new_block += "\n".join(validation_lines) + "\n\n"
            
        new_block += f"    process_fn = {process_fn}\n\n"
        
        new_block += f"    return batch_process(process_fn, file_pattern;\n"
        new_block += f"        input_dir, participant_selection, output_dir,\n"
        
        if default_dir_fn != "nothing":
            new_block += f"        default_output_dir_fn = {default_dir_fn},\n"
            
        if logging_lines:
            new_block += f"        init_log_fn = () -> begin\n"
            new_block += "\n".join(logging_lines) + "\n"
            new_block += f"        end,\n"
            
        new_block += f"        operation_name = \"{op_name}\")\n"
        
        return new_block
        
    new_content, count = re.subn(pattern, repl, content, flags=re.DOTALL)
    
    # Let's write it back
    if count > 0:
        with open(file, 'w') as f:
            f.write(new_content)
        print(f"Refactored {file}")
    else:
        print(f"Failed to refactor {file}")
