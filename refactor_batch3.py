import re
import os

files = [
    "src/analysis/erps/realign.jl",
    "src/analysis/processing/channel_average.jl",
    "src/analysis/processing/condition_average.jl",
    "src/analysis/processing/condition_difference.jl",
    "src/analysis/processing/epochs.jl",
    "src/analysis/processing/filter.jl",
    "src/analysis/processing/rereference.jl",
    "src/analysis/time_frequency/baseline.jl",
    "src/analysis/time_frequency/tf_morlet.jl",
    "src/analysis/time_frequency/tf_multitaper.jl",
    "src/analysis/time_frequency/tf_stft.jl",
]

for file in files:
    with open(file, 'r') as f:
        lines = f.readlines()
        
    new_lines = []
    i = 0
    in_try = False
    
    # We will search for setup_global_logging.
    # When we find it, we skip until `try`. Then we collect lines.
    while i < len(lines):
        line = lines[i]
        
        if "setup_global_logging" in line:
            # We found the logging setup.
            # Skip lines until try
            while not lines[i].strip() == "try":
                i += 1
            i += 1 # skip 'try'
            
            # Now we are inside the try block.
            # We collect all lines until `finally`
            body_lines = []
            while not lines[i].strip() == "finally":
                body_lines.append(lines[i])
                i += 1
                
            # skip finally and its block
            while not lines[i].strip() == "end":
                i += 1
            i += 1 # skip 'end'
            
            # if the next line is `return result`, skip it
            if i < len(lines) and "return result" in lines[i]:
                i += 1
                
            # Now we parse the body_lines
            validations = []
            logs = []
            process_fn_line = ""
            default_dir_fn = "nothing"
            op_name = ""
            
            for b in body_lines:
                b_str = b.strip()
                if "Batch" in b_str and "started at" in b_str: continue
                if b_str.startswith("@info \"Found") or "JLD2 files matching pattern" in b_str: continue
                if b_str.startswith("output_dir = something(output_dir,"):
                    m = re.search(r"something\(output_dir, (.*?)\)", b_str)
                    if m:
                        call = m.group(1).replace('input_dir', 'dir').replace('file_pattern', 'pat')
                        default_dir_fn = f"(dir, pat) -> {call}"
                    continue
                if b_str.startswith("process_fn ="):
                    process_fn_line = b
                    continue
                if b_str.startswith("mkpath(") or b_str.startswith("files = ") or b_str.startswith("if isempty") or b_str.startswith("@minimal_warning") or b_str.startswith("return nothing") or b_str.startswith("return result") or b_str == "end":
                    continue
                if b_str.startswith("results = _run_batch_operation"):
                    m = re.search(r'operation_name\s*=\s*"([^"]+)"', b_str)
                    if m: op_name = m.group(1)
                    continue
                if "_log_batch_summary" in b_str: continue
                if b_str.startswith("result = (success"): continue
                
                if b_str.startswith("@info") or b_str.startswith("@log_call"):
                    logs.append(b)
                elif b_str:
                    validations.append(b)
                    
            # emit new lines
            new_lines.extend(validations)
            if validations: new_lines.append("\n")
            
            if process_fn_line:
                new_lines.append(process_fn_line)
                new_lines.append("\n")
            
            new_lines.append("    return batch_process(process_fn, file_pattern;\n")
            new_lines.append("        input_dir, participant_selection, output_dir,\n")
            if default_dir_fn != "nothing":
                new_lines.append(f"        default_output_dir_fn = {default_dir_fn},\n")
            if logs:
                new_lines.append("        init_log_fn = () -> begin\n")
                new_lines.extend(logs)
                new_lines.append("        end,\n")
            new_lines.append(f"        operation_name = \"{op_name}\")\n")
            
            continue
            
        if "log_file =" in line and "batch" not in line and "log" in line:
            # check if it's the `log_file = "..."` right before setup
            if i+1 < len(lines) and "setup_global_logging" in lines[i+1]:
                i += 1
                continue
            
        new_lines.append(line)
        i += 1
        
    with open(file, 'w') as f:
        f.writelines(new_lines)
    print(f"Refactored {file}")

