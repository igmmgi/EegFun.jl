import re
import glob

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

for file in files:
    with open(file, 'r') as f:
        content = f.read()

    # Find the function definition
    # The batch function is typically defined after its docs, but before _process_...
    # We can match `batch_process(process_fn` to find where the batch_process is called.
    
    # We want to replace `@log_call` and `@info` that are at indentation 4 spaces inside the MAIN function.
    
    # Let's split content into lines
    lines = content.split('\n')
    
    # Find the main batch_process call. It usually starts with "    batch_process(process_fn"
    batch_idx = -1
    for i, line in enumerate(lines):
        if line.startswith("    batch_process(process_fn"):
            batch_idx = i
            break
            
    if batch_idx == -1:
        continue
        
    # Walk backward from batch_process up to find the function signature end.
    # The signature end is usually ")\n".
    # And we collect the @log_call and @info
    
    logs = []
    
    # We will search the lines BEFORE batch_idx for `    @log_call` and `    @info`
    # BUT only ones that are NOT inside other functions. 
    # Actually, process_fn is defined before batch_process.
    # So we should search from the function signature down to process_fn
    
    # Let's find process_fn
    process_fn_idx = -1
    for i in range(batch_idx - 1, -1, -1):
        if line.startswith("    process_fn ="):
            process_fn_idx = i
            break
            
    # Search for logs before process_fn
    new_lines = []
    for i in range(len(lines)):
        line = lines[i]
        
        # if we are inside the main function (approx before batch_process)
        if i < batch_idx and (line.startswith("    @log_call") or line.startswith("    @info")):
            logs.append(line)
        else:
            new_lines.append(line)
            
    # Now insert init_log_fn into batch_process call
    # batch_process spans multiple lines. We want to insert it as a kwarg.
    # We can find `batch_process(process_fn` and insert it after `file_pattern;`
    
    res = "\n".join(new_lines)
    
    if logs:
        init_code = "        init_log_fn = () -> begin\n" + "\n".join(["        " + l.strip() for l in logs]) + "\n        end,\n"
        res = re.sub(r"(batch_process\(process_fn, file_pattern;)", r"\1\n" + init_log_code, res)
        
    with open(file, 'w') as f:
        f.write(res)

print("Fixed files 2")
