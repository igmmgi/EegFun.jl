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

    # Find the function definition and the body
    # They all start with `@log_call` or `@info`
    
    # We want to find lines starting with "    @log_call" or "    @info" that are NOT inside init_log_fn already.
    # And move them into batch_process.
    
    # Let's extract the block of @log_call and @info
    def_pattern = r"(function [^\n]+\n(?:(?:    [^\n]+\n)*?\)\n))"
    
    m = re.search(def_pattern, content)
    if not m: continue
    
    # We'll just collect all @log_call and @info that are standalone
    # and remove them from their original place, then put them in init_log_fn
    
    # Find all standalone @log_call and @info
    logs = []
    def repl(match):
        logs.append(match.group(1))
        return ""
        
    new_content = re.sub(r"^(    @log_call[^\n]*\n|    @info[^\n]*\n)", repl, content, flags=re.MULTILINE)
    
    if not logs: continue
    
    # Now add init_log_fn to batch_process
    init_log_code = "        init_log_fn = () -> begin\n" + "".join(["    " + l for l in logs]) + "        end,\n"
    
    new_content = re.sub(r"(batch_process\(.*?;)", r"\1\n" + init_log_code, new_content, flags=re.DOTALL)
    
    with open(file, 'w') as f:
        f.write(new_content)

print("Fixed files")
