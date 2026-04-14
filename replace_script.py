import os, re

substitutions = [
    # _, idx = findmin(abs.(A .- B)) -> idx = find_closest_time_index(A, B)
    (r"_\,\s*([A-Za-z0-9_]+)\s*=\s*findmin\(abs\.\(([^\)]+?)\s*\.-\s*([^\)]+?)\)\)", r"\1 = find_closest_time_index(\2, \3)"),
    # argmin(abs.(A .- B)) -> find_closest_time_index(A, B)
    (r"argmin\(abs\.\(([^\) ](?:(?:[^()]*|\([^()]*\))*)\s*)\.-\s*([^\) ](?:(?:[^()]*|\([^()]*\))*)\s*)\)\)", r"find_closest_time_index(\1, \2)"),
    # argmin(abs.(A)) -> find_closest_time_index(A, 0.0)
    (r"argmin\(abs\.\(([^\)]+)\)\)", r"find_closest_time_index(\1, 0.0)"),
    # findmin(abs.(A .- B))[2]
    (r"findmin\(abs\.\(([^\) ](?:(?:[^()]*|\([^()]*\))*)\s*)\.-\s*([^\) ](?:(?:[^()]*|\([^()]*\))*)\s*)\)\)\[2\]", r"find_closest_time_index(\1, \2)")
]

modified_files = []
for root, dirs, files in os.walk("/home/ian/Documents/Julia/EegFun.jl/src"):
    for file in files:
        if file.endswith(".jl"):
            path = os.path.join(root, file)
            with open(path, "r") as f:
                content = f.read()
            
            new_content = content
            
            # Custom simple regexes because the nested parentheses one is too complex for standard python re without regex module
            new_content = re.sub(r"_\,\s*([A-Za-z0-9_]+)\s*=\s*findmin\(abs\.\(([^)]+)\s*\.-\s*([^)]+)\)\)", r"\1 = find_closest_time_index(\2, \3)", new_content)
            new_content = re.sub(r"argmin\(abs\.\(([^)]+)\s*\.-\s*([^)]+)\)\)", r"find_closest_time_index(\1, \2)", new_content)
            new_content = re.sub(r"argmin\(abs\.\(([^)]+)\)\)", r"find_closest_time_index(\1, 0.0)", new_content)
            new_content = re.sub(r"findmin\(abs\.\(([^)]+)\s*\.-\s*([^)]+)\)\)\[2\]", r"find_closest_time_index(\1, \2)", new_content)

            if new_content != content:
                with open(path, "w") as f:
                    f.write(new_content)
                modified_files.append(path)

for m in modified_files:
    print(m)

