
function check_files_exist(conditions, filetype)
    problem = false
    for condition in conditions
        fname = "$(condition)_$(filetype).jld2"
        if !isfile(fname)
            println("File: $(fname) does not exist")
            problem = true
        end
    end
    return problem
end


function check_files_exist(subjects, conditions, filetype)
    problem = false
    for subject in subjects
        for condition in conditions
            fname = "$(subject)_$(condition)_$(filetype).jld2"
            if !isfile(fname)
                println("File: $(fname) does not exist")
                problem = true
            end
        end
    end
    return problem
end


function channel_number_to_channel_label(channel_labels, channel_numbers::Int64)
    return [channel_labels[channel_numbers]]
end

function channel_number_to_channel_label(channel_labels, channel_numbers::Vector{Int64})
    return channel_labels[channel_numbers]
end

datarange(x) = -(-(extrema(x)...))

colmeans(df::DataFrame, cols) = reduce(+, eachcol(df[!, cols])) ./ length(cols)
colmeans(df::Matrix) = reduce(+, eachrow(df)) ./ size(df)[1]
colmeans(df::Matrix, cols) = reduce(+, eachrow(df[:, cols])) ./ size(df)[1]


function consecutive(f, A::AbstractVector; step = 1)
    [f(A[i+step], A[i]) for i = 1:length(A)-step]
end


function splitgroups(v)
    start = 1
    start_idx::Vector{Int64} = []
    end_idx::Vector{Int64} = []
    for stop in [findall(diff(v) .> 1); lastindex(v)]
        push!(start_idx, v[start])
        push!(end_idx, v[stop])
        start = stop + 1
    end
    start_idx, end_idx
end
