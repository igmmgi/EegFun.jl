
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
