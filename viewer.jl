# TODO:. ENV["TERM_PROGRAM"] == "vscode"
function viewer(dat)
    vscodedisplay(dat)
end

function head(dat::DataFrame, n=nothing)
    isnothing(n) && (n=5)
    viewer(dat[1:n, :])
end

function viewer(dat::SingleDataFrameEeg)
    viewer(data(dat))
end

function head(dat::SingleDataFrameEeg, n=nothing)
    isnothing(n) && (n=5)
    viewer(data(data)[1:n, :])
end

# function viewer(dat:MultiDataFrameEeg)
#     viewer(to_data_frame(dat))
# end
# 
# function head(dat::MultiDataFrameEeg, n=nothing)
#     isnothing(n) && (n=5)
#     viewer(to_data_frame(dat)[1:n, :])
# end
# 
# function viewer(dat::Vector{EpochData})
#     viewer(to_data_frame(dat))
# end
# 
# function head(dat::Vector{EpochData}, n=nothing)
#     isnothing(n) && (n=5)
#     viewer(to_data_frame(dat)[1:n, :])
# end 
