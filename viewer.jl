# TODO:. ENV["TERM_PROGRAM"] == "vscode"
function viewer(dat)
    vscodedisplay(dat)
end

function head(dat::DataFrame, n=nothing)
    isnothing(n) && (n=5)
    viewer(dat[1:n, :])
end

function viewer(dat::Union{BioSemiBDF.BioSemiData, ContinuousData})
    viewer(dat.data)
end

function head(dat::Union{BioSemiBDF.BioSemiData, ContinuousData}, n=nothing)
    isnothing(n) && (n=5)
    viewer(dat.data[1:n, :])
end

function viewer(dat::EpochData)
    viewer(to_data_frame(dat))
end

function head(dat::EpochData, n=nothing)
    isnothing(n) && (n=5)
    viewer(to_data_frame(dat)[1:n, :])
end

function viewer(dat::Vector{EpochData})
    viewer(to_data_frame(dat))
end

function head(dat::Vector{EpochData}, n=nothing)
    isnothing(n) && (n=5)
    viewer(to_data_frame(dat)[1:n, :])
end 
