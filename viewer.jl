function viewer(dat)
    ENV["TERM_PROGRAM"] == "vscode" ? vscodedisplay(dat) : display(dat)
end

function viewer(dat::EegData)
    viewer(data(dat))
end

function head(dat::EegData; n=nothing)
    isnothing(n) && (n=5)
    viewer(data(dat)[1:n, :])
end
