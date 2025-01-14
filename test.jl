using Random, VoronoiDelaunay, Plots

import Base.==

struct MyEdge{T<:AbstractPoint2D}
    _a::T
    _b::T
end

==(e1::MyEdge{T}, e2::MyEdge{T}) where {T<:AbstractPoint2D} = ((e1._a === e2._a) && (e1._b === e2._b)) || ((e1._b === e2._a) && (e2._b === e1._a))

###==(p1::T, p2::T) where {T<:AbstractPoint2D} = (getx(p1) == getx(p2)) && (gety(p1) == gety(p2))

### Create a Delaunay tesselation from random points
tess = DelaunayTessellation2D(10)

for _ in 1:5
    push!(tess, Point2D(rand()+1, rand()+1))
end

edges = MyEdge[]
function add_edge!(edges, edge)
    i = findfirst(e -> e == edge, edges)

    if isnothing(i) # not found
        push!(edges, edge)
    else # found so not outer, remove this edge
        deleteat!(edges, i) 
    end
end

for trig in tess
    a, b, c = geta(trig), getb(trig), getc(trig)
    add_edge!(edges, MyEdge(b, c))
    add_edge!(edges, MyEdge(a, b))
    add_edge!(edges, MyEdge(a, c))
end

### PLOT
x, y = Float64[], Float64[] # outer edges
for edge in edges
    push!(x, getx(edge._a))
    push!(x, getx(edge._b))
    push!(x, NaN)
    push!(y, gety(edge._a))
    push!(y, gety(edge._b))
    push!(y, NaN)
end

xall, yall = getplotxy(delaunayedges(tess)) # all the edges

plot(xall, yall, color=:blue, fmt=:svg, size=(400,400))
plot!(x, y, color=:red, linewidth=3, opacity=0.5)
