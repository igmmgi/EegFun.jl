function circle_mask!(dat, grid_scale)
    for col = 1:size(dat)[1]
        for row = 1:size(dat)[2]
            xcentre, ycenter = (grid_scale / 2) - col, (grid_scale / 2) - row
            if sqrt((xcentre^2 + ycenter^2)) > (grid_scale / 2)
                dat[col, row] = NaN
            end
        end
    end
end


# TODO: which type of interpolation method as default?
# Compare to FieldTrip/MNE?
function data_interpolation_topo(dat, points, grid_scale)
    radius = 88 # mm
    x = y = range(-radius, radius, length = grid_scale)
    X, Y = repeat(x, outer = length(x))[:], repeat(y, inner = length(y))[:]
    grid = [X Y]'
    dat = interpolate(Multiquadratic(), points, dat)
    dat = ScatteredInterpolation.evaluate(dat, grid)
    dat = reshape(dat, grid_scale, grid_scale)
    circle_mask!(dat, grid_scale)
    return dat
end
