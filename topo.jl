
function circle_mask!(dat, grid_scale)
  for col in 1:size(dat)[1]
    for row in 1:size(dat)[2]
      xcentre = (grid_scale / 2) - col
      ycenter = (grid_scale / 2) - row
      if sqrt((xcentre^2 + ycenter^2)) > (grid_scale / 2)
        dat[col, row] = NaN
      end
    end
  end
end


function data_interpolation_topo(dat, points; grid_scale=300)

  radius = 88 # mm
  x = y = range(-radius, radius, length=grid_scale)
  X, Y = repeat(x', grid_scale)[:], repeat(y', grid_scale)[:]
  grid = [X Y]'
  dat = interpolate(Multiquadratic(), points, dat)
  dat = ScatteredInterpolation.evaluate(dat, grid)
  dat = reshape(dat, grid_scale, grid_scale)
  circle_mask!(dat, grid_scale)
  return dat
end



