
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
