clear;

//==============================================================================
// Evaluetes 2D linear Lagrange basis functions. Throws error message if local
// node number num or location of interest xi is out of range
// \param num local node (basis function) number
// \param xi1 location of interest 1
// \param xi2 location of interest 2
// \return value of basis function
//==============================================================================
function value = psi2D(num, xi1, xi2)
  if xi1 >= 0 & xi1 <= 1 & xi2 >= 0 & xi2 <= 1 then
    if num == 1 then
      value = (1 - xi1) * (1 - xi2);
    elseif num == 2 then
      value = xi1 * (1 - xi2);
    elseif num == 3 then
      value = (1 - xi1) * xi2;
    elseif num == 4 then
      value = xi1 * xi2;
    else
      error('Invalid local node (basis function) number.');
    end
  else
    error('Invalid xi value.');
  end
endfunction

//==============================================================================
// Interpolate the x-, y-coordinates of point of interest within the element.
// Assumes maximum four nodes in element, does 2D bilinear interpolation.
// \param nodePos array containing position for the element in the format
//        nodePos(n, d) which gives the x (d = 1) or y (d = 2) coordinate of the
//        nth node in the element
// \param xi1 location of interest 1
// \param xi2 location of interest 2
// \return interpolated coordinate of the element
//==============================================================================
function [x, y] = interpolate2D(nodePos, xi1, xi2)
  // Number of nodes in element
  numNodes = length(nodePos(:, 1));
  x = 0;
  y = 0;
  for i = 1:numNodes
    x = x + psi2D(i, xi1, xi2) * nodePos(i, 1);
    y = y + psi2D(i, xi1, xi2) * nodePos(i, 2);
  end
endfunction

nodePos = zeros(4, 2);
nodePos(1, :) = [2, 1];
nodePos(2, :) = [5, 0];
nodePos(3, :) = [2, 3];
nodePos(4, :) = [5, 4];

[x, y] = interpolate2D(nodePos, 0.5, 0.5);
printf("%f, %f\n", x, y);
