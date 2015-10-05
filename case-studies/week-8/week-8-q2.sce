clear;

//==============================================================================
// Evaluetes 2D linear Lagrange basis functions. Throws error message if local
// node number num or location of interest xi is out of range
// \param num local node (basis function) number
// \param xi1 location of interest 1
// \param xi2 location of interest 2
// \return value of basis function
//==============================================================================
function value = Psi2D(num, xi1, xi2)
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
function [x, y] = Interpolate2D(node_pos, xi1, xi2)
  // Number of nodes in element
  num_nodes = length(node_pos) / 2;
  x = 0;
  y = 0;
  for i = 1:num_nodes
    x = x + Psi2D(i, xi1, xi2) * node_pos(i, 1);
    y = y + Psi2D(i, xi1, xi2) * node_pos(i, 2);
  end
endfunction

//==============================================================================
// 2D cross product between two vectors
// \param vec1 first vector (x1, y1)
// \param vec2 second vector (x2, y2)
// \return magnitude of resultant vector from the cross product of (x1, y1, 0)
//         and (x2, y2, 0)
//==============================================================================
function value = Cross2D(vec1, vec2)
  value = vec1(1) * vec2(2) - vec1(2) * vec2(1);
endfunction

//==============================================================================
// Finds root of the polynomial A*(1-x)^2 + 2*B*x*(1-x) + C*x^2 = 0
// \param A first coefficient
// \param B second coefficient
// \param C third coefficient
// \return the root where 0 <= x <= 1
//==============================================================================
function root = FindRoot(A, B, C)
  // denominator to find root to quadratic polynomial
  den = A - 2 * B + C;
  // check for special case when polynomial is linear
  if abs(den) < 1e-20 then
    // denominator when polynomial is linear
    lin_den = A - C;
    if abs(lin_den) < 1e-20 then
      if abs(A) < 1e-20 then
        error('All values of root contains (x, y)');
      else
        error('No values of root contains (x, y)');
      end
    else
      root = A / lin_den;
    end
  else
    ans1 = ((A - B) + sqrt(B * B - A * C)) / den;
    ans2 = ((A - B) - sqrt(B * B - A * C)) / den;
    root = ans1 * (ans1 >= 0 & ans1 <= 1) + ans2 * (ans2 >= 0 & ans2 <= 1);
  end
endfunction

//==============================================================================
// Find xi1 and xi2 when given (x, y) and coordinates of the four local nodes of
// the element. Performs findRoot twice instead of subbing xi1 back to find xi2
// to avoid division by zero error when some edges are vertical
// \param nodePos array containing position for the element in the format
//        nodePos(n, d) which gives the x (d = 1) or y (d = 2) coordinate of the
//        nth node in the element
// \param x interpolated x coordinate
// \param y interpolated y coordinate
// \return xi1 and x2 of the element
//==============================================================================
function [xi1, xi2] = ReverseInterpolate2D(node_pos, x, y)
  // for clarity
  node = [x, y];
  a = node_pos(1, :) - node;
  b1 = node_pos(1, :) - node_pos(3, :);
  b2 = node_pos(1, :) - node_pos(2, :);
  c1 = node_pos(2, :) - node_pos(4, :);
  c2 = node_pos(3, :) - node_pos(4, :);
  d1 = node_pos(2, :) - node;
  d2 = node_pos(3, :) - node;
  A1 = Cross2D(a, b1);
  B1 = (Cross2D(a, c1) + Cross2D(d1, b1)) / 2;
  C1 = Cross2D(d1, c1);
  A2 = Cross2D(a, b2);
  B2 = (Cross2D(a, c2) + Cross2D(d2, b2)) / 2;
  C2 = Cross2D(d2, c2);
  xi1 = FindRoot(A1, B1, C1);
  xi2 = FindRoot(A2, B2, C2);
endfunction

node_pos = zeros(4, 2);
node_pos(1, :) = [1, 1];
node_pos(2, :) = [5, 1];
node_pos(3, :) = [1, 5];
node_pos(4, :) = [5, 5];

xi1 = 0.35;
xi2 = 0.81;

[x, y] = Interpolate2D(node_pos, xi1, xi2);
[ans1, ans2] = ReverseInterpolate2D(node_pos, x, y);
printf("Interpolated coordinates: (%f, %f)\n", x, y);
printf("Original xi\n xi1: %f, xi2: %f\n", xi1, xi2);
printf("Revere interpolated xi\n xi1: %f, xi2: %f\n", ans1, ans2);
