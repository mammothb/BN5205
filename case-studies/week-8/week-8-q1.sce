clear;

//==============================================================================
// Evaluetes 1D linear Lagrange basis functions. Throws error message if local
// node number num or location of interest xi is out of range
// \param num local node (basis function) number
// \param xi location of interest
// \return value of basis function
//==============================================================================
function value = Psi1D(num, xi)
  if xi >= 0 & xi <= 1 then
    if num == 1 then
      value = 1 - xi;
    elseif num == 2 then
      value = xi;
    else
      error('Invalid local node (basis function) number.');
    end
  else
    error('Invalid xi value.');
  end
endfunction

//==============================================================================
// Perform 1D linear Lagrange interpolation for u = x^2 over x = [0, 2].
// \param numElems number of element in the doman
// \param pointOfInterest loation of point of interest
// \return interpolated value at point of interest
//==============================================================================
function value = Interpolate1D(num_elems, point_of_interest)
  domain_length = 2.0;
  // Length of each local space wrt global space
  dx = domain_length / num_elems;
  elem = floor(point_of_interest / dx);  // the element number which POI lies in
  u = zeros(2, 1);
  u(1) = elem * dx * elem * dx;  // first local node
  u(2) = (elem + 1) * dx * (elem + 1) * dx;  // second local node
  xi = point_of_interest / dx - elem;
  value = u(1) * Psi1D(1, xi) + u(2) * Psi1D(2, xi);
endfunction

// simulation parameters
point_of_interest = 1.2;
num_elems = 1;
exact = point_of_interest * point_of_interest;
err = 100.0;
while err >= 1.0 then
  result = Interpolate1D(num_elems, point_of_interest);
  num_elems = num_elems + 1;
  err = abs(result - exact) / exact * 100.0;
end
for i = 1:10
  result = Interpolate1D(i, point_of_interest);
  printf("Elements: %d, Error: %f\n", i, abs(result - exact) / exact * 100.0);
end
// 5 elements gives exact answer because POI happenes to be on a node
printf("%d\n", num_elems - 1);
