clear;

//==============================================================================
// Evaluetes 1D linear Lagrange basis functions. Throws error message if local
// node number num or location of interest xi is out of range
// \param num local node (basis function) number
// \param xi location of interest
// \return value of basis function
//==============================================================================
function value = psi1D(num, xi)
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
function value = interpolate1D(numElems, pointOfInterest)
  dx = 2 / numElems;  // Length of each local space wrt global space
  elem = floor(pointOfInterest / dx);  // the element number which POI lies in
  U = zeros(2, 1);
  U(1) = elem * dx * elem * dx;  // first local node
  U(2) = (elem + 1) * dx * (elem + 1) * dx;  // second local node
  xi = pointOfInterest / dx - elem;
  value = U(1) * psi1D(1, xi) + U(2) * psi1D(2, xi);
endfunction

// simulation parameters
pointOfInterest = 1.2;
numElems = 1;
exact = pointOfInterest * pointOfInterest;
err = 100.0;
while err >= 1.0 then
  result = interpolate1D(numElems, pointOfInterest);
  numElems = numElems + 1;
  err = abs(result - exact) / exact * 100.0;
end
printf("%d\n", numElems);
