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
// The function to integrate
// f(x) = x^3 + 3*x^2 + 2x - 1
// \param x x-value for the function
// \return y the value of the function evaluated at point x
//==============================================================================
function y = Func1(x)
  y = x * x * x + 3 * x * x + 2 * x - 1;
endfunction

//==============================================================================
// The function to integrate
// f(x) = x^2
// \param x x-value for the function
// \return y the value of the function evaluated at point x
//==============================================================================
function y = Func2(x)
  y = x * x;
endfunction

//==============================================================================
// The function to integrate (assumed linear) f(x) given that
// f(0) = 1 and f(1) = 2
// \param x x-value for the function
// \return y the value of the function evaluated at point x
//==============================================================================
function y = Func3(x)
  y = 1 * Psi1D(1, x) + 2 * Psi1D(2, x);
endfunction

//==============================================================================
// Returns the Gaussian Quadrature points and their weights. Up to 4 points
// \param num_points number of points
// \return position vector containing the point positions
// \return weight vector containing the point weights
//==============================================================================
function [position, weight] = GaussianQuadrature1D(num_points)
  positions = zeros(4);
  weights = zeros(4);

  positions(1, 1) = 0.5;
  positions(2, 1) = 0.5 - 0.5 / sqrt(3);
  positions(2, 2) = 0.5 + 0.5 / sqrt(3);
  positions(3, 1) = 0.5 - 0.5 * sqrt(0.6);
  positions(3, 2) = 0.5;
  positions(3, 3) = 0.5 + 0.5 * sqrt(0.6);
  positions(4, 1) = 0.5 - sqrt(525 + 70 * sqrt(30)) / 70;
  positions(4, 2) = 0.5 - sqrt(525 - 70 * sqrt(30)) / 70;
  positions(4, 3) = 0.5 + sqrt(525 - 70 * sqrt(30)) / 70;
  positions(4, 4) = 0.5 + sqrt(525 + 70 * sqrt(30)) / 70;

  weights(1, 1) = 1;
  weights(2, 1) = 1 / 2;
  weights(2, 2) = 1 / 2;
  weights(3, 1) = 5 / 18;
  weights(3, 2) = 4 / 9;
  weights(3, 3) = 5 / 18;
  weights(4, 1) = (18 - sqrt(30)) / 72;
  weights(4, 2) = (18 + sqrt(30)) / 72;
  weights(4, 3) = (18 + sqrt(30)) / 72;
  weights(4, 4) = (18 - sqrt(30)) / 72;

  position = positions(num_points, 1 : num_points);
  weight = weights(num_points, 1 : num_points);
endfunction

//==============================================================================
// Performs numerical integration with Gauss-Legendre quadrature
//==============================================================================
function value = GaussianIntegrate1D(a, b, f, num_points)
  [position, weight] = GaussianQuadrature1D(num_points);
  value = 0;
  for i = 1 : num_points
    value = value + weight(i) * f(a + (b - a) * position(i));
  end
  value = value * (b - a);
endfunction

exact_1 = intg(0, 1, Func1);  // Exact calculation
gauss1_1 = GaussianIntegrate1D(0, 1, Func1, 1);  // Using 1 point
gauss2_1 = GaussianIntegrate1D(0, 1, Func1, 2);  // Using 2 points
gauss3_1 = GaussianIntegrate1D(0, 1, Func1, 3);  // Using 3 points
gauss4_1 = GaussianIntegrate1D(0, 1, Func1, 4);  // Using 4 points

printf('f(x) = x^3 + 3*x^2 + 2x - 1\n');
printf('Exact: %.4f\nOne point: %.4f\nTwo points: %.4f\nThree points: %.4f\n'...
    + 'Four points: %.4f\n', exact_1, gauss1_1, gauss2_1, gauss3_1, gauss4_1);

exact_2 = intg(0, 2, Func2);  // Exact calculation
gauss2_2 = GaussianIntegrate1D(0, 2, Func2, 2);  // Using 2 points
printf('f(x) = x^2\n');
printf('Exact: %.4f\nTwo points: %.4f\n', exact_2, gauss2_2);

gauss2_3 = GaussianIntegrate1D(0, 1, Func3, 2);  // Using 2 points
printf('f(0) = 1 and f(1) = 2\n');
printf('Two points: %.4f\n', gauss2_3);
