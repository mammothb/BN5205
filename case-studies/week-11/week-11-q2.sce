clear;
//==============================================================================
// Evaluates 1D linear Lagrange basis functions. Throws error message if local
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
// The function to integrate (assumed linear) f(x) given that
// f(0) = 1 and f(1) = 2
// \param x x-value for the function
// \return y the value of the function evaluated at point x
//==============================================================================
function y = E11(x)
  y = (1 - Psi1D(2, x)) * (1 - Psi1D(2, x)) + (1 - Psi1D(1, x)) *...
      (1 - Psi1D(1, x));
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

function mat = GradMat(xi1, xi2)
  mat = [-1 + xi2, -1 + xi1;
         1 - xi2, -xi1;
         -xi2, 1 - xi1;
         xi2, xi1];
endfunction

num_gauss_points = 4; // 2x2
gp1 = 0.5 - 0.5 / sqrt(3);
gp2 = 0.5 + 0.5 / sqrt(3);
gauss_pos = [gp1, gp1; gp2, gp1; gp1, gp2; gp2, gp2];
gauss_weight = [0.25; 0.25; 0.25; 0.25];

answ = 0;
for n = 1 : 4
  xi1 = gauss_pos(n, 1);
  xi2 = gauss_pos(n, 2);
  Bmat = GradMat(xi1, xi2);
  answ = answ + gauss_weight(n) .* Bmat * Bmat';
end
disp(answ);

aaa = [1, 2; 3, 4];
disp(aaa(2, 1));
