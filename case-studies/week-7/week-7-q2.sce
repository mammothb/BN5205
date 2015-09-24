clear;
clf;

//==============================================================================
// This function creates the penta-diagonal matrix M used in the FD method
// For a problem with nx * ny internal nodes, matrix M has a size of n * n,
// where n = nx * ny. M is a tri-diagonal block matrix:
// A I O
// I A I
// O I A
// where I is a nx * nx identity matrix and O is a nx * nx zero matrix, A is a
// nx * nx matrix that looks like:
// -4  1  0
//  1 -4  1
//  0  1 -4
// when nx = 3
//==============================================================================
function M = MakePentaDiagonalMatrix(nx, ny)
  // main diagonal of -4's
  A_main_diag = -4 * ones(ny, 1);
  // only one created since they are identical
  A_sub_diag = ones(nx - 1, 1);
  temp_sub_diag = [A_sub_diag; 0];
  // maind diagonal of identity matrix
  I_diag = ones(nx, 1);
  // replicate A_main_diag in M_main_diag nx times, based on relation between M
  // and A
  M_main_diag = repmat(A_main_diag, nx, 1);
  M_sub_diag = repmat(temp_sub_diag, ny - 1, 1);
  M_sub_diag = [M_sub_diag; A_sub_diag];
  M_sub_sub_diag = repmat(I_diag, ny - 1, 1);
  // creates the Penta-diagonal matrix from diagonals
  M = diag(M_main_diag, 0) + diag(M_sub_diag, 1) + diag(M_sub_diag, -1) +...
      diag(M_sub_sub_diag, nx) + diag(M_sub_sub_diag, -nx);
endfunction

// model parameters
Q = 4e14;  // J/(s*m^3), Heat energy generated
k = 0.9e7;  // m^2/s, Thermal diffusivity
rho = 1e-3;  // kg/m^3, Tissue density
c = 3.7e3;  // J/(kg*K), Specific thermal capacity
L = 3.0 / 100;  // m, Length of tissue sample (AB)
H = 4.0 / 100;  // m, Height of tissue sample (AD)
z = 1.5 / 100;  // m, distance of cryoablation line (z)
d = 0.5 / 100;  // m, distance of cyroablation line (d)
T_AB = 37.0 + 273.15;  // K, temperature of along AB
T_DC = 20.0 + 273.15;  // K, temperature of along AB
T_grad = (20.0 - 37.0) / H;  // K/m, temperature gradient along AD and BC
T_side = T_AB;  // K, temperature of sides at bottom
T_peak = -20.0 + 273.15;  // K, peak minimum temperature

// simulation parameters
h = 0.1 / 100;  // m, space step
x = h:h:L;
y = h:h:H;
xplot = 0:h:L;  // x axis for plotting
yplot = 0:h:H;  // x axis for plotting
num_nodes_x = L / h + 1;  // number of nodes in x direction
num_nodes_y = H / h + 1;  // number of nodes in y direction
nx = num_nodes_x - 2;  // number of internal nodes in x direction
ny = num_nodes_y - 2;  // number of internal nodes in y direction
mat_size = nx * ny;  // var for matrix size

// matrices and vectors
solution = zeros(num_nodes_x, num_nodes_y);
q = zeros(mat_size, 1);  // source/sink term for internal nodes
u = zeros(mat_size, 1);  // vector to store solution for internal nodes
rhs = zeros(mat_size, 1);  // vector for RHS
M = MakePentaDiagonalMatrix(nx, ny);

// set temperature at top and bottom
solution(:, 1) = T_AB;
solution(:, $) = T_DC;

for j = 1:ny
  // Determind temperature at sides using Forward Euler
  T_side = T_side + h * T_grad;
  // set temperature at sides
  solution(1, j + 1) = T_side;
  solution($, j + 1) = T_side;
  for i = 1:nx
    n = (j - 1) * nx + i;
    if y(j) >= d & y(j) <= H - d & x(i) == z then
      q(n) = Q / c / rho / k;
    end
    if j == 1 then
      rhs(n) = rhs(n) - T_AB;
    end
    if j == ny then
      rhs(n) = rhs(n) - T_DC;
    end
    if i == 1 | i == nx then
      rhs(n) = rhs(n) - T_side;
    end
  end
end
rhs = h * h * q + rhs;
u = M \ rhs;
for j = 1:ny
  for i = 1:nx
    n = (j - 1) * nx + i;
    solution(i + 1, j + 1) = u(n);
  end
end

contour(xplot, yplot, solution, 20);
xlabel("x position");
ylabel("y position");
title("Temperature in tissue");
