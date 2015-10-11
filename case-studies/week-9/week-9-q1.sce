clear;
clf;

//==============================================================================
// Compute element stiffness matrix E_{mn} of a element
// \param dxidx dxi/dx of element
// \param thermal conductivity of model
//==============================================================================
function [E11, E12, E22] = StiffnessMatrix(dxidx, k)
  // for clarity
  J = abs(1 / dxidx);
  E11 = J * (dxidx * dxidx * k + 1 / 3);
  E12 = J * (-dxidx * dxidx * k + 1 / 6);
  E22 = J * (dxidx * dxidx * k + 1 / 3);
endfunction

// problem parameters
x_max = 1;  // max domain length
k = 1;

// simulation parameters
num_elem = 5;  // number of elements
dx = x_max / num_elem;
dx_ana = x_max / 100;
dxidx = 1 / (x_max / num_elem);
x = 0:dx:x_max;
x_ana = 0:dx_ana:x_max;

u_ana = exp(1) / (exp(2) - 1) * (exp(x_ana) - exp(-x_ana));

[E11, E12, E22] = StiffnessMatrix(dxidx, k);
// creating matrix K in KU = F
K_short_diag = E12 * ones(num_elem, 1);
K_main_diag = 2 * E11 * ones(num_elem + 1, 1);
K_main_diag(1) = K_main_diag(1) * 0.5;
K_main_diag($) = K_main_diag($) * 0.5;
K = diag(K_main_diag, 0) + diag(K_short_diag, -1) + diag(K_short_diag, 1);
K2 = K;
// essential boundary conditions, u(0) = 0, u(1) = 1, i.e., U_first = 0,
// U_last = 1
K(1, :) = [1, zeros(1, num_elem)];
K($, :) = [zeros(1, num_elem), 1];
F = [zeros(num_elem, 1); 1];

U = K \ F;
flux = K2 * U;

plot(x_ana, u_ana);
plot(x', U, 'k.');
xlabel("x position");
ylabel("Temperature");
title("1D heat equation (k = 1)");
legend(["Analytical solution", "FEM solution"], -1);
printf("Flux at first node: %f\nFlux at last node: %f", flux(1), flux($));
