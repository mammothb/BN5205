clear;
clf;

// model parameters
v_eff = 0;  // cm / s, effective velocity
D = 1e-4;  // cm^2 / s, diffusion coefficient
C_max = 10e6;  // mM, max leukocyte concentration
x_max = 0.2;  // cm, length of channel

// solver parameters
dx = 0.01;  // cm, space step
dt = 0.05;  // s, time step

// check stability condition
delta = D * dt / dx / dx;

if delta > 0.5 then
  disp("Failed stability condition");
else
  // RHS to simplify equation
  f = v_eff * dt / 2 / dx;
  // coefficients for discretized form of diffusion-convection equation after
  // rearranging
  a = [delta + f, 1 - 2 * delta, delta - f];
  time = [0:dt:35];
  len = [0:dx:x_max];
  // forming the diagonals of the two tri-diagonal matrices, A and B
  A_main_diag = zeros(length(len) - 2, 1);
  A_upper_diag = zeros(length(len) - 3, 1);
  A_lower_diag = zeros(length(len) - 3, 1);
  B_main_diag = zeros(length(len) - 2, 1);
  B_upper_diag = zeros(length(len) - 3, 1);
  B_lower_diag = zeros(length(len) - 3, 1);
  C = zeros(length(len) - 2, 1);  // C updates every time interval

  for i = 1:length(len) - 2
    A_main_diag(i) = 2 * (1 + delta);
    B_main_diag(i) = 2 * (1 - delta);
    if i < length(len) - 2 then
      A_upper_diag(i) = -delta;
      A_lower_diag(i) = -delta;
      B_upper_diag(i) = delta;
      B_lower_diag(i) = delta;
    end
  end  // i

  // // initial conditions
  // C(1, :) = C_max;
  //
  // for t = 1:length(time) - 1
  //   for x = 2:length(len) - 1
  //     C(x, t + 1) = a(1) * C(x - 1, t) + a(2) * C(x, t) + a(3) * C(x + 1, t);
  //   end  // x
  // end  // t
  // // normalized to C_max
  // plot(len', C(:, 5 / dt + 1) / C_max);
  // plot(len', C(:, 15 / dt + 1) / C_max, 'r-');
  // plot(len', C(:, 25 / dt + 1) / C_max, 'g-');
  // plot(len', C(:, 35 / dt + 1) / C_max, 'm-');
  xlabel("$Distance\ along\ channel\ x$", "fontsize", 3);
  ylabel("$Concentration\ normalized\ to\ C_{max}$", "fontsize", 3);
  title("Leukocyte concentration profile");
  legend(["t = 5s", "t = 15s", "t = 25s", "t = 35s"], -1);
end
