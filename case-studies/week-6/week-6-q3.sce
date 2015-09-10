clear;
clf;

// model parameters
v_eff = 1e-5;  // cm / s, effective velocity
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
  // RHS to simplify equation, 2 * tau * f
  f = v_eff * dt / dx;
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
  solution = zeros(length(len), length(time));
  solution(1, :) = C_max;

  for i = 1:length(len) - 2
    A_main_diag(i) = 2 * (1 + delta);
    B_main_diag(i) = 2 * (1 - delta);
    if i < length(len) - 2 then
      A_shorter_diag(i) = -delta;
      B_shorter_diag(i) = delta;
    end
  end  // i
  // create A and B tri-diagonal matrices
  A = diag(A_main_diag, 0) + diag(A_shorter_diag, 1) + diag(A_shorter_diag, -1);
  B = diag(B_main_diag, 0) + diag(B_shorter_diag, 1) + diag(B_shorter_diag, -1);
  for t = 1:length(time) - 1
    // update RHS at every time step
    for x = 2:length(len) - 2
      C(x - 1) = f * (solution(x + 1, t) - solution(x - 1, t));
    end
    //boundary conditions
    C(1) = C(1) + delta * (solution(1, t) + solution(1, t + 1));
    solution(2:$-1, t + 1) = A \ B * solution(2:$-1, t) + A \ C;
  end
  // normalized to C_max
  plot(len', solution(:, 5 / dt + 1) / C_max);
  plot(len', solution(:, 15 / dt + 1) / C_max, 'r-');
  plot(len', solution(:, 25 / dt + 1) / C_max, 'g-');
  plot(len', solution(:, 35 / dt + 1) / C_max, 'm-');
  xlabel("$Distance\ along\ channel\ x$", "fontsize", 3);
  ylabel("$Concentration\ normalized\ to\ C_{max}$", "fontsize", 3);
  title("Leukocyte concentration profile");
  legend(["t = 5s", "t = 15s", "t = 25s", "t = 35s"], -1);
end
