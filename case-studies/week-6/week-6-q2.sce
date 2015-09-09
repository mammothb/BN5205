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

function slope = dC_idt(C_prev, C_now, C_next)
  slope = D * (C_next - 2 * C_now + C_prev) / dx / dx - v_eff * (C_next -...
      C_prev) / 2 / dx;
endfunction

if delta > 0.5 then
  disp("Failed stability condition");
else
  f = v_eff * dt / 2 / dx;  // RHS to simplify equation
  // coefficients for discretized form of diffusion-convection equation after
  // rearranging
  a = [delta + f, 1 - 2 * delta, delta - f];
  time = [0:dt:35];
  len = [0:dx:x_max];
  C = zeros(length(len), length(time));
  // initial conditions
  C(:, 1) = C_max;
  for t = 1:length(time) - 1
    for x = 2:length(len) - 1
      C(x, t + 1) = a(1) * C(x - 1, t) + a(2) * C(x, t) + a(3) * C(x + 1, t);
    end  // x
  end  // t

  // allocate memory for lines method solution
  C_l = zeros(length(len), length(time));
  // initial conditions
  C_l(:, 1) = C_max;
  for t = 1:length(time) - 1
    for x = 2:length(len) - 1
      C_l(x, t + 1) = C_l(x, t) + dt * dC_idt(C_l(x - 1, t), C_l(x, t),...
          C_l(x + 1, t));
    end  // x
  end  // t

  // normalized to C_max
  plot(len', C(:, 15 / dt + 1) / C_max);
  plot(len', C_l(:, 15 / dt + 1) / C_max, 'r-');
  xlabel("$Distance\ along\ channel\ x$", "fontsize", 3);
  ylabel("$Concentration\ normalized\ to\ C_{max}$", "fontsize", 3);
  title("$Leukocyte\ concentration\ profile\ at\ t = 15s$");
  legend(["FTCS", "Method of lines"], -1);
end
