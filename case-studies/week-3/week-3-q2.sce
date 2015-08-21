clear;
clf;
// read from file
filename = pwd() + "\case-studies\week-3\exp_data.txt";
data = fscanfMat(filename);
// constant variables
frequency = 0.1;  // Hz
omega = 2 * %pi * frequency;
sigma_0 = 6;  // MPa
E = 3;  // elastic modulus (MPa)
eta = 5;  // viscosity (MPa/s)
dt = 0.29;  // can reach 0.4 for 4th order Runge-Kutta
function sigma_t = Stress(t)
  sigma_t = sigma_0 * sin(omega * t);
endfunction

function elastic_modulus = E_nc(epsilon)
  elastic_modulus = E + 46 * epsilon^2;
endfunction

function slope = dedt(t, epsilon)
  slope = (-E * epsilon + Stress(t)) / eta;
endfunction

function slope = dedt_nc(t, epsilon)
  slope = (-E_nc(epsilon) * epsilon + Stress(t)) / eta;
endfunction

function y = ForwardEuler(t, y_prev, h)
  y = y_prev + h * dedt(t, y_prev);
endfunction

function y = HeunsMethod(t, y_prev, h)
  y = y_prev + h / 2 * (dedt(t, y_prev) + dedt(t + h, y_prev +...
      h * dedt(t, y_prev)));
endfunction

function y = MidpointMethod(t, y_prev, h)
  y = y_prev + h * dedt(t + h / 2, y_prev + h / 2 *...
      dedt(t, y_prev));
endfunction

function y = RK4(t, y_prev, h)
  k1 = dedt(t, y_prev);
  k2 = dedt(t + h / 2, y_prev + h / 2 * k1);
  k3 = dedt(t + h / 2, y_prev + h / 2 * k2);
  k4 = dedt(t + h, y_prev + h * k3);
  y = y_prev + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
endfunction

function y = ForwardEuler_nc(t, y_prev, h)
  y = y_prev + h * dedt_nc(t, y_prev);
endfunction

function y = HeunsMethod_nc(t, y_prev, h)
  y = y_prev + h / 2 * (dedt_nc(t, y_prev) + dedt_nc(t + h,...
      y_prev + h * dedt_nc(t, y_prev)));
endfunction

function y = MidpointMethod_nc(t, y_prev, h)
  y = y_prev + h * dedt_nc(t + h / 2, y_prev + h / 2 *...
      dedt_nc(t, y_prev));
endfunction

function y = RK4_nc(t, y_prev, h)
  k1 = dedt_nc(t, y_prev);
  k2 = dedt_nc(t + h / 2, y_prev + h / 2 * k1);
  k3 = dedt_nc(t + h / 2, y_prev + h / 2 * k2);
  k4 = dedt_nc(t + h, y_prev + h * k3);
  y = y_prev + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
endfunction

time = [0:dt:100];
// Finding strain using Forward Euler
fe = zeros(time);
fe_nc = zeros(time);
// Finding strain using Heun's method
hm = zeros(time);
hm_nc = zeros(time);
// Finding strain using Midpoint method
mm = zeros(time);
mm_nc = zeros(time);
// Finding strain using 4th order Runge-Kutta
rk = zeros(time);
rk_nc = zeros(time);
for t = 1:length(time) - 1
  fe(t + 1) = ForwardEuler(time(t), fe(t), dt);
  hm(t + 1) = HeunsMethod(time(t), hm(t), dt);
  mm(t + 1) = MidpointMethod(time(t), mm(t), dt);
  rk(t + 1) = RK4(time(t), rk(t), dt);
  fe_nc(t + 1) = ForwardEuler_nc(time(t), fe_nc(t), dt);
  hm_nc(t + 1) = HeunsMethod_nc(time(t), hm_nc(t), dt);
  mm_nc(t + 1) = MidpointMethod_nc(time(t), mm_nc(t), dt);
  rk_nc(t + 1) = RK4_nc(time(t), rk_nc(t), dt);
end

plot(data(:,1), data(:,2), 'k.');
//plot(fe, Stress(time));
plot(fe_nc, Stress(time), 'm');
//plot(hm, Stress(time), 'r');
plot(hm_nc, Stress(time), 'r');
//plot(mm, Stress(time), 'c');
plot(mm_nc, Stress(time), 'c');
//plot(rk, Stress(time), 'b');
plot(rk_nc, Stress(time), 'b');
xlabel("$Strain\ \epsilon$", "fontsize", 3)
ylabel("$Stress\ \sigma$", "fontsize", 3);
title("Hysteresis Loop");
legend(["Experiment"; "FE"; "HM"; "MM"; "RK4"], -1);
