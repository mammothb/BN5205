clear;
clf;
// constant variables
frequency = 0.1;  // Hz
omega = 2 * %pi * frequency;
sigma_0 = 6;  // MPa
E = 3;  // elastic modulus (MPa)
eta = 5;  // viscosity (MPa/s)
dt = 0.25;
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
for t = 1:length(time) - 1
  fe(t + 1) = ForwardEuler(time(t), fe(t), dt);
  hm(t + 1) = HeunsMethod(time(t), hm(t), dt);
  mm(t + 1) = MidpointMethod(time(t), mm(t), dt);
  fe_nc(t + 1) = ForwardEuler_nc(time(t), fe_nc(t), dt);
  hm_nc(t + 1) = HeunsMethod_nc(time(t), hm_nc(t), dt);
  mm_nc(t + 1) = MidpointMethod_nc(time(t), mm_nc(t), dt);
end

filename = pwd() + "\week-3\exp_data.txt";
data = fscanfMat(filename);

plot(data(:,1), data(:,2));
//plot(fe, Stress(time));
plot(fe_nc, Stress(time));
//plot(hm, Stress(time), 'r');
plot(hm_nc, Stress(time), 'r');
//plot(mm, Stress(time), 'c');
plot(mm_nc, Stress(time), 'c');
