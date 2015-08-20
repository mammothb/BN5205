clear;
clf;
// constant variables
frequency = 0.1;  // Hz
omega = 2 * %pi * frequency;
sigma_0 = 6;  // MPa
E = 3;  // elastic modulus (MPa)
eta = 5;  // viscosity (MPa/s)
dt = 0.01;
dt10 = dt * 10;
dt100 = dt * 100;

function sigma_t = Stress(t)
  sigma_t = sigma_0 * sin(omega * t);
endfunction

function epsilon = Analytical(t)
  epsilon = sigma_0 / (E^2 + eta^2 * omega^2) * (omega * eta *...
      exp(-E / eta * t) + E * sin(omega * t) - omega * eta *...
      cos(omega * t));
endfunction

// Forward Euler using 3 different values of dt
time = [0:dt:100];
time10 = [0:dt10:100];
time100 = [0:dt100:100];
epsilon_fe_1 = zeros(time);
epsilon_fe_10 = zeros(time10);
epsilon_fe_100 = zeros(time100);
for t = 1:length(time) - 1
  epsilon_fe_1(t + 1) = epsilon_fe_1(t) + dt / eta * (-E *...
      epsilon_fe_1(t) + Stress(time(t)));
  if (~modulo(t, 10))
    t10 = t / 10;
    epsilon_fe_10(t10 + 1) = epsilon_fe_10(t10) + dt10 / eta *...
        (-E * epsilon_fe_10(t10) + Stress(time10(t10)));
  end
  if (~modulo(t, 100))
    t100 = t / 100;
    epsilon_fe_100(t100 + 1) = epsilon_fe_100(t100) + dt100 /...
        eta * (-E * epsilon_fe_100(t100) + Stress(time100(t100)));
  end    
end

plot(time, epsilon_fe_1);
plot(time10, epsilon_fe_10, 'b-');
plot(time100, epsilon_fe_100, 'g-');
plot(time, Analytical(time), 'r-');
