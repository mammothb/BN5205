clear;
clf;
// constant variables
frequency = 0.1;  // Hz
omega = 2 * %pi * frequency;
sigma_0 = 6;  // MPa
E = 3;  // elastic modulus (MPa)
eta = 5;  // viscosity (MPa/s)
dt1 = 0.01;
dt2 = dt1 * 10;
dt3 = dt1 * 100;

function sigma_t = Stress(t)
  sigma_t = sigma_0 * sin(omega * t);
endfunction

function epsilon = Analytical(t)
  epsilon = sigma_0 / (E^2 + eta^2 * omega^2) * (omega * eta *...
      exp(-E / eta * t) + E * sin(omega * t) - omega * eta *...
      cos(omega * t));
endfunction

// Time domains for 3 different values of dt
time1 = [0:dt1:100];
time2 = [0:dt2:100];
time3 = [0:dt3:100];
// Forward Euler using 3 different values of dt
fe_1 = zeros(time1);
fe_2 = zeros(time2);
fe_3 = zeros(time3);
// RK2 Heun's method using 3 different values of dt
hm_1 = zeros(time1);
hm_2 = zeros(time2);
hm_3 = zeros(time3);
// RK2 Midpoint method using 3 different values of dt
mm_1 = zeros(time1);
mm_2 = zeros(time2);
mm_3 = zeros(time3);
for t1 = 1:length(time1) - 1
  fe_1(t1 + 1) = fe_1(t1) + dt1 / eta * (-E * fe_1(t1) +...
      Stress(time1(t1)));

  // Heun's method uses Forward Euler as the intermediate value
  hm_1(t1 + 1) = hm_1(t1) + dt1 / 2 / eta * (-E * hm_1(t1) +...
      Stress(time1(t1)) + -E * fe_1(t1 + 1) +...
      Stress(time1(t1 + 1)));

  mm_intermediate = mm_1(t1) + dt1 / 2 / eta * (-E * mm_1(t1) +...
      Stress(time1(t1)));
  mm_1(t1 + 1) = mm_1(t1) + dt1 / eta * (-E * mm_intermediate +...
      Stress(time1(t1) + dt1 / 2));
  if (~modulo(t1, 10))
    t2 = t1 / 10;
    fe_2(t2 + 1) = fe_2(t2) + dt2 / eta * (-E * fe_2(t2) +...
        Stress(time2(t2)));

    hm_2(t2 + 1) = hm_2(t2) + dt2 / 2 / eta * (-E * hm_2(t2) +...
    Stress(time2(t2)) + -E * fe_2(t2 + 1) +...
    Stress(time2(t2 + 1)));

    mm_intermediate = mm_2(t2) + dt2 / 2 / eta * (-E * mm_2(t2) +...
        Stress(time2(t2)));
    mm_2(t2 + 1) = mm_2(t2) + dt2 / eta * (-E * mm_intermediate +...
        Stress(time2(t2) + dt2 / 2));
  end
  if (~modulo(t1, 100))
    t3 = t1 / 100;
    fe_3(t3 + 1) = fe_3(t3) + dt3 / eta * (-E * fe_3(t3) +...
        Stress(time3(t3)));

    hm_3(t3 + 1) = hm_3(t3) + dt3 / 2 / eta * (-E * hm_3(t3) +...
    Stress(time3(t3)) + -E * fe_3(t3 + 1) +...
    Stress(time3(t3 + 1)));

    mm_intermediate = mm_3(t3) + dt3 / 2 / eta * (-E * mm_3(t3) +...
        Stress(time3(t3)));
    mm_3(t3 + 1) = mm_3(t3) + dt3 / eta * (-E * mm_intermediate +...
        Stress(time3(t3) + dt3 / 2));
  end
end

plot(time1, Analytical(time1));
plot(time1, fe_1, 'm');
plot(time2, fe_2, 'g');
plot(time3, fe_3, 'k');
plot(time1, hm_1, 'm--');
plot(time2, hm_2, 'g--');
plot(time3, hm_3, 'k--');
plot(time1, mm_1, 'm:');
plot(time2, mm_2, 'g:');
plot(time3, mm_3, 'k:');
xlabel("$Time\ t$", "fontsize", 3)
ylabel("$Strain\ \epsilon$", "fontsize", 3);
title("Viscoelastic model of IVD (Analytical vs Forward Euler vs"...
    + " RK2)");
legend(["Analytical"; "FE, dt = 0.01"; "FE, dt = 0.1";...
    "FE, dt = 1"; "HM, dt = 0.01"; "HM, dt = 0.1";...
    "HM, dt = 1"; "MM, dt = 0.01"; "MM, dt = 0.1";...
    "MM, dt = 1"]);
