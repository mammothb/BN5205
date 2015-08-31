clear;
clf;
// constant variables
dt1 = 0.02;  // seconds
dt2 = dt1 * 10;  // seconds
dt3 = dt1 * 100;  // seconds

function dMsdt = slopes(t, Ms)
  if t <= 5
    K = [0.55, 0.5, 0.4, 0.1, 0.5, 0.55, 0.01]; // 1/s
  else
    K = [0.3, 0.5, 0.4, 0.1, 0.5, 0.3, 0.01];
  end
  dMsdt(1) = K(7) * Ms(2) + K(2) * Ms(4) - K(1) * Ms(1);
  dMsdt(2) = -(K(7) + K(6)) * Ms(2) + K(5) * Ms(3);
  dMsdt(3) = K(3) * Ms(4) + K(6) * Ms(2) - (K(5) + K(4)) * Ms(3);
  dMsdt(4) = K(1) * Ms(1) + K(4) * Ms(3) - (K(2) + K(3)) * Ms(4);
endfunction

function y = ForwardEuler(t, y_prev, h)
  y = y_prev + h * slopes(t, y_prev);
endfunction

function y = HeunsMethod(t, y_prev, h)
  intermediate = slopes(t, y_prev);
  y = y_prev + h / 2 * (intermediate + slopes(t + h, y_prev + h *...
      intermediate));
endfunction

function y = MidpointMethod(t, y_prev, h)
  y = y_prev + h * slopes(t + h / 2, y_prev + h / 2 * slopes(t,...
      y_prev));
endfunction

time1 = [0:dt1:60];
time2 = [0:dt2:60];
time3 = [0:dt3:60];
// [M, AM, AM_p, M_p]
Ms_fe1 = zeros(4, length(time1));
Ms_fe2 = zeros(4, length(time2));
Ms_fe3 = zeros(4, length(time3));
Ms_hm1 = zeros(4, length(time1));
Ms_hm2 = zeros(4, length(time2));
Ms_mm3 = zeros(4, length(time3));
Ms_mm1 = zeros(4, length(time1));
Ms_mm2 = zeros(4, length(time2));
Ms_hm3 = zeros(4, length(time3));
// At t = 0, M = 1, AM = 0, AM_p = 0, M_p = 0
Ms_fe1(:, 1) = [1; 0; 0; 0];
Ms_fe2(:, 1) = [1; 0; 0; 0];
Ms_fe3(:, 1) = [1; 0; 0; 0];
Ms_hm1(:, 1) = [1; 0; 0; 0];
Ms_hm2(:, 1) = [1; 0; 0; 0];
Ms_hm3(:, 1) = [1; 0; 0; 0];
Ms_mm1(:, 1) = [1; 0; 0; 0];
Ms_mm2(:, 1) = [1; 0; 0; 0];
Ms_mm3(:, 1) = [1; 0; 0; 0];

for t1 = 1:length(time1) - 1
  Ms_fe1(:, t1 + 1) = ForwardEuler(time1(t1), Ms_fe1(:, t1), dt1);
  Ms_hm1(:, t1 + 1) = HeunsMethod(time1(t1), Ms_hm1(:, t1), dt1);
  Ms_mm1(:, t1 + 1) = MidpointMethod(time1(t1), Ms_hm1(:, t1), dt1);
  if modulo(t1, 10) == 0
    t2 = t1 / 10;
    Ms_fe2(:, t2 + 1) = ForwardEuler(time2(t2), Ms_fe2(:, t2), dt2);
    Ms_hm2(:, t2 + 1) = HeunsMethod(time2(t2), Ms_hm2(:, t2), dt2);
    Ms_mm2(:, t2 + 1) = MidpointMethod(time2(t2), Ms_mm2(:, t2),...
        dt2);
  end
  if modulo(t1, 100) == 0
    t3 = t1 / 100;
    Ms_fe3(:, t3 + 1) = ForwardEuler(time3(t3), Ms_fe3(:, t3), dt3);
    Ms_hm3(:, t3 + 1) = HeunsMethod(time3(t3), Ms_hm3(:, t3), dt3);
    Ms_mm3(:, t3 + 1) = MidpointMethod(time3(t3), Ms_mm3(:, t3),...
        dt3);
  end
end  // t1
plot(time1, (Ms_fe1(2, :) + Ms_fe1(3, :)));
plot(time2, (Ms_fe2(2, :) + Ms_fe2(3, :)), 'm');
plot(time3, (Ms_fe3(2, :) + Ms_fe3(3, :)), 'g');
plot(time1, (Ms_hm1(2, :) + Ms_hm1(3, :)), ':');
plot(time2, (Ms_hm2(2, :) + Ms_hm2(3, :)), 'm:');
plot(time3, (Ms_hm3(2, :) + Ms_hm3(3, :)), 'g:');
plot(time1, (Ms_mm1(2, :) + Ms_mm1(3, :)), '--');
plot(time2, (Ms_mm2(2, :) + Ms_mm2(3, :)), 'm--');
plot(time3, (Ms_mm3(2, :) + Ms_mm3(3, :)), 'g--');
title("Stress vs time");
xlabel("$Time\ t$", "fontsize", 3);
ylabel("$Stress\ (AM+AM_p)$", "fontsize", 3);
legend(["FE, dt = 0.01"; "FE, dt = 0.1"; "FE, dt = 1";...
    "HM, dt = 0.01"; "HM, dt = 0.1"; "HM, dt = 1";...
    "MM, dt = 0.01"; "MM, dt = 0.1"; "MM, dt = 1"], -1);
