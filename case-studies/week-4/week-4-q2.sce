clear;
clf;

// constant variables
dt = 0.1;

function dMsdt = slopes(t, Ms)
  if t <= 5
    K = [0.55, 0.5, 0.4, 25, 0.5, 0.55, 0.01]; // 1/s
  else
    K = [0.3, 0.5, 0.4, 25, 0.5, 0.3, 0.01];
  end
  dMsdt(1) = K(7) * Ms(2) + K(2) * Ms(4) - K(1) * Ms(1);
  dMsdt(2) = -(K(7) + K(6)) * Ms(2) + K(5) * Ms(3);
  dMsdt(3) = K(3) * Ms(4) + K(6) * Ms(2) - (K(5) + K(4)) * Ms(3);
  dMsdt(4) = K(1) * Ms(1) + K(4) * Ms(3) - (K(2) + K(3)) * Ms(4);
endfunction

function y = MidpointMethod(t, y_prev, h)
  y = y_prev + h * slopes(t + h / 2, y_prev + h / 2 * slopes(t,...
      y_prev));
endfunction

time = [0:dt:60];
Ms = zeros(4, length(time));
Ms_mm = zeros(4, length(time));
// At t = 0, M = 1
Ms(1, 1) = 1;
Ms_mm(1, 1) = 1;

for t = 1:length(time) - 1
  if t <= 5
    K = [0.55 0.5 0.4 25 0.5 0.55 0.01]; // 1/s
  else
    K = [0.3 0.5 0.4 25 0.5 0.3 0.01];
  end
  LHS = [1 + dt * K(1), -dt * K(7), 0, -dt * K(2);
         0, 1 + dt * (K(7) + K(6)), -dt * K(5), 0;
         0, -dt * K(6), 1 + dt * (K(5) + K(4)), -dt * K(3);
         -dt * K(1), 0, -dt * K(4), 1 + dt * (K(2) + K(3))];
  Ms(:, t + 1) = LHS \ Ms(:, t);
  Ms_mm(:, t + 1) = MidpointMethod(time(t), Ms_mm(:, t), dt);
end  // t

plot(time, Ms(2, :) + Ms(3, :));
//plot(time, Ms_mm(2, :) + Ms_mm(3, :), 'r');
title("Stress vs time (with muscle relaxant)");
xlabel("$Time\ t$", "fontsize", 3);
ylabel("$Stress\ (AM+AM_p)$", "fontsize", 3);
legend(["BE"; "MM"], -1);
