clear;
clf;
// Biomechanical parameters
F = 8500;  // N, force of the tackle
T = 833;  // N, compressive force due to body weight
L = 40;  // cm, length of tibia
I = 26684 * 0.1^4;  // cm^4, moment of inertia
E = 18e9 / 100^2;  // N / cm^2, young's modulus of the tibia
a = 20;  // cm, point on tibia where the force is acting on

// parameters in simplified equation
mu = F / E / I / L;  // 1 / cm^3
alpha = T / E / I;  // 1 / cm^2

// parameters for solver
dx = 0.1;  // cm
len = [0:dx:L];

function f = rhs(x, a)
  if (x > a)
    f = mu * a * (L - x);
  else
    f = mu * x * (L - a);
  end
endfunction

// indicate first run of loop so we don't try to guess for a new
// beta
first_time = %T;
counter = 1;  // toggle to alternate between beta1 and beta2
tol = 1e-6;  // tolerance value

// auxiliary variables, y1 = y, y2 = y' = y1'
y1_s = zeros(len);
y2_s = zeros(len);

// boundary value for y(L) = 0;
bv = 0;
// set boundary condition to be false so while loop can start
y1_s($) = bv + 2 * tol;
// intial guesses for y2
// [beta1, B1;
//  beta2, B2]
shoot = [0, 0;
         1, 0];
// shooting method (secant)
while abs(y1_s($) - bv) > tol
  index = 2 - modulo(counter, 2);
  // initial condition for y(0) = 0
  y1_s(1) = 0;
  y2_s(1) = shoot(index, 1);
  for x = 1:length(len) - 1;
    // Forward Euler
    y1_s(x + 1) = y1_s(x) + dx * y2_s(x);
    y2_s(x + 1) = y2_s(x) + dx * (rhs(len(x), a) - alpha * y1_s(x));
  end  // x
  shoot(index, 2) = y1_s($);
  // Determing new beta from the previous two guesses
  if (~first_time)
    slope = (shoot(1, 2) - shoot(2, 2)) / (shoot(1, 1) -...
        shoot(2, 1));
    shoot(1 + modulo(index, 2), 1) = shoot(index, 1) + (bv -...
        shoot(index, 2)) / slope;
  end
  first_time = %F;
  // alternating between beta1 and beta2
  counter = modulo(counter + 1, 2);
end

// equilibrium method
// new length vector which starts at i = 2 and ends at i = n - 1
len_e = [len(2):dx:len($-1)];
// tri-diagonal matrix
A = zeros(length(len_e), length(len_e));
B = zeros(length(len_e));
// boundary values
bv1 = 0;
bv2 = 0;
y_e = zeros(len_e);

for y = 1:length(len_e)
  // no need if else statements for B since boundary values are zero
  B(y) = dx * dx * rhs(len_e(y), a);
  if (y == 1)
    A(y, y) = -2 + dx * dx * alpha;
    A(y, y + 1) = 1;  // P(x) is zero
  elseif y == length(len_e)
    A(y, y - 1) = 1;
    A(y, y) = -2 + dx * dx * alpha;
  else
    A(y, y - 1) = 1;
    A(y, y) = -2 + dx * dx * alpha;
    A(y, y + 1) = 1;
  end
end

y_e = A \ B;

plot(len, y1_s);
plot(len_e', y_e, 'r-');
