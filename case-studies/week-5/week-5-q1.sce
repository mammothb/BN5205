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
// indicate first run of loop so we don't try to guess for a new
// beta
first_time = %T;
counter = 1;  // toggle to alternate between beta1 and beta2
tol = 1e-6;  // tolerance value

// auxiliary variables, y1 = y, y2 = y' = y1'
y1 = zeros(len);
y2 = zeros(len);

// boundary value for y(L) = 0;
bv = 0;
// set boundary condition to be false so while loop can start
y1($) = bv + 2 * tol;
// intial guesses for y2
// [beta1, B1;
//  beta2, B2]
shoot = [0, 0;
         1, 0];
// shooting method (secant)
while abs(y1($) - bv) > tol
  index = 2 - modulo(counter, 2);
  // initial condition for y(0) = 0
  y1(1) = 0;
  y2(1) = shoot(index, 1);
  for x = 1:length(len) - 1;
    // Forward Euler
    y1(x + 1) = y1(x) + dx * y2(x);
    if (len(x) > a)
      y2(x + 1) = y2(x) + dx * (mu * a * (L - len(x)) - alpha *...
          y1(x));
    else
      y2(x + 1) = y2(x) + dx * (mu * len(x) * (L - a) - alpha *...
          y1(x));
    end
  end  // x
  shoot(index, 2) = y1($);
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

plot(len, y1)
