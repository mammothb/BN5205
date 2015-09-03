clear;
clf;
// Biomechanical parameters
F = 8500;  // N, force of the tackle
T = 833;  // N, compressive force due to body weight
L = 40;  // cm, length of tibia
I = 26684 * 0.01^4;  // cm^4, moment of inertia
E = 18e9 / 100^2;  // N / cm^2, young's modulus of the tibia

// parameters in simplified equation
mu = F / E / I / L;  // 1 / cm^3
alpha = T / E / I;  // 1 / cm^2
