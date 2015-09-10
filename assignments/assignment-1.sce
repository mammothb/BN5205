clear;
clf;
//
// model parameters
E_Leak = -60e-3;  // V, leak reversal potential
E_Na = 55e-3;  // V, sodium reversal potential
E_K = -92e-3;  // V, potassium reversal potential
g_Leak = 2e-6;  // S/cm^2, leak conductance per unit area
g_bar_NaT = 12e-6;  // S/cm^2, max transient sodium conductance per unit area
g_bar_NaP = 1.1e-6;  // S/cm^2, max persistent sodium conductance per unit area
g_bar_K = 6e-6;  // S/cm^2, maximum potassium conductance per unit area
C = 1e-9;  // C/cm^2, membrance capacitance per unit area
tau_K = 4e-3;  // s, time constant for potassium channel
tau_NaT = 3e-3;  // s, time constant for transient sodium channel

// simulation parameters
dt = 0.1e-3;  // ms, time step
stimulus = 50e-9;  // A, stimulus
time = [0:dt:400e-3];  // vector containing all time points
I_stim = zeros(length(time), 1);  // stimulus at each time point
I_noise = zeros(length(time), 1);  // noise current at each time point
// solution matrix containing only components which are not known beforehand
// [V, I_Leak, I_K, I_NaT, I_NaP, n, h_t, h_p]
soln = zeros(8, length(time));

///
// ms, time constant for presistent sodium channel
// \param voltage membrane voltage
// \return time_constant time constant for persistent sodium channel at said
//         membrane potential
//
function time_constant = tau_NaP(voltage)
  time_constant = 1000 + 10000 / (1 + exp((voltage * 1000 + 60) / 10));
endfunction

///
// (unitless), steady state of 'n' gating particle
// \param voltage membrane voltage
// \return state steady state of 'n' gating particle at said membrane voltage
//
function state = n_ss(voltage)
  state = 1 / (1 + exp(-(voltage * 1000+ 43) / 3.9));
endfunction

///
// (unitless), steady state of transient 'm' gating particle
// \param voltage membrane voltage
// \return state steady state of transient 'm' gating particle at said membrane
//         voltage
//
function state = m_t_ss(voltage)
  state = 1 / (1 + exp(-(voltage * 1000+ 31.3) / 4.3));
endfunction

///
// (unitless), steady state of transient 'h' gating particle
// \param voltage membrane voltage
// \return state steady state of transient 'h' gating particle at said membrane
//         voltage
//
function state = h_t_ss(voltage)
  state = 1 / (1 + exp((voltage * 1000+ 55) / 7.1));
endfunction

///
// (unitless), steady state of persistent 'm' gating particle
// \param voltage membrane voltage
// \return state steady state of persistent 'm' gating particle at said membrane
//         voltage
//
function state = m_p_ss(voltage)
  state = 1 / (1 + exp(-(voltage * 1000 + 50) / 6.4));
endfunction

///
// (unitless), steady state of persistent 'h' gating particle
// \param voltage membrane voltage
// \return state steady state of persistent 'h' gating particle at said membrane
//         voltage
//
function state = h_p_ss(voltage)
  state = 1 / (1 + exp((voltage * 1000 + 52) / 14));
endfunction

function output = slopes(vals, I, noise)
  // reassigning variables for clarity
  V = vals(1);
  I_Leak = vals(2);
  I_K = vals(3);
  I_NaT = vals(4);
  I_NaP = vals(5);
  n = vals(6);
  h_t = vals(7);
  h_p = vals(8);
  output(1) = (I - I_Leak - I_K - I_NaT - I_NaP + noise) / C;
  output(2) = g_Leak * (V - E_Leak);
  output(3) = g_bar_K * n * (V - E_K);
  output(4) = g_bar_NaT * m_t_ss(V) * h_t * (V - E_Na);
  output(5) = g_bar_NaP * m_p_ss(V) * h_p * (V - E_Na);
  output(6) = (n_ss(V) - n) / tau_K;
  output(7) = (h_t_ss(V) - h_t) / tau_NaT;
  output(8) = (h_p_ss(V) - h_p) / tau_NaP(V);
endfunction

// Initial conditions
soln(1, 1) = -60e-3;  // V

for t = 1:length(time) - 1
  if time(t) >= 50 & time(t) <= 250 then
    I(t) = stimulus;  // nA, stimulus current
  end
end

for t = 1:length(time) - 1
  soln(:, t + 1) = soln(:, t) + dt * slopes(soln(:, t), I_stim(t), I_noise(t));
end

plot(time, soln(1, :), '-');
