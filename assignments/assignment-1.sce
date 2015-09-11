clear;
clf;
//==============================================================================
// ms, time constant for presistent sodium channel
// \param voltage membrane voltage
// \return time_constant time constant for persistent sodium channel at said
//         membrane potential
//==============================================================================
function time_constant = tau_NaP(voltage)
  time_constant = 1000 + 10000 / (1 + exp((voltage + 60) / 10));
endfunction

//==============================================================================
// (unitless), steady state of 'n' gating particle
// \param voltage membrane voltage
// \return state steady state of 'n' gating particle at said membrane voltage
//==============================================================================
function state = n_ss(voltage)
  state = 1 / (1 + exp(-(voltage + 43) / 3.9));
endfunction

//==============================================================================
// (unitless), steady state of transient 'm' gating particle
// \param voltage membrane voltage
// \return state steady state of transient 'm' gating particle at said membrane
//         voltage
//==============================================================================
function state = m_t_ss(voltage)
  state = 1 / (1 + exp(-(voltage + 31.3) / 4.3));
endfunction

//==============================================================================
// (unitless), steady state of transient 'h' gating particle
// \param voltage membrane voltage
// \return state steady state of transient 'h' gating particle at said membrane
//         voltage
//==============================================================================
function state = h_t_ss(voltage)
  state = 1 / (1 + exp((voltage + 55) / 7.1));
endfunction

//==============================================================================
// (unitless), steady state of persistent 'm' gating particle
// \param voltage membrane voltage
// \return state steady state of persistent 'm' gating particle at said membrane
//         voltage
//==============================================================================
function state = m_p_ss(voltage)
  state = 1 / (1 + exp(-(voltage + 50) / 6.4));
endfunction

//==============================================================================
// (unitless), steady state of persistent 'h' gating particle
// \param voltage membrane voltage
// \return state steady state of persistent 'h' gating particle at said membrane
//         voltage
//==============================================================================
function state = h_p_ss(voltage)
  state = 1 / (1 + exp((voltage + 52) / 14));
endfunction

//==============================================================================
// Calculates the d()/dt of the necessary variables, V, n, h_t, h_p and noise at
// a particular time point
// \param vals a column of the soln matrix at the specified time point
// \param I stimulus current at the specified time point
// \param xi random variable at specified time point for the Ornstein-Uhlenbeck
//        process, the variable has normal distribution, zero mean and unitary
//        variance generated at the beginning of the simulation
// \return a vector containing [dV/dt; dn/dt; dh_t/dt; dh_p/dt; dnoise/dt]
//==============================================================================
function output = slopes(vals, I, rand_var)
  // reassigning variables for clarity
  V = vals(1);
  n = vals(2);
  h_t = vals(3);
  h_p = vals(4);
  noise = vals(5);
  // values of leak, potassium, transient sodium and persistent sodium currents
  I_Leak = g_Leak * (V - E_Leak);
  I_K = g_bar_K * n * (V - E_K);
  I_NaT = g_bar_NaT * m_t_ss(V) * h_t * (V - E_Na);
  I_NaP = g_bar_NaP * m_p_ss(V) * h_p * (V - E_Na);
  // calculating the d()/dt values
  output(1) = (I - I_Leak - I_K - I_NaT - I_NaP + noise) / C;
  output(2) = (n_ss(V) - n) / tau_K;
  output(3) = (h_t_ss(V) - h_t) / tau_NaT;
  output(4) = (h_p_ss(V) - h_p) / tau_NaP(V);
  output(5) = -noise + 3 * rand_var;
endfunction

// model parameters
E_Leak = -60;  // mV, leak reversal potential
E_Na = 55;  // mV, sodium reversal potential
E_K = -92;  // mV, potassium reversal potential
g_Leak = 2;  // uS/cm^2, leak conductance per unit area
g_bar_NaT = 12;  // uS/cm^2, max transient sodium conductance per unit area
g_bar_NaP = 1.1;  // uS/cm^2, max persistent sodium conductance per unit area
g_bar_K = 6;  // uS/cm^2, maximum potassium conductance per unit area
C = 1;  // nF/cm^2, membrance capacitance per unit area
tau_K = 4;  // ms, time constant for potassium channel
tau_NaT = 3;  // ms, time constant for transient sodium channel

// simulation parameters
dt = 0.1;  // ms, time step
stimulus = 50;  // nA, stimulus
time = [0:dt:400];  // vector containing all time points
I_stim = zeros(length(time), 1);  // stimulus at each time point
// random variable for Ornstein-Uhlenbeck process, both solutions share the same
// set of random variables so that results are comparable
xi = zeros(length(time), 1);
// solution matrix containing only components that requires solving numerically
// [V, n, h_t, h_p, noise]
soln_fe = zeros(5, length(time));  // solution using forward euler
soln_hm = zeros(5, length(time));  // solution using heun's method

// Initial conditions
soln(1, 1) = -60;  // mV, initial membrane potential

rand('normal');
for t = 1:length(time)
  // setting up random variable
  xi(t) = rand();
  // setting up stimulus current from 50ms to 250ms
  if time(t) >= 50 & time(t) <= 250 then
    I_stim(t) = stimulus;
  end
end

// solving using Forward Euler and Heun's Method (RK2)
for t = 1:length(time) - 1
  // forward euler
  soln_fe(:, t + 1) = soln_fe(:, t) + dt * slopes(soln_fe(:, t), I_stim(t),...
      xi(t));
  // heun's method
  intermediate = slopes(soln_hm(:, t), I_stim(t), xi(t));
  soln_hm(:, t + 1) = soln_hm(:, t) + dt / 2 * (intermediate +...
      slopes(soln_hm(:, t) + dt * intermediate, I_stim(t), xi(t)));
end

plot(time, soln_fe(1, :), 'b-');
plot(time, soln_hm(1, :), 'r-');
title("$Membrane\ potential\ V\ vs\ time\ t$", "fontsize", 4);
ylabel("$Membrane\ potential\ V\ (mV)$", "fontsize", 4);
xlabel("$Time\ t\ (ms)$", "fontsize", 4);
legend(['Forward Euler'; 'Heun''s Method (RK2)'], -1);
