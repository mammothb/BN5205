clear;
clf;

// model parameters
E_Leak = -60;  // mV, leak reversal potential
E_Na = 55;  // mV, sodium reversal potential
E_K = -92;  // mV, potassium reversal potential
g_Leak = 2;  // uS/cm^2, leak conductance per unit area
g_bar_NaT = 12;  // uS/cm^2, max transient sodium conductance per unit area
g_bar_NaP = 1.1;  // uS/cm^2, max persistent sodium conductance per unit area
g_bar_K = 6;  // uS/cm^2, maximum potassium conductance per unit area
C = 1;  // nC/cm^2, membrance capacitance per unit area
tau_K = 4;  // ms, time constant for potassium channel
tau_NaT = 3;  // ms, time constant for transient sodium channel

// ms, time constant for presistent sodium channel
// \param voltage membrane voltage
// \return time_constant time constant for persistent sodium channel at said
//         membrane potential
function time_constant = tau_NaP(voltage)
  time_constant = 1000 + 10000 / (1 + exp((voltage + 60) / 10));
endfunction

// (unitless), steady state of 'n' gating particle
// \param voltage membrane voltage
// \return state steady state of 'n' gating particle at said membrane voltage
function state = n_ss(voltage)
  state = 1 / (1 + exp(-(voltage + 43) / 3.9));
endfunction

// (unitless), steady state of transient 'm' gating particle
// \param voltage membrane voltage
// \return state steady state of transient 'm' gating particle at said membrane
//         voltage
function state = m_t_ss(voltage)
  state = 1 / (1 + exp(-(V + 32) / 4.3));
endfunction

// (unitless), steady state of transient 'h' gating particle
// \param voltage membrane voltage
// \return state steady state of transient 'h' gating particle at said membrane
//         voltage
function state = h_t_ss(voltage)
  state = 1 / (1 + exp((voltage + 55) / 7.1));
endfunction

// (unitless), steady state of persistent 'm' gating particle
// \param voltage membrane voltage
// \return state steady state of persistent 'm' gating particle at said membrane
//         voltage
function state = m_p_ss(voltage)
  state = 1 / (1 + exp(-(voltage + 50) / 6.4));
endfunction

// (unitless), steady state of persistent 'h' gating particle
// \param voltage membrane voltage
// \return state steady state of persistent 'h' gating particle at said membrane
//         voltage
function state = h_p_ss(voltage)
  state = 1 / (1 + exp((voltage + 52) / 14));
endfunction
