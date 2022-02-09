




























% Simulation constants
global dt; dt = p_dt;
global maxTimeSimulation; maxTimeSimulation = p_maxTimeSimulation;

% Simulation parameters
global nInhNeurons; nInhNeurons = nNeurons - nExcNeurons;                            % Amount of inhibitory neurons

% Parameters needed to detect episodes
global thA; thA = p_thA;
global thDA; thDA = p_thDA;
global activePhase; activePhase = false;
global burstCount; burstCount = 0;

% Parameters needed to detect spikes
global depolarization; depolarization = [];

% Cellular parameters
global deli; deli = p_deli;
global vna; vna = p_vna;
global vk; vk = p_vk;
global vl; vl = p_vl;
global gnabar; gnabar = p_gnabar;
global gkbar; gkbar = p_gkbar;
global gl; gl = p_gl;
global h0; h0 = p_h0;
global vExc; vExc = p_vExc;
global vInh; vInh = p_vInh;
global vsyn; vsyn = p_vsyn;

% Synaptic parameters
global taus; taus = p_taus;
global tauf; tauf = p_tauf;
global gsyn; gsyn = p_gsyn / nNeurons;
global alphad; alphad = p_alphad;
global betad; betad = p_betad;
global Vthresh; Vthresh = p_Vthresh;
global kv; kv = p_kv;