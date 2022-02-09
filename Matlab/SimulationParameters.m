








% Input/Output parameters
% RAND, ASC, DES
global iappOrder; iappOrder = "DES";
global SAVE_SIMULATION; SAVE_SIMULATION = true;
global n_precision; n_precision = 14;
global outputDir; outputDir = './results';

% Simulation constants
global p_dt; p_dt = 0.01;                                  % Time step: ms
global p_maxTimeSimulation; p_maxTimeSimulation = 8000.01;    % Maximum simulation time: ms

% Simulation parameters
global nNeurons; nNeurons = 100;
global nExcNeurons; nExcNeurons = 80;                      % Amount of excitatory neurons
global maxNumBurst; maxNumBurst = 200;                     % Amount of burst to be reached to stop the simulation

% Parameters needed to detect episodes
global p_thA; p_thA = 0.1730;                              % thA = 0.25*(maxAt-minAt); based on Patrick paper. Previous calculated in Matlab
global p_thDA; p_thDA = 0.1490;

% Cellular parameters
% p deli=15
% p vna=115  vk=-12  vl=10.6  gnabar=36  gkbar=12  gl=0.1
% p h0=0.8
global p_deli; p_deli = 15.0;
global p_vna; p_vna = 115.0;
global p_vk; p_vk = -12.0;
global p_vl; p_vl = 10.6;
global p_gnabar; p_gnabar = 36.0;
global p_gkbar; p_gkbar = 12.0;
global p_gl; p_gl = 0.1;
global p_h0; p_h0 = 0.8;

% p vexc=70
% p vinh = -12
global p_vExc; p_vExc = 70.0;
global p_vInh; p_vInh = 70.0;
global p_vsyn; p_vsyn = 70.0;

% Synaptic parameters
% p taus=10 tauf=1
% p gsyn=3.6 vsyn=70
% p alphad=0.0015 betad=0.12
% p Vthresh=40 kv=1
global p_taus; p_taus = 10.0;
global p_tauf; p_tauf = 1.0;
global p_gsyn; p_gsyn = 3.6;
global p_alphad; p_alphad = 0.0015;
global p_betad; p_betad = 0.12;
global p_Vthresh; p_Vthresh = 40.0;
global p_kv; p_kv = 1.0;