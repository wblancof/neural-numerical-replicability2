#ifndef SIMULATION_PARAMETERS_H_
#define SIMULATION_PARAMETERS_H_

/* Real number type representation */
/* Leave ONLY ONE of the lines bellow uncommented */
// #define DoublePrecision
#define HighPrecision


// Input/Output parameters
// RAND, ASC, DES
const char iappOrder[] = "DES";
const bool SAVE_SIMULATION = true;
const int n_precision = 14;
const std::string outputDir  = "./results";

// Simulation constants
const char p_dt[] = "0.01";                         // Time step: ms
const char p_maxTimeSimulation[] = "8000";          // Maximum simulation time: ms

// Simulation parameters
const unsigned int nNeurons = 100;
const unsigned int nExcNeurons = 80;                // Amount of excitatory neurons
const unsigned int maxNumBurst = 200;               // Amount of burst to be reached to stop the simulation

// Parameters needed to detect episodes
const char p_thA[] = "0.1730";                      // thA = 0.25*(maxAt-minAt); based on Patrick paper. Previous calculated in Matlab
const char p_thDA[] = "0.1490";

// Cellular parameters
// p deli=15
// p vna=115  vk=-12  vl=10.6  gnabar=36  gkbar=12  gl=0.1
// p h0=0.8
const char p_deli[] = "15.0";
const char p_vna[] = "115.0";
const char p_vk [] = "-12.0";
const char p_vl[] = "10.6";
const char p_gnabar[] = "36.0";
const char p_gkbar[] = "12.0";
const char p_gl[] = "0.1";
const char p_h0[] = "0.8";

// p vexc=70
// p vinh = -12
const char p_vExc[] = "70.0";
const char p_vInh[] = "70.0";
const char p_vsyn[] = "70.0";

// Synaptic parameters
// p taus=10 tauf=1
// p gsyn=3.6 vsyn=70
// p alphad=0.0015 betad=0.12
// p Vthresh=40 kv=1
const char p_taus[] = "10.0";
const char p_tauf[] = "1.0";
const char p_gsyn[] = "3.6";
const char p_alphad[] = "0.0015";
const char p_betad[] = "0.12";
const char p_Vthresh[] = "40.0";
const char p_kv[] = "1.0";

#endif