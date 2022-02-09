#ifndef SIMULATION_INITIALIZATION_H_
#define SIMULATION_INITIALIZATION_H_

#include "SimulationParameters.h"
#include <string>

#ifdef HighPrecision
    #include <boost/multiprecision/cpp_dec_float.hpp>
    typedef boost::multiprecision::cpp_dec_float_100 actualRealType;

    inline actualRealType to_actualRealType(const char *str) {
        return actualRealType(str);
    }
    std::string const actualPrecisionType = "HighPrecision";


#else
    typedef double actualRealType;

    inline actualRealType to_actualRealType(const char *str) {
        return atof(str);
    }
    std::string const actualPrecisionType = "DoublePrecision";
#endif





// Simulation constants
const actualRealType dt = to_actualRealType(p_dt);
const actualRealType maxTimeSimulation = to_actualRealType(p_maxTimeSimulation);

// Simulation parameters
const unsigned int nInhNeurons = nNeurons - nExcNeurons;                            // Amount of inhibitory neurons

// Parameters needed to detect episodes
const actualRealType thA = to_actualRealType(p_thA);
const actualRealType thDA = to_actualRealType(p_thDA);
bool activePhase = false;
unsigned int burstCount = 0;

// Parameters needed to detect spikes
bool depolarization[nNeurons];

// Cellular parameters
const actualRealType deli = to_actualRealType(p_deli);
const actualRealType vna = to_actualRealType(p_vna);
const actualRealType vk = to_actualRealType(p_vk);
const actualRealType vl = to_actualRealType(p_vl);
const actualRealType gnabar = to_actualRealType(p_gnabar);
const actualRealType gkbar = to_actualRealType(p_gkbar);
const actualRealType gl = to_actualRealType(p_gl);
const actualRealType h0 = to_actualRealType(p_h0);
const actualRealType vExc = to_actualRealType(p_vExc);
const actualRealType vInh = to_actualRealType(p_vInh);
const actualRealType vsyn = to_actualRealType(p_vsyn);

// Synaptic parameters
const actualRealType taus = to_actualRealType(p_taus);
const actualRealType tauf = to_actualRealType(p_tauf);
const actualRealType gsyn = to_actualRealType(p_gsyn) / nNeurons;
const actualRealType alphad = to_actualRealType(p_alphad);
const actualRealType betad = to_actualRealType(p_betad);
const actualRealType Vthresh = to_actualRealType(p_Vthresh);
const actualRealType kv = to_actualRealType(p_kv);

// Constants inside functions
const actualRealType dot1 = to_actualRealType("0.1");
const actualRealType dot01 = to_actualRealType("0.01");
const actualRealType dot25 = to_actualRealType("0.25");
const actualRealType dot125 = to_actualRealType("0.125");
const actualRealType zero = to_actualRealType("0");
const actualRealType int1 = to_actualRealType("1");
const actualRealType int2 = to_actualRealType("2");
const actualRealType int3 = to_actualRealType("3");
const actualRealType int4 = to_actualRealType("4");
const actualRealType int6 = to_actualRealType("6");
const actualRealType int10 = to_actualRealType("10");
const actualRealType int18 = to_actualRealType("18");
const actualRealType int25 = to_actualRealType("25");
const actualRealType int80 = to_actualRealType("80");
#endif