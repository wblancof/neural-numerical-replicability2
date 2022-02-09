#include <cmath>        // pow
#include <array>
#include <vector>
#include <fstream>      // std::ofstream
#include <iostream>     // std::cout
#include <filesystem>

#include "SimulationParameters.h"
#include "SimulationInitialization.h"
#include "iappInit.hpp"
#include "rk4.hpp"

// Create a table of applied currents Ii
// Table iapp % 100 0 99 i0+t*deli/99
std::vector<actualRealType> iapp(nNeurons);                   // for Excitatory neurons
actualRealType iappj;

// Synaptic drive from excitatory cells
actualRealType atotExc, atotExcj;
// Synaptic drive from inhibitory cells
actualRealType atotInh, atotInhj;

// v[i] membrane potential of cell i
// n[i] activation of K+ conductance for cell i
// a[i] synaptic drive from cell i
// s[i] synaptic recovery for terminals from cell i

// Vector of neurons -> The Network
std::array<std::array<actualRealType, 4>, nNeurons> network;

// SpikeTrain Class to manage/store/save the spike times for each neuron
// SpikeTrain *mSpikeTrain = new SpikeTrain(nNeurons);

// Auxiliary variables for output
std::vector<std::array<actualRealType, 4> > ave;               // Average of V N A and S
std::vector<std::array<actualRealType, 4> > networkOutputs;    // Average Ainh and Aexc, episodes time,

//==========
// Functions
//==========

// Rate constants and steady state activation function for Na+ and K+ currents
// am(v)=.1*(25-v)/(exp(.1*(25-v))-1)
inline actualRealType am (const actualRealType v) {
    return dot1*(int25-v)/(exp(dot1*(int25-v))-int1);
}

// bm(v)=4.0*exp(-v/18)
inline actualRealType bm (const actualRealType v) {
    return int4*exp(-v/int18);
}

// minf(v)=am(v)/(am(v)+bm(v))
inline actualRealType minf (const actualRealType v) {
    return am(v)/(am(v)+bm(v));
}

// bn(v)=.125*exp(-v/80)
inline actualRealType bn (const actualRealType v) {
    return dot125*exp(-v/int80);
}

// an(v)=.01*(10-v)/(exp(.1*(10-v))-1.0)
inline actualRealType an (const actualRealType v) {
    return dot01*(int10-v)/(exp(dot1*(int10-v))-int1);
}

// Synaptic activation
// fsyn(V) = 1./(1+exp((Vthresh-V)/kv))
inline actualRealType fsyn (const actualRealType v) {
    //return 1.0/(1+exp((Vthresh-v)/kv));
    return int1/(int1+exp(Vthresh-v)); // if kv = 1;
}

// System ODEs per Neuron
// This function uses the following global variables
//  vna=115.0;
//  vk = -12.0;
//	vl=10.6;
//	gnabar=36.0;
//	gkbar=12.0;
//	gl=0.1;
//	h0=0.8;
//  iappj;
//  gsyn=3.6/nNeurons;
void HH_NeuronModel(actualRealType t, const std::array<actualRealType, 4> network, std::array<actualRealType, 4> &dxdt) {
    actualRealType vj = network[0];
    actualRealType nj = network[1];
    actualRealType aj = network[2];
    actualRealType sj = network[3];
    
    // Differential equations (XPP original version)
    // v[i] membrane potential of cell i
    // v[0..99]'= -gl*(v[j]-vl)
    //              -gnabar*minf(v[j])^3*(h0-n[j])*(v[j]-vna)
    //              -gkbar*n[j]^4*(v[j]-vk)
    //              -gsyn*(atot-a[j]*s[j]/100)*(v[j]-vsyn)
    //              +iapp([j])
    dxdt[0] = -gl * (vj - vl)
                    - gnabar * pow(minf(vj), int3) * (h0 - nj) * (vj - vna)
                    - gkbar * pow(nj, int4) * (vj - vk)
                    - gsyn * atotExcj * (vj - vExc)
                    - gsyn * atotInhj * (vj - vInh)
                    + iappj;
    
    // n[i] Activation of K+ conductance for cell i
    // n[0..99]'= an(v[j])-(an(v[j])+bn(v[j]))*n[j]
    dxdt[1] = an(vj) - (an(vj) + bn(vj)) * nj;

    // a[i] Synaptic drive from cell i
    // a[0..99]'= fsyn(v[j])*(1-a[j])/tauf - a[j]/taus
    dxdt[2] = fsyn(vj) * (int1 - aj) / tauf - aj / taus;

    // s[i] Synaptic recovery for terminals from cell i
    // s[0..99]'=alphad*(1-s[j])-betad*fsyn(v[j])*s[j]
    dxdt[3] = alphad * (int1 - sj) - betad * fsyn(vj) * sj;
}

void initEachNeuron() {


    // Initial conditions
    // init v[0..99]=0 n[j]=0 a[j]=0.01 s[j]=.25
    for (unsigned int n = 0; n < nNeurons; n++) {
    	network[n][0] = zero;        // v
    	network[n][1] = zero;        // n
    	network[n][2] = dot01;       // a
    	network[n][3] = dot25;       // s
        depolarization[n] = false;   // no spike found
    }
}

// Function to write a .txt file with the output parameters of the network state
void writeToFile (std::string const fileNameStr, std::vector<std::array<actualRealType, 4> > v) {
    std::ofstream myfile;
    std::cout << "Writing in file: " << fileNameStr << std::endl;
    try {
        myfile.open(fileNameStr);
        myfile.precision(n_precision);                              // Adjust precision
        myfile.setf(std::ios::fixed, std::ios::floatfield);
        for(unsigned int i = 0; i < v.size(); i++)
            myfile << v[i][0] << '\t' << v[i][1] << '\t' << v[i][2] << '\t' << v[i][3] << std::endl;
        myfile.close();
    } catch(std::exception &e) {
        std::cout << "Unable to open file: " << fileNameStr << std::endl;
    }
}

// Function to write a .txt with the iapp values
void writeIappToFile(std::string const fileNameStr) {
    std::ofstream myfile;
    std::cout << "Writing in file: " << fileNameStr << std::endl;

    try {
        myfile.open(fileNameStr);
        myfile.precision(n_precision);                              // Adjust precision
        myfile.setf(std::ios::fixed, std::ios::floatfield);
        for(unsigned int i = 0; i < nNeurons; i++) {
            myfile << "iapp[" << i << "]=" << iapp[i] << ";" << std::endl;
        }
        myfile.close();
    } catch(std::exception &e) {
        std::cout << "Unable to open file: " << fileNameStr << std::endl;
    }
}

void showParameters() {
    std::cout << "GCC version: " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__ << std::endl;
    std::cout << "Running with the parameters:" << std::endl;
    std::cout << "Total Neurons = "<< nNeurons << std::endl;
    std::cout << "Exc Neurons = "<< nExcNeurons << std::endl;
    std::cout << "Inh Neurons = "<< nInhNeurons << std::endl;
    std::cout << "maxTimeSimulation = " << maxTimeSimulation / 1000 << "s" << std::endl;
    std::cout << "nBurst = " << maxNumBurst << std::endl;
    std::cout << "dt = " << dt << "ms" << std::endl;
    std::cout << "vInh = " << vInh << std::endl;
    std::cout << "Set precision = "<< n_precision << std::endl;
    std::cout << "Save Simulation = " << SAVE_SIMULATION << std::endl;
    std::cout << "Organization = " << iappOrder << std::endl;
    std::cout << std::endl;
}

// Check if the directory exist, case not, it will be created 
// =============================================================
void directory_existsCreate(std::string s){
    if(!std::filesystem::exists(s)) {
        try {
            std::filesystem::create_directory(s);
    	} catch (const std::exception& ex) {
        	std::cout << "Error creating the " << s << " directory" << std::endl;
        	exit(0);
    	}
	}
}

int main(int argc, char* argv[]) {
    actualRealType t = 0;
    actualRealType tspan[2];
    std::array<actualRealType, 4> r; /* Vetor de retorno da funcao rk4, o tamanho depende da quantidade de passos (n + 1) */
    actualRealType actTotalExc;
    actualRealType actTotalInh;
    actualRealType v_1;
    actualRealType sTotalExc;
    actualRealType sTotalInh;
    std::array<actualRealType, 4> sN_1;
    std::array<actualRealType, 4> sN;
    std::array<actualRealType, 4> sA;
    std::vector<actualRealType *> outputAct;
    clock_t simStartTime;

    // Creating output logs directory
    // ============================
    directory_existsCreate(outputDir);

    // Initial conditions applied current for every neuron
    // The same distribution for all simulations
    // =========================================
    iappInit(iapp, iappOrder);

    // Simulation
    // ==================================
    std::cout << "==> Simulation <==" << std::endl;
    showParameters();

    // Initial conditions for every neuron
    // ===================================
    initEachNeuron();

    simStartTime = clock();     // To get the simulation time
    std::cout.precision(4);     // Adjust precision to cout
    std::cout.setf(std::ios::fixed, std::ios::floatfield);

    while(burstCount < maxNumBurst && t <= maxTimeSimulation) {
        // Total synaptic drive exc/inh for all cells
        atotExc = zero;
        atotInh = zero;

        // Synaptic drive from excitatory cells
        // atotexc=sum(0,100)of(shift(s0,i')*shift(a0,i'))/100
        // remembering the order: v, n, a, s -> 0, 1, 2, 3
        for (unsigned n = 0; n < nExcNeurons; n++) {
            atotExc += network[n][3] * network[n][2];
        }
        // Synaptic drive from inhibitory cells
        for (unsigned n = nExcNeurons; n < nNeurons; n++) {
            atotInh += network[n][3] * network[n][2];
        }

        // Saving the previous state of the network to detecting episodes
        sN_1[0] = zero;
        sN_1[1] = zero;
        sN_1[2] = zero;
        sN_1[3] = zero;
        for(unsigned n = 0; n < nNeurons; n++) {
            sN_1[0] += network[n][0];
            sN_1[1] += network[n][1];
            sN_1[2] += network[n][2];
            sN_1[3] += network[n][3];
        }

        sN_1[0] = sN_1[0] / nNeurons;
        sN_1[1] = sN_1[1] / nNeurons;
        sN_1[2] = sN_1[2] / nNeurons;
        sN_1[3] = sN_1[3] / nNeurons;

        // Initialize the Network state
        sN[0] = zero;
        sN[1] = zero;
        sN[2] = zero;
        sN[3] = zero;

        // For each neuron
        for(unsigned n = 0; n < nNeurons; n++) {
            // Synaptic drive exc/inh
            atotExcj = atotExc;
            atotInhj = atotInh;
            // Taking out the contribution of the own cell i
            iappj = iapp[n];                                // Get the original distribution [-10 .. 5] iapp dist
            if (n < nExcNeurons) { // Excitatory
                atotExcj -= network[n][3] * network[n][2];    // remove its own synapse
            } else {             // Inhibitory
                atotInhj -= network[n][3] * network[n][2];    // remove its own synapse
            }

            v_1 = network[n][0];                            // Save previous voltage, used to detect spikes
            tspan[0] = t;
            tspan[1] = t + dt;
            // Integration, solving the ODE
            rk4(HH_NeuronModel, tspan, network[n], 4, r);
            network[n][0] = r[0];
            network[n][1] = r[1];
            network[n][2] = r[2];
            network[n][3] = r[3];

            // ==========================================================================================
            //        Detecting spikes and Depolarization
            // It is not already on Depolarization period and Voltage >= Vth and Positive slope
            if (!depolarization[n] && network[n][0] >= Vthresh && (network[n][0] - v_1) / dt > zero) {
                depolarization[n] = true;
            }
            // ==========================================================================================
            //        Detecting Repolarization
            //  It was on Depolarization period and Voltage <= Vth and Negative slope
            if ( depolarization[n] && network[n][0] <= Vthresh && (network[n][0] - v_1) / dt < zero) {
                depolarization[n] = false;
                // mSpikeTrain->addSpikeTimeToNeuron(n,t);
            }

            // Accumulating/Sum for each neuron state
            sN[0] += network[n][0];
            sN[1] += network[n][1];
            sN[2] += network[n][2];
            sN[3] += network[n][3];
        }
        // Normalize the current state of the network for V A N and S
        sN[0] = sN[0] / nNeurons;
        sN[1] = sN[1] / nNeurons;
        sN[2] = sN[2] / nNeurons;
        sN[3] = sN[3] / nNeurons;

        // Population activity and synaptic recovery for excitatory neurons
        actTotalExc = zero;
        sTotalExc = zero;
        for(unsigned n = 0; n < nExcNeurons; n++) {
            actTotalExc += network[n][2];
            sTotalExc += network[n][3];
        }

        // Population activity and synaptic recovery for inhibitory neurons
        actTotalInh = zero;
        sTotalInh = zero;
        for(unsigned n = nExcNeurons; n < nNeurons; n++) {
            actTotalInh += network[n][2];
            sTotalInh += network[n][3];
        }

        // Detecting episode
        // It is not already on active phase and Activity > thA and Activity derivative > thDA
        if (activePhase == false && sN[2] >= thA && (sN[2] - sN_1[2]) / dt > thDA) {
            activePhase = true;

            // Save values of Aexc, Ainh, time episode, 1
            sA[0] = actTotalExc / nExcNeurons;
            sA[1] = actTotalInh / nInhNeurons;
            sA[2] = t;
            sA[3] = int1;
            networkOutputs.push_back(sA);     // Aexc, Ainh, time episode, 1
        // It was on active phase and Activity < thA
        } else if (activePhase == true && sN[2] < thA) {
            activePhase = false;
            burstCount++;                           // Counting the episode
            std::cout << "<" << nExcNeurons << "," << vInh << ">" << "- burst:" << burstCount << ", time: "<< t << std::endl;

            // Save values of Aexc, Ainh, time episode, 1
            sA[0] = actTotalExc / nExcNeurons;
            sA[1] = actTotalInh / nInhNeurons;
            sA[2] = -t;
            sA[3] = int1;
            networkOutputs.push_back(sA);        // Aexc, Ainh, -time episode, 1
        }

        // Saving the current state of the network
        ave.push_back(sN);

        // increment time
        t += dt;
    }

    std::cout << "Simulation duration: " << ((float)(clock() - simStartTime)) / CLOCKS_PER_SEC <<" seconds"<< std::endl;
    std::cout << "=========================================" << std::endl;

    if (SAVE_SIMULATION) {
        std::string solver = "rk4_";
        std::string sdt = "dt0" + std::to_string((int)(dt*to_actualRealType("10000"))) + "_";

        std::string fileNameStr = outputDir + "/HH_BBT_" +
                solver + sdt +
                std::to_string(nNeurons) + "," + std::to_string(nInhNeurons) + ",vI" +
                std::to_string((int)vInh) + ",t=" +
                p_maxTimeSimulation +"ms_" +
                actualPrecisionType + "_Iapp" + iappOrder;
        writeToFile(fileNameStr + ".txt", ave);

        // Saving the Episodes start/end times
        writeToFile( fileNameStr + ",Epis.txt", networkOutputs);

        // Saving the applied currents values
        writeIappToFile(fileNameStr + ",Iapp.txt");
    }
    return 0;
}
