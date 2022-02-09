import os
import sys
import time





from SimulationInitialization import *
from iappInit import iappInit
from rk4 import rk4

# Create a table of applied currents Ii
# Table iapp % 100 0 99 i0+t*deli/99
iapp = [zero] * nNeurons                            # for Excitatory neurons
iappj = []

# Synaptic drive from excitatory cells
# atotExc, atotExcj;
# Synaptic drive from inhibitory cells
# atotInh, atotInhj;

# v[i] membrane potential of cell i
# n[i] activation of K+ conductance for cell i
# a[i] synaptic drive from cell i
# s[i] synaptic recovery for terminals from cell i

# Vector of neurons -> The Network
network = [[zero] * 4 for _ in range(nNeurons)]

# SpikeTrain Class to manage/store/save the spike times for each neuron
# mSpikeTrain = []

# Auxiliary variables for output
ave = []               # Average of V N A and S
networkOutputs = []    # Average Ainh and Aexc, episodes time,

#==========
# Functions
#==========

# Rate constants and steady state activation function for Na+ and K+ currents
# am(v)=.1*(25-v)/(exp(.1*(25-v))-1)
def am (v):
    return dot1*(int25-v)/(expFun(dot1*(int25-v))-int1)


# bm(v)=4.0*exp(-v/18)
def bm (v):
    return int4*expFun(-v/int18)


# minf(v)=am(v)/(am(v)+bm(v))
def minf (v):
    return am(v)/(am(v)+bm(v))


# bn(v)=.125*exp(-v/80)
def bn (v):
    return dot125*expFun(-v/int80)


# an(v)=.01*(10-v)/(exp(.1*(10-v))-1.0)
def an (v):
    return dot01*(int10-v)/(expFun(dot1*(int10-v))-int1)


# Synaptic activation
# fsyn(V) = 1./(1+exp((Vthresh-V)/kv))
def fsyn (v):
    #return 1.0/(1+exp((Vthresh-v)/kv));
    return int1/(int1+expFun(Vthresh-v)) # if kv = 1


# System ODEs per Neuron
# This function uses the following global variables
#  vna=115.0;
#  vk = -12.0;
#  vl=10.6;
#  gnabar=36.0;
#  gkbar=12.0;
#  gl=0.1;
#  h0=0.8;
#  iappj;
#  gsyn=3.6/nNeurons;
def HH_NeuronModel(t, network, dxdt):
    vj = network[0]
    nj = network[1]
    aj = network[2]
    sj = network[3]

    # Differential equations (XPP original version)
    # v[i] membrane potential of cell i
    # v[0..99]'= -gl*(v[j]-vl)
    #              -gnabar*minf(v[j])^3*(h0-n[j])*(v[j]-vna)
    #              -gkbar*n[j]^4*(v[j]-vk)
    #              -gsyn*(atot-a[j]*s[j]/100)*(v[j]-vsyn)
    #              +iapp([j])
    dxdt[0] = -gl * (vj - vl) \
                    - gnabar * (minf(vj) ** int3) * (h0 - nj) * (vj - vna) \
                    - gkbar * (nj ** int4) * (vj - vk) \
                    - gsyn * atotExcj * (vj - vExc) \
                    - gsyn * atotInhj * (vj - vInh) \
                    + iappj

    # n[i] Activation of K+ conductance for cell i
    # n[0..99]'= an(v[j])-(an(v[j])+bn(v[j]))*n[j]
    dxdt[1] = an(vj) - (an(vj) + bn(vj)) * nj

    # a[i] Synaptic drive from cell i
    # a[0..99]'= fsyn(v[j])*(1-a[j])/tauf - a[j]/taus
    dxdt[2] = fsyn(vj) * (int1 - aj) / tauf - aj / taus
    
    # s[i] Synaptic recovery for terminals from cell i
    # s[0..99]'=alphad*(1-s[j])-betad*fsyn(v[j])*s[j]
    dxdt[3] = alphad * (int1 - sj) - betad * fsyn(vj) * sj


def initEachNeuron():
    global network
    global depolarization
    # Initial conditions
    # init v[0..99]=0 n[j]=0 a[j]=0.01 s[j]=.25
    for n in range (nNeurons):
        network[n][0] = zero;         # v
        network[n][1] = zero;         # n
        network[n][2] = dot01;        # a
        network[n][3] = dot25;        # s
        depolarization[n] = False     # no spike found



# Function to write a .txt file with the output parameters of the network state
def writeToFile(fileNameStr, v):

    print("Writing in file: {}".format(fileNameStr))
    try:
        myfile = open(fileNameStr, "w")


        for i in range(len(v)):
            myfile.write("{:.{prec}f}\t{:.{prec}f}\t{:.{prec}f}\t{:.{prec}f}\n".format(v[i][0], v[i][1], v[i][2], v[i][3], prec=n_precision))
        myfile.close()
    except:
        print("Unable to open file: {}".format(fileNameStr))



# Function to write a .txt
def writeIappToFile(fileNameStr):

    print("Writing in file: {}".format(fileNameStr))
    try:
        myfile = open(fileNameStr, "w")


        for i in range(nNeurons):
            myfile.write("iapp[{}]={:.{prec}f};\n".format(i, iapp[i], prec=n_precision))
        myfile.close()
    except:
        print("Unable to open file: {}".format(fileNameStr))





def showParameters():
    print("Python version: {}".format(sys.version.replace('\n', '')))
    print("Running with the parameters:")
    print("Total Neurons = {}".format(nNeurons))
    print("Exc Neurons = {}".format(nExcNeurons))
    print("Inh Neurons = {}".format(nInhNeurons))
    print("maxTimeSimulation = {}s".format(maxTimeSimulation/1000))
    print("nBurst = {}".format(maxNumBurst))
    print("dt = {}ms".format(dt))
    print("vInh = {}".format(vInh))
    print("Set precision = {}".format(n_precision))
    print("Save Simulation = {}".format(SAVE_SIMULATION))
    print("Organization = {}".format(iappOrder))
    print("")


# Check if the directory exist, case not, it will be created 
# =============================================================
def directory_existsCreate(s):
    if not os.path.exists(s):
        try:
            os.makedirs(s)
        except:
            print("Error creating the {} directory".format(s))
            exit(0)




if __name__ == '__main__':
    t = zero
    tspan = [zero, zero]
    r = [zero, zero, zero, zero]





    sN_1 = [zero, zero, zero, zero]
    sN = [zero, zero, zero, zero]
    sA = [zero, zero, zero, zero]



    # Creating output logs directory
    # ============================
    directory_existsCreate(outputDir);

    # Initial conditions applied current for every neuron
    # The same distribution for all simulations
    # =========================================
    iappInit(iapp, iappOrder);

    # Simulation
    # ==================================
    print("==> Simulation <==")
    showParameters()

    # Initial conditions for every neuron
    # ===================================
    initEachNeuron()

    simStartTime = time.time()        # To get the simulation time



    while(burstCount < maxNumBurst and t <= maxTimeSimulation):
        # Total synaptic drive exc/inh for all cells
        atotExc = zero
        atotInh = zero

        # Synaptic drive from excitatory cells
        # atotexc=sum(0,100)of(shift(s0,i')*shift(a0,i'))/100
        # remembering the order: v, n, a, s -> 0, 1, 2, 3
        for n in range(nExcNeurons):
            atotExc += network[n][3] * network[n][2]

        # Synaptic drive from inhibitory cells
        for n in range(nExcNeurons, nNeurons):
            atotInh += network[n][3] * network[n][2]


        # Saving the previous state of the network to detecting episodes
        sN_1[0] = zero
        sN_1[1] = zero
        sN_1[2] = zero
        sN_1[3] = zero
        for n in range(nNeurons):
            sN_1[0] += network[n][0]
            sN_1[1] += network[n][1]
            sN_1[2] += network[n][2]
            sN_1[3] += network[n][3]


        sN_1[0] = sN_1[0] / nNeurons
        sN_1[1] = sN_1[1] / nNeurons
        sN_1[2] = sN_1[2] / nNeurons
        sN_1[3] = sN_1[3] / nNeurons

        # Initialize the Network state
        sN[0] = zero
        sN[1] = zero
        sN[2] = zero
        sN[3] = zero

        # For each neuron
        for n in range(nNeurons):
            # Synaptic drive exc/inh
            atotExcj = atotExc
            atotInhj = atotInh
            # Taking out the contribution of the own cell i
            iappj = iapp[n]                               # Get the original distribution [-10 .. 5] iapp dist
            if (n < nExcNeurons):  # Excitatory
                atotExcj -= network[n][3] * network[n][2]   # remove its own synapse
            else:                # Inhibitory
                atotInhj -= network[n][3] * network[n][2]   # remove its own synapse


            v_1 = network[n][0]                           # Save previous voltage, used to detect spikes
            tspan[0] = t
            tspan[1] = t + dt
            # Integration, solving the ODE
            rk4(HH_NeuronModel, tspan, network[n], 4, r)
            network[n][0] = r[0]
            network[n][1] = r[1]
            network[n][2] = r[2]
            network[n][3] = r[3]

            # ==========================================================================================
            #        Detecting spikes and Depolarization
            # It is not already on Depolarization period and Voltage >= Vth and Positive slope
            if not depolarization[n] and network[n][0] >= Vthresh and (network[n][0] - v_1) / dt > zero:
                depolarization[n] = True

            # ==========================================================================================
            #        Detecting Repolarization
            #  It was on Depolarization period and Voltage <= Vth and Negative slope
            if depolarization[n] and network[n][0] <= Vthresh and (network[n][0] - v_1) / dt < zero:
                depolarization[n] = False
                # mSpikeTrain = addSpikeTimeToNeuron(n,t)


            # Accumulating/Sum for each neuron state
            sN[0] += network[n][0]
            sN[1] += network[n][1]
            sN[2] += network[n][2]
            sN[3] += network[n][3]

        # Normalize the current state of the network for V A N and S
        sN[0] = sN[0] / nNeurons
        sN[1] = sN[1] / nNeurons
        sN[2] = sN[2] / nNeurons
        sN[3] = sN[3] / nNeurons

        # Population activity and synaptic recovery for excitatory neurons
        actTotalExc = zero
        sTotalExc = zero
        for n in range(nExcNeurons):
            actTotalExc += network[n][2]
            sTotalExc += network[n][3]


        # Population activity and synaptic recovery for inhibitory neurons
        actTotalInh = zero
        sTotalInh = zero
        for n in range(nExcNeurons, nNeurons):
            actTotalInh += network[n][2]
            sTotalInh += network[n][3]


        # Detecting episode
        # It is not already on active phase and Activity > thA and Activity derivative > thDA
        if (activePhase is False and sN[2] >= thA and (sN[2] - sN_1[2]) / dt > thDA):
            activePhase = True

            # Save values of Aexc, Ainh, time episode, 1
            sA[0] = actTotalExc / nExcNeurons
            sA[1] = division(actTotalInh, nInhNeurons) # sA[1] = actTotalInh / nInhNeurons if nInhNeurons > 0 else zero
            sA[2] = t
            sA[3] = int1
            networkOutputs.append(sA[:]);        # Aexc, Ainh, time episode, 1
        # It was on active phase and Activity < thA
        elif (activePhase is True and sN[2] < thA):
            activePhase = False
            burstCount+=1                     # Counting the episode
            print("<{},{:.4f}>-burst:{}, time: {:.4f}".format(nExcNeurons, vInh, burstCount, t))

            # Save values of Aexc, Ainh, time episode, 1
            sA[0] = actTotalExc / nExcNeurons
            sA[1] = division(actTotalInh, nInhNeurons) # sA[1] = actTotalInh / nInhNeurons if nInhNeurons > 0 else zero
            sA[2] = -t
            sA[3] = int1
            networkOutputs.append(sA[:])         # Aexc, Ainh, -time episode, 1


        # Saving the current state of the network
        ave.append(sN[:])

        # Increment time
        t += dt


    print("Simulation duration: {:.4f} seconds".format(time.time() - simStartTime))
    print("=========================================")

    if (SAVE_SIMULATION):
        solver = "rk4_"
        sdt = "dt0{}_".format(int(dt*to_actualRealType("10000")))

        fileNameStr = "{}/HH_BBT_{}{}{},{},vI{},t={}ms_{}_Iapp{}".format(outputDir, \
                solver, sdt, \
                nNeurons, nInhNeurons, \
                int(vInh), \
                p_maxTimeSimulation, \
                actualPrecisionType, iappOrder)
        writeToFile(fileNameStr + ".txt", ave)
        
        # Saving the Episodes start/end times
        writeToFile(fileNameStr + ",Epis.txt", networkOutputs);
            
        # Saving the applied currents values
        writeIappToFile(fileNameStr + ",Iapp.txt");
