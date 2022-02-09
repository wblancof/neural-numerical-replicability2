

import math
from SimulationParameters import *


if(actualPrecisionType == "HighPrecision"):
    from decimal import *
    actualRealType = Decimal
    getcontext().prec = 100
    def to_actualRealType(str):
        return actualRealType(str)


    def expFun(x):
        return(x).exp()
else:
    actualRealType = float
    from math import exp as expFun
    def to_actualRealType(str):
        return actualRealType(str)



def division(x, y):
    if(y == 0):
        return math.copysign(actualRealType('inf'), x) if x != 0 else actualRealType('-nan')
    return x / y

# Simulation constants
dt = to_actualRealType(p_dt)
maxTimeSimulation = to_actualRealType(p_maxTimeSimulation)

# Simulation parameters
nInhNeurons = nNeurons - nExcNeurons                           # Amount of inhibitory neurons

# Parameters needed to detect episodes
thA = to_actualRealType(p_thA)
thDA = to_actualRealType(p_thDA)
activePhase = False
burstCount = 0

# Parameters needed to detect spikes
depolarization = [None] * nNeurons

# Cellular parameters
deli = to_actualRealType(p_deli)
vna = to_actualRealType(p_vna)
vk = to_actualRealType(p_vk)
vl = to_actualRealType(p_vl)
gnabar = to_actualRealType(p_gnabar)
gkbar = to_actualRealType(p_gkbar)
gl = to_actualRealType(p_gl)
h0 = to_actualRealType(p_h0)
vExc = to_actualRealType(p_vExc)
vInh = to_actualRealType(p_vInh)
vsyn = to_actualRealType(p_vsyn)

# Synaptic parameters
taus = to_actualRealType(p_taus)
tauf = to_actualRealType(p_tauf)
gsyn = to_actualRealType(p_gsyn) / nNeurons
alphad = to_actualRealType(p_alphad)
betad = to_actualRealType(p_betad)
Vthresh = to_actualRealType(p_Vthresh)
kv = to_actualRealType(p_kv)

# Constants inside functions
dot1 = to_actualRealType("0.1")
dot01 = to_actualRealType("0.01")
dot25 = to_actualRealType("0.25")
dot125 = to_actualRealType("0.125")
zero = to_actualRealType("0")
int1 = to_actualRealType("1")
int2 = to_actualRealType("2")
int3 = to_actualRealType("3")
int4 = to_actualRealType("4")
int6 = to_actualRealType("6")
int10 = to_actualRealType("10")
int18 = to_actualRealType("18")
int25 = to_actualRealType("25")
int80 = to_actualRealType("80")
