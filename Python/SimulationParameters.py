


# Real number type representation
# Leave ONLY ONE of the lines bellow uncommented
# actualPrecisionType = "DoublePrecision"
actualPrecisionType = "HighPrecision"


# Input/Output parameters
# RAND, ASC, DES
iappOrder = "DES"
SAVE_SIMULATION = True
n_precision = 14
outputDir = "./results"

# Simulation constants
p_dt = "0.01"                         # Time step: ms
p_maxTimeSimulation = "8000"          # Maximum simulation time: ms

# Simulation parameters
nNeurons = 100
nExcNeurons = 80                      # Amount of excitatory neurons
maxNumBurst = 200                     # Amount of burst to be reached to stop the simulation

# Parameters needed to detect episodes
p_thA = "0.1730"                      # thA = 0.25*(maxAt-minAt) based on Patrick paper. Previous calculated in Matlab
p_thDA = "0.1490"

# Cellular parameters
# p deli=15
# p vna=115  vk=-12  vl=10.6  gnabar=36  gkbar=12  gl=0.1
# p h0=0.8
p_deli = "15.0"
p_vna = "115.0"
p_vk  = "-12.0"
p_vl = "10.6"
p_gnabar = "36.0"
p_gkbar = "12.0"
p_gl = "0.1"
p_h0 = "0.8"

# p vexc=70
# p vinh = -12
p_vExc = "70.0"
p_vInh = "70.0"
p_vsyn = "70.0"

# Synaptic parameters
# p taus=10 tauf=1
# p gsyn=3.6 vsyn=70
# p alphad=0.0015 betad=0.12
# p Vthresh=40 kv=1
p_taus = "10.0"
p_tauf = "1.0"
p_gsyn = "3.6"
p_alphad = "0.0015"
p_betad = "0.12"
p_Vthresh = "40.0"
p_kv = "1.0"
