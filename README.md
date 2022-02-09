# Mitigating computer´s limitations in replicating numerical simulations of a neural network model with Hodgkin-Huxley-type neurons
This project continues the study of Blanco 2020 (https://github.com/wblancof/neural-numerical-replicability) and brings a new analyses of the influence of numerical imprecision under different operating systems (*Windows and Ubuntu*) and across languages/softwares (*XPP/XPPAUT, C/C++, Python and Matlab*); in the context of neural computational replicability. The present study case is an 8 seconds biological simulation where a small neuronal network composed of 100 Hodgkin-Huxley-type neurons was modeled under Windows and Ubuntu OSs. The original code in C++ (Blanco 2020) was carefully altered and refactored into *XPP/XPPAUT, Python and Matlab*; aiming the replicability and mitigate non-replicability issues. The model parameters are described in the article methods section. The two main scenarios of Blanco 2020 were executed in all languages/softwares, i.e., simulations composed of (1) neurons 100% excitatory and (2) simulations composed of neurons 80% excitatory and 20% inhibitory. The Ordinary Differential Equations (ODEs) were solved by Runge–Kutta fourth-order method (RK4) and this method was implemented in C/C++, Python and Matlab languages, XPP/XPPAUT has its own RK4 and does not have support to outside implementations. To accurate results, was utilized the library Boost on C/C++ and Decimal on Python. XPP/XPPAUT has no support to high floating-point precision (HFPP). The HFPP in Matlab was impracticable due to the high time demand.

## Software requirements
- XPP/XPPAUT version 6.9.0.31
    - http://www.math.pitt.edu/~bard/xpp/xpp.html
- C/C++ version 9.3.0
    - https://www.cplusplus.com/
    - Boost C++ libraries version 1.72.0 (https://boostorg.jfrog.io/artifactory/main/release/1.72.0/source/)
- Python version 3.8.0
    - https://www.python.org/
- Matlab version 9.5.0 2017b
    - https://www.mathworks.com/products/matlab.html
- Make builder

## Project Structure
For better organization, the directories structure of this project is divided into languages/software with their respective codes. However, the code files performed equally for the Operating Systems utilized in this project (*Windows and Ubuntu*).
- XPP/XPPAUT
    * 1\. iappASC.txt
    * 2\. iappDESC.txt
    * 3\. iappRND.txt
    * 4\. xppExc.ode
    * 5\. xppInh.ode
- C/C++
    * 1\. HH_BBT2020_allP.cpp
    * 2\. Makefile
    * 3\. SimulationInitialization.h
    * 4\. SimulationParameters.h
    * 5\. iappInit.hpp
    * 6\. rk4.hpp
- Python
    * 1\. HH_BBT2020_allP.py
    * 2\. SimulationInitialization.py
    * 3\. SimulationParameters.py
    * 4\. iappInit.py
    * 5\. rk4.py
- Matlab
    * 1\. HH_BBT2020_allP.m
    * 2\. SimulationInitialization.m
    * 3\. SimulationParameters.m
    * 4\. iappInit.m
    * 5\. rk4.m

## How to execute simulations
- XPP/XPPAUT
    >>Note: XPP/XPPAUT does not have High floating-point precision.
    * 1\. Once the installation of the software was done, copy all files (*xppInh.ode, xppExc.ode, iappASC.txt, iappDESC.txt and iappRND.txt*) to the root directory XPP/XPPAUT. 
    * 2\. The neuron organization order (Iapp commutation) is set in the files *iappASC.txt, iappDESC.txt and iappRND.txt* and can be selected by comment/uncomment the lines 28-30 in executer file (*xppInh.ode and xppExc.ode*).
    * 3\. Open the XPP/XPPAUT software and browse the system to find *xppExc.ode* file and run simulation composed with neurons 100% Excitatory; or *xppInh.ode* to run simulation composed with neurons 80% Excitatory and 20% Inhibitory. It is also possible to drag and drop these files to XPP/XPPAUT to run the simulation. After the application initialization, click on `Initialconds` in the left vertical panel or the <kbd>I</kbd> key, then `Go` or <kbd>G</kbd> key.
    * 4\. After the execution, click in `Graphic stuff` <kbd>G</kbd> then `Export data` <kbd>O</kbd> to export simulation outcomes.
- C/C++
    * 1\. Edit lines 6 (LIBS) and 7 (CXXINCS) of make file *Makefile* to set the boost libs and boost library directory.
    * 2\. To execute simulations with Double or Boost precision, comment/uncomment the lines 6 (*#define DoublePrecision*) or 7 (*#define HighPrecision*) in file *SimulationParameters.h*.
    * 3\. Execute the make command and certify that HH_BBT2017_allP.exe file was created
    * 4\. The neuron organization order (Iapp commutation) is set in the file *SimulationParameters.h* and can be selected by changing line 12 (*const char iappOrder[] = "DES";*). The options are shown in line 11. DESC for Descendant organization, ASC for Ascendant and RAND for Random organization.
    * 5\. The number of Excitatory neurons can be set in line 23 (*const unsigned int nExcNeurons = 80;*). When this number is inferior to the total of neurons (total of 100, setting on line 22), the system define this difference as Inhibitory neurons.
    * 6\. Execute the following command in the terminal: ```./HH_BBT2020_allP.exe```
    * 7\. The simulation will create a folder named *results* on the same directory files and place the outcomes there
- Python
    * 1\. To execute simulations with Double or Decimal precision, comment/uncomment the lines 6 (*actualPrecisionType = "DoublePrecision"*) or 7 (*actualPrecisionType = "HighPrecision"*) in file *SimulationParameters.py*.
    * 2\. The neuron organization order (Iapp commutation) is set in the file *SimulationParameters.h* and can be selected by changing line 12 (*iappOrder = "DES"*). The options are shown in line 11. DESC for Descendant organization, ASC for Ascendant and RAND for Random organization.
    * 3\. The number of Excitatory neurons can be set in line 23 (*nExcNeurons=80*). When this number is inferior to the total of neurons (total of 100, setting on line 22), the system define this difference as Inhibitory neurons.
    * 4\. Execute the following command in the terminal: ```python HH_BBT2020_allP.py```
    * 5\. The simulation will create a folder named *results* on the same directory files and place the outcomes there

- Matlab
    >>Note: Matlab does not have High floating-point precision.
    * 1\. The neuron organization order (Iapp commutation) is set in the file *SimulationParameters.m* and can be selected by changing line 12 (*iappOrder = "DES"*). The options are shown in line 11. DESC for Descendant organization, ASC for Ascendant and RAND for Random organization.
    * 3\. The number of Excitatory neurons can be set in line 23 (*nExcNeurons=80*). When this number is inferior to the total of neurons (total of 100, setting on line 22), the system define this difference as Inhibitory neurons.
    * 4\. Open the files directory in Matlab program and execute the main file: ```HH_BBT2020_allP.m```
    * 5\. The simulation will create a folder named *results* on the same directory files and place the outcomes there
