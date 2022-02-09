






SimulationParameters;
SimulationInitialization;

clc; clear;
global nNeurons;
% Create a table of applied currents Ii
% Table iapp % 100 0 99 i0+t*deli/99
global iapp; iapp = zeros(1,nNeurons);                   % for Excitatory neurons
global iappj; iappj = 0;

% Synaptic drive from excitatory cells
global atotExc; atotExc = 0.0;   global atotExcj; atotExcj = 0.0;
% Synaptic drive from inhibitory cells
global atotInh; atotInh = 0.0;   global atotInhj; atotInhj = 0.0;

% v[i] membrane potential of cell i
% n[i] activation of K+ conductance for cell i
% a[i] synaptic drive from cell i
% s[i] synaptic recovery for terminals from cell i

% Vector of neurons -> The Network
global network; network = zeros(nNeurons,4);

% SpikeTrain Class to manage/store/save the spike times for each neuron
% global spikes; spikes = struct('mSpikeTrain', []);

% Auxiliary variables for output
global ave; ave = zeros(800000,4);                             % Average of V N A and S
global networkOutputs; networkOutputs = zeros(800000,4);       % Average Ainh and Aexc, episodes time,
main;                                   % call the main function (for matlab only)
% ==========
% Functions
% ==========

% Rate constants and steady state activation function for Na+ and K+ currents
% am(v)=.1*(25-v)/(exp(.1*(25-v))-1)
function a = am(v)                
    a = 0.1*(25-v)/(exp(.1*(25-v))-1);
end

% bm(v)=4.0*exp(-v/18)
function b = bm(v)                
    b = 4.0*exp(-v/18);
end

% minf(v)=am(v)/(am(v)+bm(v))
function m = minf(v)              
    m = am(v)/(am(v)+bm(v));
end

% bn(v)=.125*exp(-v/80)
function b = bn(v)                
    b = .125*exp(-v/80);
end

% an(v)=.01*(10-v)/(exp(.1*(10-v))-1.0)
function a = an(v)                
    a = .01*(10-v)/(exp(.1*(10-v))-1.0);
end

% Synaptic activation
% fsyn(V) = 1./(1+exp((Vthresh-V)/kv))
function b = fsyn(V)              
    global Vthresh;
    b = 1./(1+exp((Vthresh-V))); 
end

% System ODEs per Neuron
% This function uses the following global variables
%  vna=115.0;
%  vk = -12.0;
%  vl=10.6;
%  gnabar=36.0;
%  gkbar=12.0;
%  gl=0.1;
%  h0=0.8;
%  iappj;
%  gsyn=3.6/nNeurons;
function dxdt = HH_NeuronModel(network)
    global vna; global vk; global vl; global gnabar; global gkbar; global gl; global h0; global gsyn; global atotExcj; 
    global atotInhj; global vExc; global vInh; global iappj; global tauf; global taus; global alphad; global betad;
    vj = network(1);
    nj = network(2);
    aj = network(3);
    sj = network(4);
   
    % Differential equations (XPP original version)
    % v[i] membrane potential of cell i
    % v[0..99]'= -gl*(v[j]-vl)-gnabar*minf(v[j])^3*(h0-n[j])*(v[j]-vna)-gkbar*n[j]^4*(v[j]-vk)
                    %-gsyn*(atot-a[j]*s[j]/100)*(v[j]-vsyn)
                    % + iapp([j])
    dxdt(1) = -gl*(vj-vl) ...
                  - gnabar*(minf(vj)^3)*(h0-nj)*(vj-vna) ...
                  - gkbar*(nj^4)*(vj-vk) ...
                  - gsyn*atotExcj*(vj-vExc) ...
                  - gsyn*atotInhj*(vj-vInh) ...
                  + iappj;
                                    
    % n[i] Activation of K+ conductance for cell i
    % n[0..99]'= an(v[j])-(an(v[j])+bn(v[j]))*n[j]
    dxdt(2) = an(vj)-(an(vj)+bn(vj))*nj;
    
    % a[i] Synaptic drive from cell i
    % a[0..99]'= fsyn(v[j])*(1-a[j])/tauf - a[j]/taus
    dxdt(3) = fsyn(vj)*(1-aj)/tauf - aj/taus;
    
    % s[i] Synaptic recovery for terminals from cell i
    % s[0..99]'=alphad*(1-s[j])-betad*fsyn(v[j])*s[j]
    dxdt(4) = alphad*(1-sj)-betad*fsyn(vj)*sj;
end

function [net depol] = initEachNeuron(net, depol) 
    
    global nNeurons;
    %Initial conditions
    % init v[0..99]=0 n[j]=0 a[j]=0.01 s[j]=.25
    for n = 1:nNeurons
        net(n,1) = 0.0;      % v
        net(n,2) = 0.0;      % n
        net(n,3) = 0.01;     % a
        net(n,4) = 0.25;     % s
        depol(n) = false;     % no spike found
    end
end

% Function to write a .txt file with the output parameters of the network state
function writeToFile(fileNameStr, v) 

   sprintf("Writing in file: " + fileNameStr + "\n")
   
   myfile = fopen(fileNameStr,'w');
   
   
   
   fprintf(myfile,'%2.14f\t%2.14f\t%2.14f\t%2.14f\n',v.');
   fclose(myfile);



end

% Function to write a .txt with the iapp values
function writeIappToFile(fileNameStr) 
   global iapp;
   sprintf("Writing in file: " + fileNameStr + "\n")

   
   myfile = fopen(fileNameStr, 'w');
        
   
   for i=1:nNeurons
        fprintf(myfile,'iapp[%d]= %2.14f\n', i, iapp(i));
   end
   fclose(myfile);



end

function par = showParameters( ) 
    global nNeurons; global nExcNeurons; global nInhNeurons; global maxTimeSimulation; 
    global maxNumBurst; global dt; global vInh; global n_precision; global iappOrder; 
    [v d] = version;  sprintf("Matlab version: " + v + "\n")
    sprintf("Running with the parameters:\n")
    sprintf("Total Neurons = %d \n", nNeurons)
    sprintf("Exc Neurons = %d \n", nExcNeurons)
    sprintf("Inh Neurons = %d \n", nInhNeurons)
    sprintf("maxTimeSimulation = %d s \n", maxTimeSimulation/1000)
    sprintf("nBurst = %d \n", maxNumBurst)
    sprintf("dt = %d ms \n", dt)
    sprintf("vInh = %d \n", vInh)
    sprintf("Set precision = %d \n", n_precision)
    sprintf("Organization = " + iappOrder + "\n")
end

% Check if the directory exist, case not, it will be created 
% =============================================================
function directory_existsCreate(s)
        
    [status, msg, msgID] = mkdir(s);
    
    if ~status
        sprintf("Error creating the " +  s + " directory! \n")
    end
    
    
end

function main
    t = 0;
    tspan = [0.0, 0.0];
    r = zeros(2,4);
    actTotalExc = 0;
    actTotalInh = 0;
    v_1 = 0.0;
    sTotalExc = 0;
    sTotalInh = 0;
    sN_1 = zeros(1, 4);
    sN = zeros(1, 4);
    sA = zeros(1, 4);
    %outputAct;
    simStartTime = 0;
    global outputDir iappOrder network depolarization SAVE_SIMULATION iapp ave networkOutputs;
    % Creating output logs directory
    % ============================
    directory_existsCreate(outputDir);

    % Initial conditions applied current for every neuron
    % The same distribution for all simulations
    % =========================================
    iapp = iappInit(iappOrder);

    % Simulation
    % ==================================
    sprintf("==> Simulation <== \n")
    showParameters()

    % Initial conditions for every neuron
    % ===================================
    [network, depolarization] = initEachNeuron(network, depolarization);
    
    simStartTime = tic;         % To get the simulation time 
    k = 1;                      % specific of Matlab
    global burstCount maxNumBurst maxTimeSimulation activePhase nInhNeurons nExcNeurons; 
    global nNeurons dt Vthresh thA thDA vInh atotExc atotInh atotExcj atotInhj iappj;
    while (burstCount < maxNumBurst) && (t <= maxTimeSimulation)
        % Total synaptic drive exc/inh for all cells
        atotExc = 0;
        atotInh = 0;

        % Synaptic drive from excitatory cells
        % atotexc=sum(0,100)of(shift(s0,i')*shift(a0,i'))/100
        % remembering the order: v, n, a, s -> 0, 1, 2, 3
        for n=1:nExcNeurons 
            atotExc = atotExc + network(n,4)*network(n,3);  % Order: v, n, a, s -> 1, 2, 3, 4
        end  
        % Synaptic drive from inhibitory cells
        for n=nExcNeurons+1:nNeurons  
            atotInh = atotInh + network(n,4)*network(n,3);  % Order: v, n, a, s -> 1, 2, 3, 4
        end

        % Saving the previous state of the network to detecting episodes
        sN_1(1) = 0.0; 
        sN_1(2) = 0.0;
        sN_1(3) = 0.0;
        sN_1(4) = 0.0;   
        for n=1:nNeurons 
            sN_1(1) = sN_1(1) + network(n,1);
            sN_1(2) = sN_1(2) + network(n,2);
            sN_1(3) = sN_1(3) + network(n,3);
            sN_1(4) = sN_1(4) + network(n,4);
        end

        sN_1(1) = sN_1(1)/nNeurons;
        sN_1(2) = sN_1(2)/nNeurons;
        sN_1(3) = sN_1(3)/nNeurons;
        sN_1(4) = sN_1(4)/nNeurons;
      
        % Initialize the Network state
        sN(1) = 0.0; 
        sN(2) = 0.0;
        sN(3) = 0.0; 
        sN(4) = 0.0;  

        % For each neuron
        for  n=1:nNeurons 
            % Synaptic drive Exc/Inh
            atotExcj = atotExc; 
            atotInhj = atotInh;    
            iappj = iapp(n);                           % Get the original distribution [-10 .. 5] iapp dist
            % Taking out the contribution of the own cell i
            if n <= nExcNeurons  % [Excitatory]
                atotExcj = atotExcj - (network(n,4)*network(n,3));    % remove its own synapse
            else                 % [Excitatory]
                atotInhj = atotInhj - (network(n,4)*network(n,3));    % remove its own synapse
            end
            
            v_1 = network(n,1);   % Save previous voltage, used to detect spikes
            tspan(1) = t;
            tspan(2) = t + dt; 
            % Integration, solving the ODE
            r = rk4 (@HH_NeuronModel, tspan, network(n,:), 1 );
            network(n,:) = r(:,2);

            
            
            
            % ==========================================================================================
            %        Detecting Depolarization
            %  It is not already on Depolarization period and Voltage >= Vth and Positive slope
            if ~depolarization(n) && network(n,1) >= Vthresh && (network(n,1)-v_1)/dt > 0 
                depolarization(n) = true;
            end
            % ==========================================================================================
            %        Detecting Repolarization
            %  It was on Depolarization period and Voltage <= Vth and Negative slope
            if depolarization(n) && network(n,1) <= Vthresh && (network(n,1)-v_1)/dt < 0  
                depolarization(n) = false;
                %addSpikeTimeToNeuron(n,t);
            end
 
            % Accumulating/Sum for each neuron state
            sN(1)= sN(1) + network(n,1);
            sN(2)= sN(2) + network(n,2);
            sN(3)= sN(3) + network(n,3);
            sN(4)= sN(4) + network(n,4);          
        end
        % Normalize the current state of the network for V A N and S
        sN(1) = sN(1)/nNeurons;
        sN(2) = sN(2)/nNeurons; 
        sN(3) = sN(3)/nNeurons; 
        sN(4) = sN(4)/nNeurons;
      
        % Population activity and synaptic recovery for excitatory neurons
        actTotalExc = 0; 
        sTotalExc = 0;
        for n=1:nNeurons 
            actTotalExc = actTotalExc + network(n,3);
            sTotalExc = sTotalExc + network(n,4);
        end
    
        % Population activity and synaptic recovery for Inhibitory neurons
        actTotalInh = 0; 
        sTotalInh = 0;
        for n = nExcNeurons+1:nNeurons 
            actTotalInh = actTotalInh + network(n,3);
            sTotalInh = sTotalInh + network(n,4);
        end
     
        % Detecting episode    
        % It is not already on active phase, activity > thA, activity derivative > thDA
        if ~activePhase && sN(3) >= thA && (sN(3)-sN_1(3))/dt > thDA
            activePhase = true;

            % Save values of Aexc, Ainh, time episode, 1
            sA(1) = actTotalExc/nExcNeurons;
            if (nInhNeurons > 0) sA(2) = actTotalInh/nInhNeurons; end
            sA(3) = t;
            sA(4) = 1;
            networkOutputs(k,:) = sA;          % Aexc, Ainh, time episode, 1
            % It was on active phase, activity < thA    
        elseif activePhase && sN(3) < thA              
            activePhase = false;
            burstCount = burstCount + 1;                  % Counting the episode

            
            %Save values of Aexc, Ainh, time episode, 1
            sA(1) = actTotalExc/nExcNeurons;
            if (nInhNeurons > 0) sA(2) = actTotalInh/nInhNeurons;   end
            sA(3) = -t;
            sA(4) = 1;
            networkOutputs(k,:) = sA;          % Aexc, Ainh, -time episode, 1
        end

        % Saving the current state of the network
        ave(k,:)= sN;

        t= t+dt;    % increment time
        k=k+1;
    end
    
    sprintf("Simulation duration: %10.10f seconds \n", toc(simStartTime))
    sprintf("========================================= \n")
    global p_maxTimeSimulation;
    if SAVE_SIMULATION
        solver = "rk4_";
        sdt = "dt0" + cast(dt*10000, 'int8') + "_";
        
        fileNameStr = outputDir + "/HH_BBT_" ...
                + solver + sdt ...
                + nNeurons + "," + nInhNeurons + ",vI" ...
                + cast(vInh, 'int8') + ",t=" ...
                + p_maxTimeSimulation + "ms_" ...
                + "double_Iapp" + iappOrder;
        writeToFile(fileNameStr + ".txt", ave);

        % Saving the Episodes start/end times
        writeToFile( fileNameStr + ",Epis.txt", networkOutputs);

        % Saving the applied currents values
        writeIappToFile(fileNameStr + ",Iapp.txt");
    end
    
end

