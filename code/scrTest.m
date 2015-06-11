% This is a script to run the polymer simulation

% add Utils folder to the working path 
addpath(genpath(fullfile(pwd,'..','..','Utils')));

dimension           = 3;
numParticles        = 32;
dt                  = 0.01;
diffusionConst      = 10;
numSteps            = Inf;
numSimulations      = 1000;
frictionCoefficient = 1;
connectedBeads      = [];
fixedBeads          = [];
metBeadNum          =[1 32];
metBeadNum1         =[1 16 32];
%metBeadNum1          = [ones(1,15);2:16;32*ones(1,15)]';
b                   = 1;
encounterDistance   = 3*b./5;
%initialize the class

rouseModelBis = RouseModelMetTime(dimension,numParticles,dt,diffusionConst,...
                numSteps,frictionCoefficient,connectedBeads, fixedBeads ,...
                numSimulations,metBeadNum,metBeadNum1,b,encounterDistance);
            
        
% run the simulation                    
rouseModelBis.Simulation;

% plot 
 rouseModelBis.Plot