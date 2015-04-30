% This is a script to run the polymer simulation

% add Utils folder to the working path 
addpath(genpath(fullfile(pwd,'..','..','Utils')));

dimension           = 3;
numParticles        = 32;
dt                  = 0.01;
diffusionConst      = 1;
numSteps            = Inf;
numSimulations      = 1;
frictionCoefficient = 1;
connectedBeads      = [];
fixedBeads          = [];
metBeedNum          = [1 16 32];
b                   = 1;
encounterDistance   = b./2;
%initialize the class

rouseModel = RouseModelMetTime(dimension,numParticles,dt,diffusionConst,...
                numSteps,frictionCoefficient,connectedBeads, fixedBeads ,...
                numSimulations,metBeedNum,b,encounterDistance);
            
        
% run the simulation                    
rouseModel.Simulation;

% plot 
 rouseModel.Plot

 
 
%  rouseModelBis = RouseModel(dimension,numParticles,dt,diffusionConst,...
%                        numSteps,frictionCoefficient,connectedBeads,fixedBeads);
%                        
% 
%  rouseModelBis.Simulation;
%  
%  rouseModelBis.Plot