% This is a script to run the polymer simulation

% add Utils folder to the working path 
addpath(genpath(fullfile(pwd,'..','..','Utils')));

dimension           = 3;
numParticles        = 100;
dt                  = 0.001;
diffusionConst      = 1;
numSteps            = Inf;
numSimulations      = 1;
frictionCoefficient = 1;
connectedBeads      = [];
fixedBeads          = [];
<<<<<<< HEAD
metBeedNum          = [1 16 32];
=======
metBeedNum          = [2 14 27];
>>>>>>> 0749a1fe5c463e410ed2a63f45107acc9a3d5e8a
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
%                        numSteps,frictionCoefficient,connectedBeads,fixedBeads,b);
%                        
% 
%  rouseModelBis.Simulation;
%  
%  rouseModelBis.Plot