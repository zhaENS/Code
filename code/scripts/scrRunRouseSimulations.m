% This is a script to run the polymer simulation

% add Utils folder to the working path 
addpath(genpath(fullfile(pwd,'..','..','Utils')));

dimension           = 3;
numParticles        = 100;
dt                  = 0.01;
diffusionConst      = 10;
numSteps            = 100;
numSimulations      = 1;
frictionCoefficient = 1;
connectedBeads      = [];
fixedBeads          = [];
%metBeadNum1          = [1 16 32];
metBeadNum1          = [ones(1,30);[2:31];32*ones(1,30)]';
b                     = 1;
encounterDistance     = 3*b./5;
%initialize the class

rouseModel = RouseModelMetTime(dimension,numParticles,dt,diffusionConst,...
                numSteps,frictionCoefficient,connectedBeads, fixedBeads ,...
                numSimulations,metBeadNum1,b,encounterDistance);
            
        
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