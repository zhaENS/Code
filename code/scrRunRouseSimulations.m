% This is a script to run the polymer simulation

% add Utils folder to the working path 
addpath(genpath(fullfile(pwd,'..','..','Utils')));

dimension           = 3;
numParticles        = 32;
dt                  = 0.01;
diffusionConst      = 1;
numSteps            = 100;
numSimulations      = 1000;
frictionCoefficient = 1;
connectedBeads      = [];
fixedBeads          = [];
%metBead          = [1 16 32];
metBeadNum          = [ones(1,14);[2:15];32*ones(1,14)]';
b                   = 1;
encounterDistance   = b./5;
%initialize the class

% rouseModel = RouseModelMetTime(dimension,numParticles,dt,diffusionConst,...
%                 numSteps,frictionCoefficient,connectedBeads, fixedBeads ,...
%                 numSimulations,metBeadNum,b,encounterDistance);
%             
%         
% % run the simulation                    
% rouseModel.Simulation;
% 
% % plot 
%  rouseModel.Plot

 
 
 rouseModelBis = RouseModel(dimension,numParticles,dt,diffusionConst,...
                       numSteps,frictionCoefficient,connectedBeads,fixedBeads,b);
                       

 rouseModelBis.Simulation;
 
 rouseModelBis.Plot