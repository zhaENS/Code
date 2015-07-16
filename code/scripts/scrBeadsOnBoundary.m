
%addpath(genpath(fullfile(pwd,'..','..','PolymerChainDynamics')));
%addpath(genpath(fullfile(pwd,'..','..','Utils')));% import Utils 

close all 
dp                = DomainHandlerParams;
dp.dimension      = 3;
dp.domainShape    = 'sphere';
dp.domainWidth    = 2;
dp.domainCenter   = [0 0 0];
dp.showDomain     = true;
dp.dt             = 0.01;
dp.diffusionConst = 1;
numBeads          = 100;

domainClass           = DomainHandler(dp);
beadsOnBoundary       = [5 50 100];% beads on the boundary

initialPoint = domainClass.GetRandomBoundarySample(1);
%chainPath= BrownianBridgeSim(domainClass,dp,beads,numBeads);
chainPath = BrownianBridgeSim(initialPoint,domainClass,dp,beadsOnBoundary,numBeads)


 
 

