
addpath(genpath(fullfile(pwd,'..','..','PolymerChainDynamics')));
addpath(genpath(fullfile(pwd,'..','..','Utils')));% import Utils 

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


%points = DiffuseOnSphere(initialPoint,numSteps,radius, dt, diffusionConst);
domainClass = DomainHandler(dp);
beads       = [1 50 90 100];% beads on the boundary
%points      = DiffuseOnSphere(initialPoint,numSteps,dp.domainWidth, dp.dt, dp.diffusionConst)
%points      = BeadsOnBoundary(dp.domainWidth, dp.dt, dp.diffusionConst,beads);

chainPath= BrownianBridgeSim(domainClass,dp,beads,numBeads);



 
 

