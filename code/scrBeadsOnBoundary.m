
addpath(genpath(fullfile(pwd,'..','..','PolymerChainDynamics')));
close all 
dp = DomainHandlerParams;
dp.domainShape    = 'sphere';
dp.domainWidth    = sqrt(14);
dp.domainCenter   = [0 0 0];
dp.showDomain     = true;
dp.dt             = 0.01;
dp.diffusionConst = 1;
%initialPoint      = [2 1 3];
numBeads          = 200;


%points = DiffuseOnSphere(initialPoint,numSteps,radius, dt, diffusionConst);
domainClass = DomainHandler(dp);
beads       = [1 3 60 90];% beads on the boundary
%points      = DiffuseOnSphere(initialPoint,numSteps,dp.domainWidth, dp.dt, dp.diffusionConst)
%points      = BeadsOnBoundary(dp.domainWidth, dp.dt, dp.diffusionConst,beads);

BrownianBridgeSim(domainClass,dp,beads,numBeads);



 
 

