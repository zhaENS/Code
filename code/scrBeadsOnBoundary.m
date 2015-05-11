
addpath(genpath(fullfile(pwd,'..','..','PolymerChainDynamics')));
%addpath(genpath(fullfile(pwd,'..','..','utils')));
dp = DomainHandlerParams;
dp.domainShape    = 'sphere';
dp.domainWidth    = sqrt(14);
dp.domainCenter   = [0 0 0];
dp.showDomain     = true;
dp.dt             = 0.01;
dp.diffusionConst = 0.1;
initialPoint      = [2 1 3];
numSteps          = 20;


%points = DiffuseOnSphere(initialPoint,numSteps,radius, dt, diffusionConst);
domainClass = DomainHandler(dp);
beads       = [ 2 5 7 10];% beads on the boundary
%points      = DiffuseOnSphere(initialPoint,numSteps,dp.domainWidth, dp.dt, dp.diffusionConst)
%points      = BeadsOnBoundary(dp.domainWidth, dp.dt, dp.diffusionConst,beads);
numSteps    = 20;

path = BrownianBridgeSim(numSteps,domainClass,dp,beads);



 
 

