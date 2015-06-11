%This function is to test the fixed points on domain wih several chains;


%Initialisation the framework parameters;

frameWorksParams = SimulationFrameworkParams('numSteps',Inf,'showSimulation',true,'dt',0.1,'numChains',3,'objectInteraction',false);

%definit a domain

diffConst = 0.01;
sphereForce = ForceManagerParams('diffusionForce',true,'diffusionConst',diffConst,...
                                 'lennardJonesForce',false,'LJPotentialWidth',0.01,'LJPotentialDepth',0.01,'dt',...
                                  frameWorksParams.simulator.dt);
dp          = DomainHandlerParams('domainShape','sphere','forceParams',sphereForce,'showDomain',true,'domainWidth',...
                                   15);
%definit a chain;

chainForce  = ForceManagerParams('springForce',true,'diffusionForce',true,'diffusionConst',diffConst,...
                                 'bendingElasticityForce',false,'bendingConst',0.1,'fixedParticleNum',3,...
                                 'springConst',0.05,...
                                 'minParticleEqDistance',0);
                             
cp(1)= ChainParams('numBeads',10,'b',1,'initializeInDomain',1,'springForce',true,...
                   'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,'fixedBeadNum',...
                   floor(1+(10-1)*rand),'fixedBeadsPosition',randn(1,3));

cp(2)= ChainParams('numBeads',15,'b',1,'initializeInDomain',1,'springForce',true,...
                   'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,'fixedBeadNum',...
                   floor(1+(15-1)*rand),'fixedBeadsPosition',cp(1).fixedBeadsPosition);

cp(3)= ChainParams('numBeads',20,'b',1,'initializeInDomain',1,'springForce',true,...
                   'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,'fixedBeadNum',...
                   floor(1+(20-1)*rand),'fixedBeadsPosition',cp(1).fixedBeadsPosition);
%registre parameters;
frameWorksParams.SetDomainParams(dp);
frameWorksParams.SetChainParams(cp);

% initialize simulator framework
r = RouseSimulatorFramework(frameWorksParams);
r.Run;
% numSteps = 1000;
% for sIdx = 1:numSteps
% 
% r.Step;
% end

