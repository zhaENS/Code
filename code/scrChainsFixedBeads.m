%This function is to test the fixed points on domain wih several chains;


%Initialisation the framework parameters;

frameWorksParams = SimulationFrameworkParams('numSteps',Inf,'showSimulation',true,'dt',0.1,'numChains',3,'objectInteraction',false);

%define a domain

diffConst = 0.01;
sphereForce = ForceManagerParams('diffusionForce',true,'diffusionConst',diffConst,...
                                 'lennardJonesForce',true,'LJPotentialWidth',0.01,'LJPotentialDepth',0.01,'dt',...
                                  frameWorksParams.simulator.dt);
dp          = DomainHandlerParams('domainShape','sphere','forceParams',sphereForce,'showDomain',true,'domainWidth',...
                                  25,'dt',frameWorksParams.simulator.dt,'diffusionConst',diffConst);
                              
dh          =DomainHandler(dp);                             
%define a chain;

chainForce  = ForceManagerParams('springForce',true,'diffusionForce',true,'diffusionConst',diffConst,...
                                 'bendingElasticityForce',false,'bendingConst',0.1,...
                                 'springConst',0.05,...
                                 'minParticleEqDistance',0);
 initialPt = dh.GetRandomBoundarySample(1);
cp(1)= ChainParams('numBeads',30,'b',1,'initializeInDomain',1,'springForce',true,...
                   'fixedBeadNum',[15],'fixedBeadsPosition',initialPt,...
                   'forceParams',chainForce,...
                   'dt',frameWorksParams.simulator.dt,...
                   'beadsOnBoundary',[1 15 30]);
               

     

cp(2)= ChainParams('numBeads',20,'b',1,'initializeInDomain',1,'springForce',true,...
                   'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,...
                    'fixedBeadNum',[10],'fixedBeadsPosition',initialPt,...
                    'beadsOnBoundary',[1 10 20]);

cp(3)= ChainParams('numBeads',15,'b',1,'initializeInDomain',1,'springForce',true,...
                   'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,'fixedBeadsPosition',initialPt,...
                   'fixedBeadNum',[8],'beadsOnBoundary',[1 8 15]);

                   


%registre parameters;
frameWorksParams.SetDomainParams(dp);
frameWorksParams.SetChainParams(cp);

% initialize simulator framework
r = RouseSimulatorFramework(frameWorksParams);
%run;
r.Run;


