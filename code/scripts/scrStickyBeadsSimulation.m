%this test is used to test sticky beads connected and disconncted on sphere

frameWorksParams = SimulationFrameworkParams('numSteps',Inf,'showSimulation',true,'dt',0.1,...
                    'numChains',10,'objectInteraction',false,'dimension',3,...
                     'encounterDist',0.6); 

%define a domain
diffConst = 1;
sphereForce = ForceManagerParams('diffusionForce',true,'diffusionConst',diffConst,...
                                 'lennardJonesForce',true,'LJPotentialWidth',0.01,'LJPotentialDepth',0.01,'dt',...
                                  frameWorksParams.simulator.dt);
dp          = DomainHandlerParams('domainShape','sphere','forceParams',sphereForce,'showDomain',true,'domainWidth',...
                                  25,'dt',frameWorksParams.simulator.dt,'diffusionConst',diffConst);
                              
dh          =DomainHandler(dp);                             
%define a chain;

chainForce  = ForceManagerParams('springForce',true,'diffusionForce',true,'diffusionConst',diffConst,...
                                 'bendingElasticityForce',false,'bendingConst',0.1,...
                                 'springConst',1,...
                                 'minParticleEqDistance',0);
 initialPt = dh.GetRandomBoundarySample(1);
 

 
cp(1)= ChainParams('numBeads',64,'initializeInDomain',1,'springForce',true,...
                   'fixedBeadNum',32,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
                   'forceParams',chainForce,...
                   'b', sqrt(frameWorksParams.simulator.dimension),...
                   'dt',frameWorksParams.simulator.dt,...
                   'beadsOnBoundary',[1 64],'stickyBeads',[1,64]); 

cp(2)= ChainParams('numBeads',64,'b',sqrt(frameWorksParams.simulator.dimension),'initializeInDomain',1,'springForce',true,...
                   'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,...
                    'fixedBeadNum',32,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
                    'beadsOnBoundary',[1 64],'stickyBeads',[1,64]);

cp(3)= ChainParams('numBeads',64,'b',sqrt(frameWorksParams.simulator.dimension),'initializeInDomain',1,'springForce',true,...
                   'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,'fixedBeadsPosition',initialPt,...
                   'diffusionConst',diffConst,'fixedBeadNum',32,'beadsOnBoundary',[1 64],'stickyBeads',[1 64]);

     
               
cp(4)= ChainParams('numBeads',64,'b',sqrt(frameWorksParams.simulator.dimension),'initializeInDomain',1,'springForce',true,...
                   'diffusionConst',0.1,'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,'fixedBeadsPosition',initialPt,...
                   'fixedBeadNum',32,'beadsOnBoundary',[1 64],'stickyBeads',[1,64]);

 cp(5)= ChainParams('numBeads',64,'b',sqrt(frameWorksParams.simulator.dimension),'initializeInDomain',1,'springForce',true,...
                   'diffusionConst',diffConst,'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,'fixedBeadsPosition',initialPt,...
                   'fixedBeadNum',32,'beadsOnBoundary',[1 64]);

 
 cp(6)= ChainParams('numBeads',64,'b',sqrt(frameWorksParams.simulator.dimension),'initializeInDomain',1,'springForce',true,...
                   'diffusionConst',diffConst,'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,'fixedBeadsPosition',initialPt,...
                   'fixedBeadNum',32,'beadsOnBoundary',[1 64]);

 cp(7)= ChainParams('numBeads',64,'b',sqrt(frameWorksParams.simulator.dimension),'initializeInDomain',1,'springForce',true,...
                   'diffusionConst',diffConst,'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,'fixedBeadsPosition',initialPt,...
                   'fixedBeadNum',32,'beadsOnBoundary',[1 64]);

 cp(8)= ChainParams('numBeads',64,'b',sqrt(frameWorksParams.simulator.dimension),'initializeInDomain',1,'springForce',true,...
                   'diffusionConst',diffConst,'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,'fixedBeadsPosition',initialPt,...
                   'fixedBeadNum',32,'beadsOnBoundary',[1 64]);
%registre parameters;
frameWorksParams.SetDomainParams(dp);
frameWorksParams.SetChainParams(cp);

% initialize simulator framework
r = RouseSimulatorFramework(frameWorksParams);
%run;
r.Run;


