%this test is used to test the number of clusters of sticky beads and the
%time when they form;

%Initialisation the framework parameters;

frameWorksParams = SimulationFrameworkParams('numSimulations',1,'numSteps',6000,'showSimulation',true,'dt',0.1,...
                    'numChains',10,'objectInteraction',false,'dimension',3,'encounterDist',sqrt(3)/5,...
                    'recipeFileName','rcpTestNumClustersStickyBeads');
                    

%define a domain
diffConst = 0.1;
sphereForce = ForceManagerParams('diffusionForce',true,'diffusionConst',diffConst,...
                                 'lennardJonesForce',false,'LJPotentialWidth',0.01,'LJPotentialDepth',0.01,'dt',...
                                  frameWorksParams.simulator.dt);
dp          = DomainHandlerParams('domainShape','sphere','forceParams',sphereForce,'showDomain',true,'domainWidth',...
                                  25,'dt',frameWorksParams.simulator.dt,'diffusionConst',diffConst);
                              
dh          = DomainHandler(dp);                             
%define a chain;

chainForce  = ForceManagerParams('springForce',true,'diffusionForce',true,'diffusionConst',diffConst,...
                                 'bendingElasticityForce',false,'bendingConst',0.1,...
                                 'springConst',1,...
                                 'minParticleEqDistance',0);
 initialPt = dh.GetRandomBoundarySample(1);
 

 
 cp(1)= ChainParams('numBeads',65,'initializeInDomain',1,'springForce',true,...
                   'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
                   'forceParams',chainForce,...
                   'b', sqrt(frameWorksParams.simulator.dimension),...
                   'dt',frameWorksParams.simulator.dt,...
                   'beadsOnBoundary',[1 65],'stickyBeads',[1,65],'allowAttachToBoundary',true);

 cp(2)= ChainParams('numBeads',65,'initializeInDomain',1,'springForce',true,...
                   'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
                   'forceParams',chainForce,...
                   'b', sqrt(frameWorksParams.simulator.dimension),...
                   'dt',frameWorksParams.simulator.dt,...
                   'beadsOnBoundary',[1 65],'stickyBeads',[1,65],'allowAttachToBoundary',true); 

 cp(3)= ChainParams('numBeads',65,'initializeInDomain',1,'springForce',true,...
                   'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
                   'forceParams',chainForce,...
                   'b', sqrt(frameWorksParams.simulator.dimension),...
                   'dt',frameWorksParams.simulator.dt,...
                   'beadsOnBoundary',[1 65],'stickyBeads',[1,65],'allowAttachToBoundary',true);

 cp(4)= ChainParams('numBeads',65,'initializeInDomain',1,'springForce',true,...
                   'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
                   'forceParams',chainForce,...
                   'b', sqrt(frameWorksParams.simulator.dimension),...
                   'dt',frameWorksParams.simulator.dt,...
                   'beadsOnBoundary',[1 65],'stickyBeads',[1,65],'allowAttachToBoundary',true); 

 cp(5)= ChainParams('numBeads',65,'initializeInDomain',1,'springForce',true,...
                   'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
                   'forceParams',chainForce,...
                   'b', sqrt(frameWorksParams.simulator.dimension),...
                   'dt',frameWorksParams.simulator.dt,...
                   'beadsOnBoundary',[1 65],'stickyBeads',[1,65],'allowAttachToBoundary',true); 

 cp(6)= ChainParams('numBeads',65,'initializeInDomain',1,'springForce',true,...
                   'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
                   'forceParams',chainForce,...
                   'b', sqrt(frameWorksParams.simulator.dimension),...
                   'dt',frameWorksParams.simulator.dt,...
                   'beadsOnBoundary',[1 65],'stickyBeads',[1,65],'allowAttachToBoundary',true);     

 cp(7)= ChainParams('numBeads',65,'initializeInDomain',1,'springForce',true,...
                   'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
                   'forceParams',chainForce,...
                   'b', sqrt(frameWorksParams.simulator.dimension),...
                   'dt',frameWorksParams.simulator.dt,...
                   'beadsOnBoundary',[1 65],'stickyBeads',[1,65],'allowAttachToBoundary',true); 

 cp(8)= ChainParams('numBeads',65,'initializeInDomain',1,'springForce',true,...
                   'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
                   'forceParams',chainForce,...
                   'b', sqrt(frameWorksParams.simulator.dimension),...
                   'dt',frameWorksParams.simulator.dt,...
                   'beadsOnBoundary',[1 65],'stickyBeads',[1,65],'allowAttachToBoundary',true); 
               
%   cp(9)= ChainParams('numBeads',65,'initializeInDomain',1,'springForce',true,...
%                    'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'beadsOnBoundary',[1 65],'stickyBeads',[1,65]);
% 
%  cp(10)= ChainParams('numBeads',65,'initializeInDomain',1,'springForce',true,...
%                    'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'beadsOnBoundary',[1 65],'stickyBeads',[1,65]); 
% 
%  cp(11)= ChainParams('numBeads',65,'initializeInDomain',1,'springForce',true,...
%                    'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'beadsOnBoundary',[1 65],'stickyBeads',[1,65]);                
% 
%  cp(12)= ChainParams('numBeads',65,'initializeInDomain',1,'springForce',true,...
%                    'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'beadsOnBoundary',[1 65],'stickyBeads',[1,65]); 
% 
%  cp(13)= ChainParams('numBeads',65,'initializeInDomain',1,'springForce',true,...
%                    'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'beadsOnBoundary',[1 65],'stickyBeads',[1,65]); 
% 
%  cp(14)= ChainParams('numBeads',65,'initializeInDomain',1,'springForce',true,...
%                    'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'beadsOnBoundary',[1 65],'stickyBeads',[1,65]);     
% 
%  cp(15)= ChainParams('numBeads',65,'initializeInDomain',1,'springForce',true,...
%                    'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'beadsOnBoundary',[1 65],'stickyBeads',[1,65]); 
% 
%  cp(16)= ChainParams('numBeads',65,'initializeInDomain',1,'springForce',true,...
%                    'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'beadsOnBoundary',[1 65],'stickyBeads',[1,65]); 
               
%registre parameters;
frameWorksParams.SetDomainParams(dp);
frameWorksParams.SetChainParams(cp);

% initialize simulator framework
r = RouseSimulatorFramework(frameWorksParams);
%run;
% profile on
r.Run;
% profile viewer

