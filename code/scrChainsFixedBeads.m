%This function is to test the movement of several chains with fixed beads
%on sphere and serveral beads moving on the boundary;


%Initialisation the framework parameters;

frameWorksParams = SimulationFrameworkParams('numSteps',Inf,'showSimulation',true,'dt',0.1,...
                    'numChains',10,'objectInteraction',false,...
                     'encounterDist',0.6,'stickyTime',zeros(1,40)); 

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
                                 'springConst',sqrt(frameWorksParams.domain.dimension),...
                                 'minParticleEqDistance',0);
 initialPt = dh.GetRandomBoundarySample(1);
 

 
 cp(1)= ChainParams('numBeads',64,'b',sqrt(frameWorksParams.domain.dimension),'initializeInDomain',1,'springForce',true,...
                   'fixedBeadNum',32,'fixedBeadsPosition',initialPt,'springConst',0.1,...
                   'forceParams',chainForce,...
                   'dt',frameWorksParams.simulator.dt,...
                   'beadsOnBoundary',[1 64],'stickyBeads',[1:10:64]);
               

     

cp(2)= ChainParams('numBeads',64,'b',sqrt(frameWorksParams.domain.dimension),'initializeInDomain',1,'springForce',true,...
                   'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,...
                    'fixedBeadNum',32,'fixedBeadsPosition',initialPt,'springConst',0.1,...
                    'beadsOnBoundary',[1 64],'stickyBeads',[1:10:64]);

cp(3)= ChainParams('numBeads',64,'b',sqrt(frameWorksParams.domain.dimension),'initializeInDomain',1,'springForce',true,...
                   'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,'fixedBeadsPosition',initialPt,...
                   'springConst',0.1,'fixedBeadNum',32,'beadsOnBoundary',[1 64],'stickyBeads',[1:10:64]);

     
               
cp(4)= ChainParams('numBeads',64,'b',sqrt(frameWorksParams.domain.dimension),'initializeInDomain',1,'springForce',true,...
                   'springConst',0.1,'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,'fixedBeadsPosition',initialPt,...
                   'fixedBeadNum',32,'beadsOnBoundary',[1 64],'stickyBeads',[1:10:64]);

 
cp(5)= ChainParams('numBeads',64,'b',sqrt(frameWorksParams.domain.dimension),'initializeInDomain',1,'springForce',true,...
                    'springConst',0.1,'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,'fixedBeadsPosition',initialPt,...
                   'fixedBeadNum',32,'beadsOnBoundary',[1 64],'stickyBeads',[1:10:64]);

 
cp(6)= ChainParams('numBeads',64,'b',sqrt(frameWorksParams.domain.dimension),'initializeInDomain',1,'springForce',true,...
                 'springConst',0.1,'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,'fixedBeadsPosition',initialPt,...
                   'fixedBeadNum',32,'beadsOnBoundary',[1 64],'stickyBeads',[1:10:64]);
% 
%  cp(7)= ChainParams('numBeads',64,'b',1,'initializeInDomain',1,'springForce',true,...
%                    'diffusionConst',diffConst,'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,'fixedBeadsPosition',initialPt,...
%                    'fixedBeadNum',32,'beadsOnBoundary',[1 64],'stickyBeads',[64]);
% 
%  cp(8)= ChainParams('numBeads',64,'b',1,'initializeInDomain',1,'springForce',true,...
%                    'diffusionConst',diffConst,'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,'fixedBeadsPosition',initialPt,...
%                    'fixedBeadNum',32,'beadsOnBoundary',[1 64],'stickyBeads',[64]);
% 
%  cp(9)= ChainParams('numBeads',64,'b',1,'initializeInDomain',1,'springForce',true,...
%                    'diffusionConst',diffConst,'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,'fixedBeadsPosition',initialPt,...
%                    'fixedBeadNum',32,'beadsOnBoundary',[1 64],'stickyBeads',[64]);
% 
%  cp(10)= ChainParams('numBeads',64,'b',1,'initializeInDomain',1,'springForce',true,...
%                    'diffusionConst',diffConst,'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,'fixedBeadsPosition',initialPt,...
%                    'fixedBeadNum',32,'beadsOnBoundary',[1 64],'stickyBeads',[64]);

 
%registre parameters;
frameWorksParams.SetDomainParams(dp);
frameWorksParams.SetChainParams(cp);

% initialize simulator framework
r = RouseSimulatorFramework(frameWorksParams);
%run;
r.Run;


