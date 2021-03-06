% scr test connect/disconnect beads from boundary 


%This function is to test the movement of several chains with fixed beads
%on sphere and several beads moving on the boundary;


%Initialisation the framework parameters;

frameWorksParams = SimulationFrameworkParams('numSteps',Inf,'showSimulation',true,'dt',0.1,...
                    'numChains',10,'objectInteraction',false,'dimension',3,'encounterDist',0.2);
                    

%define a domain
diffConst = 0.1;
sphereForce = ForceManagerParams('diffusionForce',true,'diffusionConst',diffConst,...
                                 'lennardJonesForce',true,'LJPotentialWidth',0.01,'LJPotentialDepth',0.01,'dt',...
                                  frameWorksParams.simulator.dt);
dp          = DomainHandlerParams('domainShape','sphere','forceParams',sphereForce,'showDomain',true,'domainWidth',...
                                  25,'dt',frameWorksParams.simulator.dt,'diffusionConst',diffConst,'moveDomainType',...,
                                  'rotate','moveAngles',[linspace(0,2*pi,400);linspace(0,2*pi,400);zeros(1,400)]);
                              
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
                   'beadsOnBoundary',[1 64],'allowAttachToBoundary',false,'maxStepsOnBoundaryPerTime',10,...
                   'oscillationMag',0,...
                   'oscillationAng',linspace(0,2*pi,100));

%  cp(2)= ChainParams('numBeads',64,'initializeInDomain',1,'springForce',true,...
%                    'fixedBeadNum',32,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'beadsOnBoundary',[1 64],'allowAttachToBoundary',true); 
% 
%  cp(3)= ChainParams('numBeads',64,'initializeInDomain',1,'springForce',true,...
%                    'fixedBeadNum',32,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'beadsOnBoundary',[1 64],'allowAttachToBoundary',true);                
% 
%  cp(4)= ChainParams('numBeads',64,'initializeInDomain',1,'springForce',true,...
%                    'fixedBeadNum',32,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'beadsOnBoundary',[1 64],'allowAttachToBoundary',true); 
% 
%  cp(5)= ChainParams('numBeads',64,'initializeInDomain',1,'springForce',true,...
%                    'fixedBeadNum',32,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'beadsOnBoundary',[1 64],'allowAttachToBoundary',true); 
% 
%  cp(6)= ChainParams('numBeads',64,'initializeInDomain',1,'springForce',true,...
%                    'fixedBeadNum',32,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'beadsOnBoundary',[1 64],'allowAttachToBoundary',true);     
% 
%  cp(7)= ChainParams('numBeads',64,'initializeInDomain',1,'springForce',true,...
%                    'fixedBeadNum',32,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'beadsOnBoundary',[1 64],'allowAttachToBoundary',true); 
% 
%  cp(8)= ChainParams('numBeads',64,'initializeInDomain',1,'springForce',true,...
%                    'fixedBeadNum',32,'fixedBeadsPosition',initialPt,'diffusionConst',0.1,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'beadsOnBoundary',[1 64],'allowAttachToBoundary',true); 
%                
 
%registre parameters;
frameWorksParams.SetDomainParams(dp);
frameWorksParams.SetChainParams(cp);

% initialize simulator framework
r = RouseSimulatorFramework(frameWorksParams);
%run;
%profile on 
r.Run;
%profile viewer

