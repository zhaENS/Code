%this test is used to test the number of clusters of sticky beads and the
%time when they form;

%Initialisation the framework parameters;

frameWorksParams = SimulationFrameworkParams('numSimulations',100,'numSteps',Inf,'showSimulation',true,'dt',0.1,...
                    'numChains',8,'objectInteraction',false,'dimension',3,'encounterDist',sqrt(3)/5,...,
                    'recipeFileName','rcpRelaxationTimeStickyBeads');
                    

%define a domain
diffConst = 0.01;
sphereForce = ForceManagerParams('diffusionForce',true,'diffusionConst',diffConst,...
                                 'lennardJonesForce',false,'LJPotentialWidth',0.01,'LJPotentialDepth',0.01,'dt',...
                                  frameWorksParams.simulator.dt);
dp          = DomainHandlerParams('domainShape','sphere','forceParams',sphereForce,'showDomain',true,'domainWidth',...
                                  5,'dt',frameWorksParams.simulator.dt,'diffusionConst',diffConst);
                              
dh          = DomainHandler(dp);                             
%define a chain;

chainForce  = ForceManagerParams('springForce',true,'diffusionForce',false,'diffusionConst',diffConst,...
                                 'bendingElasticityForce',false,'bendingConst',0.1,...
                                 'bendingOpeningAngle',pi,...
                                 'springConst',1,...
                                 'minParticleEqDistance',0,...
                                 'stickyParticlesSpringConst',1);
 initialPt = dh.GetRandomBoundarySample(1);
 

%  for cIdx=1:frameWorksParams.simulator.numChains
%  cp(cIdx)= ChainParams('numBeads',65,'initializeInDomain',1,...
%                    'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',diffConst,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'stickyBeads',[],'beadsOnBoundary',[1 65],'maxStepsOnBoundaryPerTime',10);
%   
%  end
for cIdx = 1:frameWorksParams.simulator.numChains
 cp(cIdx)= ChainParams('numBeads',65,'initializeInDomain',1,...
                   'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',diffConst,...
                   'forceParams',chainForce,...
                   'b', sqrt(frameWorksParams.simulator.dimension),...
                   'dt',frameWorksParams.simulator.dt,...
                   'stickyBeads',[]);
end
% 
% 
%  cp(2)= ChainParams('numBeads',65,'initializeInDomain',1,...
%                    'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',diffConst,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'stickyBeads',[]); 
% 
%  cp(3)= ChainParams('numBeads',65,'initializeInDomain',1,...
%                    'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',diffConst,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'stickyBeads',[]);
% 
%  cp(4)= ChainParams('numBeads',65,'initializeInDomain',1,...
%                    'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',diffConst,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'stickyBeads',[]); 
% 
%  cp(5)= ChainParams('numBeads',65,'initializeInDomain',1,...
%                    'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',diffConst,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'stickyBeads',[]); 
% 
%  cp(6)= ChainParams('numBeads',65,'initializeInDomain',1,...
%                    'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',diffConst,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'stickyBeads',[]);     
% 
%  cp(7)= ChainParams('numBeads',65,'initializeInDomain',1,...
%                    'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',diffConst,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'stickyBeads',[]); 
% 
%  cp(8)= ChainParams('numBeads',65,'initializeInDomain',1,...
%                    'fixedBeadNum',33,'fixedBeadsPosition',initialPt,'diffusionConst',diffConst,...
%                    'forceParams',chainForce,...
%                    'b', sqrt(frameWorksParams.simulator.dimension),...
%                    'dt',frameWorksParams.simulator.dt,...
%                    'stickyBeads',[]); 
               
               
%registre parameters;
frameWorksParams.SetDomainParams(dp);
frameWorksParams.SetChainParams(cp);

% initialize simulator framework
r = RouseSimulatorFramework(frameWorksParams);
%run;
% profile on
r.Run;
% profile viewer

