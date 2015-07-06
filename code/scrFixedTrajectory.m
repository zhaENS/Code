%This script is used to test the path of severals beads in the sphere when
%one bead is fixed on the sphere


%Initialisation the framework parameters;

frameWorksParams = SimulationFrameworkParams('numSteps',1000,'showSimulation',true,'dt',0.1,...
                    'numChains',1,'objectInteraction',false,'recipeFileName',...
                    'rcpTestParticleDetachmentFromDomain','recipesFolder','recipes',...
                     'encounterDist',0.6,'beadsPos',[]); 

%define a domain

diffConst = 0.1;
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

 
 cp(1)= ChainParams('numBeads',100,'b',1,'initializeInDomain',1,'springForce',true,...
                   'fixedBeadNum',[],'fixedBeadsPosition',[],...
                   'diffusionConst',diffConst,'forceParams',chainForce,...
                   'dt',frameWorksParams.simulator.dt,...
                   'beadsOnBoundary',[1],'stickyBeads',[]);

%registre parameters;
frameWorksParams.SetDomainParams(dp);
frameWorksParams.SetChainParams(cp);

% initialize simulator framework
r = RouseSimulatorFramework(frameWorksParams);
%run;
r.Run;




%the trajectory of the bead 10 and 50 during the simulation;
posBeads1 = r.params.simulator.beadsPos(1,:,:);
posBeads2 = r.params.simulator.beadsPos(2,:,:);

posBeads1 = reshape(posBeads1,3,r.params.simulator.numSteps-1);
posBeads1 = posBeads1';

posBeads2 = reshape(posBeads2,3,r.params.simulator.numSteps-1);
posBeads2 = posBeads2';



figure,
cameratoolbar
c = [rand rand rand];
plot3(posBeads1(:,1), posBeads1(:,2),posBeads1(:,3),'LineWidth',2);
legend('the trajectory of bead 10');
figure,
cameratoolbar
plot3(posBeads2(:,1), posBeads2(:,2),posBeads2(:,3),'Color',c,'LineWidth',2);
legend('the trajectory of bead 50');

daspect([1 1 1])