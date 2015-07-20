% test simulation of piolymers in two different domains 

frameWorksParams = SimulationFrameworkParams('numSteps',Inf,'showSimulation',true,'dt',0.1,...
                    'numChains',10,'objectInteraction',false,'dimension',3,...
                     'encounterDist',0.6); 

%define a domain
            
diffConst = 1.0;
sphereForce = ForceManagerParams('diffusionForce',true,'diffusionConst',diffConst,...
                                 'lennardJonesForce',true,'LJPotentialWidth',0.01,'LJPotentialDepth',0.01,'dt',...
                                  frameWorksParams.simulator.dt);
dp(1)          = DomainHandlerParams('domainShape','sphere','forceParams',sphereForce,'showDomain',true,'domainWidth',...
                                  25,'dt',frameWorksParams.simulator.dt,'diffusionConst',diffConst,...
                                  'domainCenter',[0 0 0]);
                              
dp(2)          = DomainHandlerParams('domainShape','sphere','forceParams',sphereForce,'showDomain',true,'domainWidth',...
                                  25,'dt',frameWorksParams.simulator.dt,'diffusionConst',diffConst,...
                                  'domainCenter',[2*dp(1).domainWidth+5, 0 0]);
%define a chain;
dh        =DomainHandler(dp);      
initialPt1   =dh.GetRandomBoundarySample(1,1);
initialPt2   =dh.GetRandomBoundarySample(1,2);

chainForce  = ForceManagerParams('springForce',true,'diffusionForce',true,'diffusionConst',diffConst,...
                                 'bendingElasticityForce',false,'bendingConst',0.1,...
                                 'springConst',1,...
                                 'minParticleEqDistance',0);
 
 
 cp(1)= ChainParams('numBeads',64,'initializeInDomain',1,'springForce',true,...
                   'diffusionConst',0.1,...
                   'forceParams',chainForce,...
                   'b', sqrt(frameWorksParams.simulator.dimension),...
                   'dt',frameWorksParams.simulator.dt, 'beadsOnBoundary',[1 64], 'fixedBeadNum',32,'fixedBeadsPosition',initialPt1,...
                   'stickyBeads',[1 64]);


cp(2)= ChainParams('numBeads',64,'b',sqrt(frameWorksParams.simulator.dimension),'initializeInDomain',1,'springForce',true,...
                   'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,...
                   'diffusionConst',0.1, 'beadsOnBoundary',[1 64], 'fixedBeadNum',32,'fixedBeadsPosition',initialPt1,...
                   'stickyBeads',[1 64]);
                    

cp(3)= ChainParams('numBeads',64,'b',sqrt(frameWorksParams.simulator.dimension),'initializeInDomain',2,'springForce',true,...
                   'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,...
                   'diffusionConst',diffConst, 'beadsOnBoundary',[1 64], 'fixedBeadNum',32,'fixedBeadsPosition',initialPt2,...
                   'stickyBeads',[1 64]);

     
               
cp(4)= ChainParams('numBeads',64,'b',sqrt(frameWorksParams.simulator.dimension),'initializeInDomain',2,'springForce',true,...
                   'diffusionConst',0.1,'forceParams',chainForce,'dt',frameWorksParams.simulator.dt,...
                    'beadsOnBoundary',[1 64], 'fixedBeadNum',32,'fixedBeadsPosition',initialPt2,...
                    'stickyBeads',[1 64]);
                   

 
%registre parameters;
frameWorksParams.SetDomainParams(dp);
frameWorksParams.SetChainParams(cp);

% initialize simulator framework
r = RouseSimulatorFramework(frameWorksParams);
%run;
r.Run;


