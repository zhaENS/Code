function HistoneReflectionTest
% This test function moves histones on a Rouse chain
addpath(genpath(fullfile(pwd,'..','PolymerChainDynamics','FrameWork')))
close all
profile on
%
% % create chain and domain and register them in the ObjectManager

% Initialize simulator framework parameters
simulatorParams = SimulationFrameworkParams('showSimulation',true,'numSteps',1,'dt',0.1);

% create a spherical domain
% assign a force to the domain
sphereForces = ForceManagerParams('lennardJonesForce',false,'diffusionForce',true,'diffusionConst',0.001,...
                                  'LJPotentialWidth',0.01,'LJPotentialDepth',0.01,'dt',simulatorParams.simulator.dt);
dp(1)        = DomainHandlerParams('domainShape','sphere','forceParams',sphereForces,...
    'domainWidth',3,'dimension',simulatorParams.simulator.dimension);
chainForces = ForceManagerParams('dt',simulatorParams.simulator.dt,'springForce',false,...
    'bendingElasticityForce',false,'bendingConst',1,'springConst',0.05);
cp          = ChainParams('numBeads',4,'initializeInDomain',1,'forceParams',chainForces,'b',1);
% cp(2)     = ChainParams('numBeads',100,'initializeInDomain',1,'forceParams',chainForces);

% register the object parameters in the simulator framework
simulatorParams.SetDomainParams(dp);
simulatorParams.SetChainParams(cp);

% initialize simulator framework
r = RouseSimulatorFramework(simulatorParams);

% get the chain position to initialize the histone on it
[~,initialChainPosition] = r.objectManager.GetMembersPosition(1);
initialChainPosition     = initialChainPosition{1};

% Initialize histones with the chain position
histoneForce = ForceManagerParams('dt',simulatorParams.simulator.dt,'diffusionConst',0.001,...
           'lennardJonesForce',false,'diffusionForce',true,'LJPotentialWidth',0.01,'LJPotentialDepth',0.01);
           
histoneParams = HistoneParams('numHistones',8,'forceParams',histoneForce);
h             = HistoneBis(histoneParams);

h.Initialize(initialChainPosition);
%beads = [1:histoneParams.numHistones];
% get axes
if simulatorParams.simulator.showSimulation
    mAxes = r.simulationGraphics.handles.graphical.mainAxes;
    
    % initialize histone graphics %TODO: incorporate histone graphics in the simulationGraphics class
    histHandle = line('XData',h.curPos(:,1),...
        'YData',h.curPos(:,2),...
        'ZData',h.curPos(:,3),...
        'marker','o',...
        'MarkerFaceColor','y',...
        'MarkerSize',10,...
        'Parent',mAxes,...
        'LineStyle','none');
    daspect([1 1 1])
end
numSteps = [100:200:1000];
t        =[];
for nIdx =1:numel(numSteps)
r.Run;% run initial simulator step
curPosVertex1 = zeros(numSteps(nIdx),histoneParams.numHistones);
curPosVertex2 = zeros(numSteps(nIdx),histoneParams.numHistones);
curPosSlope   = zeros(numSteps(nIdx),histoneParams.numHistones);
edgeLength    = zeros(numSteps(nIdx),cp.numBeads-1);
posHistone    = zeros(numSteps(nIdx),histoneParams.numHistones);
curPos        = zeros(histoneParams.numHistones,3,numSteps(nIdx));
%set the number of times that 2 histones are colliding ;
t(nIdx)       = 1;

 for sIdx = 1:numSteps(nIdx)
        % advance one simulation step
        r.Step;
        [~,chainPos] = r.objectManager.GetMembersPosition(1);
        chainPos     = chainPos{1};
    
        % move the histones
        h.Step(chainPos,simulatorParams.simulator.dt);
    
   

    
     numHistones = histoneParams.numHistones;
    %transform the polygon from 3D to 1D
    for i = 1:cp.numBeads-1
     edgeLength(sIdx,i) = sqrt(sum((chainPos(i+1,:)-chainPos(i,:)).^2));
    end
    edgeLengthTotal = [zeros(numSteps(nIdx),1),cumsum(edgeLength,2)];
   
    %stock the current postion of histones in each
    %step;
    curPos(:,:,sIdx)        =  h.curPos; %the coordinates 3D of histones in each step;
    curPosVertex1(sIdx,:)   =  h.curPosVertex1;%vertices left
    curPosVertex2(sIdx,:)   =  h.curPosVertex2;%vertices right 
    curPosSlope(sIdx,:)     = h.curPosSlope;%position ratio;
    posHistone(sIdx,:)      = edgeLengthTotal(sIdx,curPosVertex1(sIdx,:))+edgeLength(sIdx,curPosVertex1(sIdx,:))...
                            .*curPosSlope(sIdx,:);%the position 1D of histones in each step;
    dt                      = simulatorParams.simulator.dt;
    Ite                     = 3;  

%function detect collision and reflection;    
[curPos,t(nIdx)] = curReflection(numHistones,dt,Ite,curPos,curPosVertex1, curPosVertex2,...
                            chainPos,posHistone,edgeLength,edgeLengthTotal,t(nIdx),sIdx);   
                                           
  if simulatorParams.simulator.showSimulation

         set(histHandle,'XData',curPos(:,1,sIdx),'YData',curPos(:,2,sIdx),'ZData',curPos(:,3,sIdx));
      
       
        drawnow
  end
 end
end
figure
[h,x]=hist(t-1,numSteps);
bar(x,h);
end

