function HistoneReflectionTest
% This test function moves histones on a Rouse chain
addpath(genpath(fullfile(pwd,'..','PolymerChainDynamics','FrameWork')))
close all
profile on
numSteps = 1000;
%
% % create chain and domain and register them in the ObjectManager

% Initialize simulator framework parameters
simulatorParams = SimulationFrameworkParams('showSimulation',true,'numSteps',1,'dt',0.01);

% create a spherical domain
% assign a force to the domain
sphereForces = ForceManagerParams('lennardJonesForce',false,'diffusionForce',false,'diffusionConst',0.1,...
                                  'LJPotentialWidth',[],'LJPotentialDepth',[],'dt',simulatorParams.simulator.dt);
dp(1)        = DomainHandlerParams('domainShape','sphere','forceParams',sphereForces,...
    'domainWidth',3,'dimension',simulatorParams.simulator.dimension);
chainForces = ForceManagerParams('dt',simulatorParams.simulator.dt,'springForce',false,...
    'bendingElasticityForce',false,'bendingConst',1,'springConst',0.5);
cp          = ChainParams('numBeads',3,'initializeInDomain',1,'forceParams',chainForces);
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
histoneForce = ForceManagerParams('dt',simulatorParams.simulator.dt,'diffusionConst',0.1,...
           'lennardJonesForce',false,'diffusionForce',true,'LJPotentialWidth',0.01,'LJPotentialDepth',0.01);
           
histoneParams = HistoneParams('numHistones',2,'forceParams',histoneForce);
h             = HistoneBis(histoneParams);

h.Initialize(initialChainPosition);
beads = [1:histoneParams.numHistones];
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
r.Run;% run initial simulator step
curPosVertex1 = zeros(numSteps,histoneParams.numHistones);
curPosVertex2 = zeros(numSteps,histoneParams.numHistones);
curPosSlope   = zeros(numSteps,histoneParams.numHistones);
edgeLength    = zeros(numSteps,cp.numBeads-1);
for sIdx = 1:numSteps
    % advance one simulation step
    r.Step;
    [~,chainPos] = r.objectManager.GetMembersPosition(1);
    chainPos     = chainPos{1};
    
    % move the histones
    h.Step(chainPos,simulatorParams.simulator.dt);
    
    %transform the polygon from 3D to 1D
    for i = 1:cp.numBeads-1
     edgeLength(sIdx,i) = sqrt(sum((chainPos(i+1,:)-chainPos(i,:)).^2));
    end
    edgeLengthTotal = [zeros(numSteps,1), cumsum(edgeLength,2)];
    
    %stock the previous position and current postion of histones in each
    %step;
    curPosition (:,:,sIdx)  = h.curPos;%coordinate in 3D;
    curPosVertex1(sIdx,:)   =  h.curPosVertex1;
    curPosVertex2(sIdx,:)   =  h.curPosVertex2;
    curPosSlope(sIdx,:) = h.curPosSlope;%position ratio;
    if sIdx >= 2
      for vIdx = 1:histoneParams.numHistones
          %the velocity of each histone during deltaT;
        v(sIdx,vIdx) = norm(curPosition(vIdx,:,sIdx)-curPosition(vIdx,:,sIdx-1))/simulatorParams.simulator.dt;
      end
   
    %detect collision
    cIdx = 1;
    k = zeros(histoneParams.numHistones);
   for vIdx = 1:histoneParams.numHistones
       for kIdx = (vIdx+1):histoneParams.numHistones
           if min(curPosVertex1(sIdx-1,vIdx)-1,curPosVertex1(sIdx-1,kIdx)-1)~=0
         Dist   = edgeLengthTotal(sIdx-1,max(curPosVertex1(sIdx-1,vIdx),curPosVertex1(sIdx-1,kIdx)))...
                 -edgeLengthTotal(sIdx-1,min(curPosVertex1(sIdx-1,vIdx)-1,curPosVertex1(sIdx-1,kIdx)-1));
           else
         Dist   = edgeLengthTotal(sIdx-1,max(curPosVertex1(sIdx-1,vIdx),curPosVertex1(sIdx-1,kIdx)));      
           end
         if curPosVertex1(sIdx-1,vIdx)>curPosVertex1(sIdx-1,kIdx)
             k(vIdx,kIdx) = (Dist-edgeLength(sIdx-1,kIdx)*curPosSlope(sIdx-1,kIdx)-edgeLength(sIdx-1,vIdx)*...
             (1-curPosSlope(sIdx-1,vIdx)))/(v(sIdx,vIdx)-v(sIdx,kIdx));
         else
             k(vIdx,kIdx) = (Dist-edgeLength(sIdx-1,vIdx)*curPosSlope(sIdx-1,vIdx)-edgeLength(sIdx-1,kIdx)*...
             (1-curPosSlope(sIdx-1,kIdx)))/(v(sIdx,vIdx)-v(sIdx,kIdx));

         end
   %k=(x2(t)-x1(t))/(v1(t+deltaT)-v2(t+deltaT))
       if (k(vIdx,kIdx) < simulatorParams.simulator.dt) &&(k(vIdx,kIdx)>0)
            sprintf('point %s and point %s are colliding ',num2str(vIdx),num2str(kIdx))
        end
        cIdx = cIdx+1; 
       end
   end
    end
% find colliding histones
        if sIdx >1
            [hist1,hist2] = find(k>0 & k<=simulatorParams.simulator.dt);
            if ~isempty(hist1)
                disp('collision')
            end
        end
        
%   update histone graphics    
    if simulatorParams.simulator.showSimulation
          if sIdx >1
         delete(b(1:histoneParams.numHistones));
         end

        set(histHandle,'XData',h.curPos(:,1),'YData',h.curPos(:,2),'ZData',h.curPos(:,3));
        text(h.curPos(beads,1)+0.001,h.curPos(beads,2)-0.01,h.curPos(beads,3)+0.002,num2str(beads'),'FontSize',15);
      
       
        drawnow
    end
  b = findobj(gca,'Type','Text');
end
 
end