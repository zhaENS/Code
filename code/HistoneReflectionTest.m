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
histoneForce = ForceManagerParams('dt',simulatorParams.simulator.dt,'diffusionConst',0.01,...
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
posHistone    = zeros(numSteps,histoneParams.numHistones);
curPos        = zeros(histoneParams.numHistones,3,numSteps);
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
  edgeLengthTotal = [zeros(numSteps,1),cumsum(edgeLength,2)];
    
    %stock the previous position and current postion of histones in each
    %step;
    curPos(:,:,sIdx)        =h.curPos; %the coordinates 3D of histones in each step;
    curPosVertex1(sIdx,:)   =  h.curPosVertex1;%vertices left
    curPosVertex2(sIdx,:)   =  h.curPosVertex2;%vertices right 
    curPosSlope(sIdx,:)     = h.curPosSlope;%position ratio;
    posHistone(sIdx,:)      = edgeLengthTotal(sIdx,curPosVertex1(sIdx,:))+edgeLength(sIdx,curPosVertex1(sIdx,:))...
                            .*curPosSlope(sIdx,:);%the position 1D of histones in each step;
    if sIdx >1
    %the velocity of each histone during deltaT;
        v(sIdx,:) = (posHistone(sIdx,:)-posHistone(sIdx-1,:))/simulatorParams.simulator.dt;
    %detect collision
    
   
   for vIdx = 1:histoneParams.numHistones
       for kIdx = (vIdx+1):histoneParams.numHistones
           
           k(vIdx,kIdx)  = (posHistone(sIdx-1,vIdx)-posHistone(sIdx-1,kIdx))/(v(sIdx,vIdx)-v(sIdx,kIdx));
               % find colliding histones
            [hist1,hist2] = find(k(vIdx,kIdx)>0 & k(vIdx,kIdx)<=simulatorParams.simulator.dt);
           if ~isempty(hist1)
                disp('collision')
            %collision point
            posHistone(sIdx,vIdx) = posHistone(sIdx-1,vIdx)+k(vIdx,kIdx)*v(sIdx,vIdx);
            posHistone(sIdx,kIdx) = posHistone(sIdx-1,kIdx)+k(vIdx,kIdx)*v(sIdx,kIdx);
             %reflection
             timeRest = simulatorParams.simulator.dt-k(vIdx,kIdx);
             if sign(v(sIdx,vIdx))~=sign(v(sIdx,kIdx))
                 %if 2 collision points have velocity opponent,change their
                 %direction ;
                posHistone(sIdx,vIdx) = posHistone(sIdx,vIdx)-timeRest*v(sIdx,vIdx);
                posHistone(sIdx,kIdx) = posHistone(sIdx,kIdx)-timeRest*v(sIdx,kIdx);            
             else if abs(v(sIdx,vIdx)) > abs(v(sIdx,kIdx))
                 %if not,change the direction of velocity which is
                 %faster than the other
                posHistone(sIdx,vIdx) = posHistone(sIdx,vIdx)-timeRest*v(sIdx,vIdx);
                posHistone(sIdx,kIdx) = posHistone(sIdx,kIdx)+timeRest*v(sIdx,kIdx);
                 else
                posHistone(sIdx,vIdx) = posHistone(sIdx,vIdx)+timeRest*v(sIdx,vIdx);
                posHistone(sIdx,kIdx) = posHistone(sIdx,kIdx)-timeRest*v(sIdx,kIdx);
                 end
             end
             %transforme the coordinate of particles collision to 3D;
             colliIndex = [vIdx kIdx];
            
              curPos(colliIndex,:,sIdx) =  (posHistone(sIdx,colliIndex)-edgeLengthTotal(sIdx,curPosVertex1(sIdx,colliIndex)))...
                                      /edgeLength(sIdx,curPosVertex1(sIdx,colliIndex))...
                                      *(chainPos(curPosVertex2(sIdx,colliIndex),:)-chainPos(curPosVertex1(sIdx,colliIndex),:))...
                                      +chainPos(curPosVertex1(sIdx,colliIndex),:);
              %verify after reflection,the points collisions are always in
              %the same segement;
              
              %the case that one point passes the right vertice;
            Idx1 = size(max(posHistone(sIdx,vIdx),posHistone(sIdx,kIdx)),2);
            Idx2 = size(min(posHistone(sIdx,vIdx),posHistone(sIdx,kIdx)),2);
            if posHistone(sIdx,Idx1)>edgeLengthTotal(sIdx,curPosVertex2(sIdx,Idx2)) 
            curPos(Idx1,:,sIdx) = (posHistone(sIdx,Idx1)-edgeLengthTotal(sIdx,curPosVertex1(sIdx,Idx1+1)))...
                                  /edgeLength(sIdx,curPosVertex1(sIdx,Idx1+1))...
                                  *(chainPos(curPosVertex2(sIdx,Idx1+1),:)-chainPos(curPosVertex1(sIdx,Idx1+1),:))...
                                      +chainPos(curPosVertex1(sIdx,Idx1+1),:);
                
            
            end
            %the case that one point passes the left vertice;
           
            if posHistone(sIdx,Idx2)<edgeLengthTotal(sIdx,curPosVertex1(sIdx,Idx1))
                 curPos(Idx2,:,sIdx) = (posHistone(sIdx,Idx2)-edgeLengthTotal(sIdx,curPosVertex1(sIdx,Idx2)))...
                                  /edgeLength(sIdx,curPosVertex1(sIdx,max(Idx2-1,1)))...
                                  *(chainPos(curPosVertex2(sIdx,Idx2),:)-chainPos(curPosVertex1(sIdx,Idx2),:))...
                                      +chainPos(curPosVertex1(sIdx,Idx2),:);
                  
                  
            end
           end
           
       end
   end
    end
        
%   update histone graphics    

    if simulatorParams.simulator.showSimulation
          if sIdx >1
         delete(b(1:histoneParams.numHistones));
         end

       % set(histHandle,'XData',h.curPos(:,1),'YData',h.curPos(:,2),'ZData',h.curPos(:,3));
        %text(h.curPos(beads,1)+0.001,h.curPos(beads,2)-0.01,h.curPos(beads,3)+0.002,num2str(beads'),'FontSize',15);
        set(histHandle,'XData',curPos(:,1,sIdx),'YData',curPos(:,2,sIdx),'ZData',curPos(:,3,sIdx));
        text(curPos(beads,1,sIdx)+0.001,curPos(beads,2,sIdx)-0.01,curPos(beads,3,sIdx)+0.002,num2str(beads'),'FontSize',15);
        b = findobj(gca,'Type','Text');
       
        drawnow
    end
end
end
 
