function HistoneReflectionTest
% This test function moves histones on a Rouse chain
addpath(genpath(fullfile(pwd,'..','PolymerChainDynamics','FrameWork')))
close all
profile on
numSteps = 1000;
%
% % create chain and domain and register them in the ObjectManager

% Initialize simulator framework parameters
simulatorParams = SimulationFrameworkParams('showSimulation',true,'numSteps',1,'dt',0.1);

% create a spherical domain
% assign a force to the domain
sphereForces = ForceManagerParams('lennardJonesForce',false,'diffusionForce',false,'diffusionConst',0.1,...
                                  'LJPotentialWidth',[],'LJPotentialDepth',[],'dt',simulatorParams.simulator.dt);
dp(1)        = DomainHandlerParams('domainShape','sphere','forceParams',sphereForces,...
    'domainWidth',3,'dimension',simulatorParams.simulator.dimension);
chainForces = ForceManagerParams('dt',simulatorParams.simulator.dt,'springForce',false,...
    'bendingElasticityForce',false,'bendingConst',1,'springConst',0.5);
cp          = ChainParams('numBeads',7,'initializeInDomain',1,'forceParams',chainForces);
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
           
histoneParams = HistoneParams('numHistones',5,'forceParams',histoneForce);
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
   
    
Ite = 3;  
if sIdx >1
    %the velocity of each histone;
   v(sIdx,:) = (posHistone(sIdx,:)-posHistone(sIdx-1,:))/simulatorParams.simulator.dt;  
 

           k    = zeros(histoneParams.numHistones,histoneParams.numHistones);
           ite  = 1;
           cols =zeros(1,histoneParams.numHistones);
      timeRem   = simulatorParams.simulator.dt*ones(1,histoneParams.numHistones);
    
    while ite <=Ite 
       for vIdx = 1:histoneParams.numHistones
      
            for kIdx = (vIdx+1):histoneParams.numHistones
                 
           %stock all of time k that each 2 histones are colliding for each
           %row;
           k(vIdx,kIdx)         = (posHistone(sIdx-1,vIdx)-posHistone(sIdx-1,kIdx))/(v(sIdx,kIdx)-v(sIdx,vIdx));
         
            end
   
       
          column                = find(k(vIdx,:)>0 & k(vIdx,:)< timeRem(vIdx));
           
           %detect collision
          if ~isempty(column)
            disp('collision')
            
           %find the min time that 2 histones are colliding in each row;
          
            val   = min(k(vIdx,column(1)));
            cols(vIdx)  = column(1);
            for i =1:length(column)
                if k(vIdx,column(i))<val
                val=k(vIdx,column(i));
                cols(vIdx) = column(i);
                end
            end
           
            
                                                                
            %collision point
            posHistoneCol1  = posHistone(sIdx-1,vIdx)+val*v(sIdx,vIdx);
            posHistoneCol2  = posHistone(sIdx-1,cols(vIdx))+val*v(sIdx,cols(vIdx));
           
            %update the v after they are colliding;
            timeRem(vIdx)      = timeRem(vIdx)-val;
            v(sIdx,vIdx)       = (posHistone(sIdx,vIdx)-posHistoneCol1)/timeRem(vIdx);
            v(sIdx,cols(vIdx)) = (posHistone(sIdx,cols(vIdx))-posHistoneCol2)/timeRem(vIdx);
            
           
       
            %reflection
            %if 2 collision points have velocity opponent,change their direction ;
           
            if sign(v(sIdx,vIdx))~=sign(v(sIdx,cols(vIdx)))
              posHistone(sIdx,vIdx) =  posHistoneCol1-timeRem(vIdx)*v(sIdx,vIdx);
              posHistone(sIdx,cols(vIdx)) =  posHistoneCol2-timeRem(vIdx)*v(sIdx,cols(vIdx));            
             else if abs(v(sIdx,vIdx)) > abs(v(sIdx,cols(vIdx)))
             %if not,change the direction of velocity which is
             %faster than the other
                posHistone(sIdx,vIdx)       =  posHistoneCol1-timeRem(vIdx)*v(sIdx,vIdx);
                posHistone(sIdx,cols(vIdx)) =  posHistoneCol2+timeRem(vIdx)*v(sIdx,cols(vIdx));
                 else
                posHistone(sIdx,vIdx)       = posHistoneCol1+timeRem(vIdx)*v(sIdx,vIdx);
                posHistone(sIdx,cols(vIdx)) = posHistoneCol2-timeRem(vIdx)*v(sIdx,cols(vIdx));
                 end
            end
            
            
            %update the position
          
            posHistone(sIdx-1,vIdx)       = posHistone(sIdx,vIdx);
            posHistone(sIdx-1,cols(vIdx)) = posHistone(sIdx,cols(vIdx));
             
         

            %transforme the coordinate of particles collision to 3D;
               colliIndex = [vIdx cols(vIdx)];
               curPos(colliIndex,:,sIdx) =  (posHistone(sIdx,colliIndex)-edgeLengthTotal(sIdx,curPosVertex1(sIdx,colliIndex)))...
                                       /edgeLength(sIdx,curPosVertex1(sIdx,colliIndex))...
                                       *(chainPos(curPosVertex2(sIdx,colliIndex),:)-chainPos(curPosVertex1(sIdx,colliIndex),:))...
                                       +chainPos(curPosVertex1(sIdx,colliIndex),:);

%              %verify after reflection,the points collisions are always in
%               %the same segement;
%    
%             
% %            
%              Idx1 = max(posHistone(sIdx,vIdx),posHistone(sIdx,cols(vIdx)));
% %             
%              if Idx1==posHistone(sIdx,vIdx)
%                 Idx11 = vIdx;%the index of histone which position is bigger;
%                 Idx12 = cols(vIdx);
%              else
%                 Idx11 = cols(vIdx);
%                 Idx12 = vIdx;
%              end
% %            the case that one point passes the right vertice;            
%             if posHistone(sIdx,Idx11)>edgeLengthTotal(sIdx,curPosVertex2(sIdx,Idx12)) 
%             curPos(Idx11,:,sIdx) = (posHistone(sIdx,Idx11)-edgeLengthTotal(sIdx,curPosVertex1(sIdx,Idx12+1)))...
%                                   /edgeLength(sIdx,curPosVertex1(sIdx,Idx12+1))...
%                                   *(chainPos(curPosVertex2(sIdx,Idx12+1),:)-chainPos(curPosVertex1(sIdx,Idx12+1),:))...
%                                       +chainPos(curPosVertex1(sIdx,Idx12+1),:);
%                 
%             
%             end
%             %the case that one point passes the left vertice;
%            
%             if posHistone(sIdx,Idx12)<edgeLengthTotal(sIdx,curPosVertex1(sIdx,Idx11))
%                  curPos(Idx12,:,sIdx) = (posHistone(sIdx,Idx12)-edgeLengthTotal(sIdx,curPosVertex1(sIdx,Idx11-1)))...
%                                   /edgeLength(sIdx,curPosVertex1(sIdx,Idx11-1))...
%                                   *(chainPos(curPosVertex2(sIdx,Idx11-1),:)-chainPos(curPosVertex1(sIdx,Idx11-1),:))...
%                                       +chainPos(curPosVertex1(sIdx,Idx11-1),:);
%                   
%                   
%             end
            
          end
       
       end
  
       ite = ite+1;
     end
    
end      
    


     
%   update histone graphics    

    if simulatorParams.simulator.showSimulation
          if sIdx >1
         delete(b(1:histoneParams.numHistones));
         end

%         set(histHandle,'XData',h.curPos(:,1),'YData',h.curPos(:,2),'ZData',h.curPos(:,3));
%         text(h.curPos(beads,1)+0.001,h.curPos(beads,2)-0.01,h.curPos(beads,3)+0.002,num2str(beads'),'FontSize',15);
         set(histHandle,'XData',curPos(:,1,sIdx),'YData',curPos(:,2,sIdx),'ZData',curPos(:,3,sIdx));
         text(curPos(beads,1,sIdx)+0.001,curPos(beads,2,sIdx)-0.01,curPos(beads,3,sIdx)+0.002,num2str(beads'),'FontSize',15);
        b = findobj(gca,'Type','Text');
       
        drawnow
    end
end
end
 
