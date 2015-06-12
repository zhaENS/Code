%Create a brownian bridge from beadStart to beadEnd
function chainPath = BrownianBridgeSim(initialPoint,domainClass,dp,beads,numBeads)
tic
beads     = sort(beads);
%obtain the coordinates of the points on the sphere.
points    = BeadsOnBoundary(initialPoint,dp.domainWidth, dp.dt, dp.diffusionConst,beads);
chainPath = points(1,:);
% initialize the BrownianBridge class 
%noiseSTD = sqrt(2*dp.dt*dp.diffusionConst);
%bb = BrownianBridge('realizations',1,'dimension',3,'constructionType','normal','noiseSTD',noiseSTD);
% for bIdx = 1:(numel(beads)-1)
%     % build a Brownian bridge between any two points on the boundary
%     numPoints            = (beads(bIdx+1)-beads(bIdx)+1);
%     startPoint           = points(bIdx,:);
%     endPoint             = points(bIdx+1,:);
%     bb.params.numPoints  = numPoints;
%     bb.params.startPoint = startPoint;
%     bb.params.endPoint   = endPoint;
%     bb.GetBridge;
%     tempPath = bb.paths{1};
%     % check that the points are inside the domain 
%         inDomain = domainClass.InDomain(tempPath);
%         f        = find(~inDomain);
%         for fIdx = 1:numel(f)
%             pInDomain = false;
%             while ~pInDomain 
%                 tempPath(f(fIdx),:) = tempPath(f(fIdx)-1,:)+randn(1,dp.dimension).*noiseSTD;
%                 pInDomain = domainClass.InDomain(tempPath(f(fIdx),:));
%             end
%         end               
%     n = size(chainPath,1);
%     chainPath(n+1:n+numPoints-1,:)= tempPath(2:end,:);
%     bb.Reset;    
% end
%     
% % link the the ends of the chain 
% numPoints = beads(1);    
% startPath = [];
% if beads(1)~=1
%    
%     pos1      = chainPath(1,:);
%     for bIdx = 1:numPoints        
%         inDomain = false;
%         while ~inDomain 
%              pos2 = pos1 + randn(1,dp.dimension).*noiseSTD;
%              inDomain = domainClass.InDomain(pos2);            
%         end
%         startPath(bIdx,:) = pos1;
%         pos1 = pos2;
%         
%     end
%     startPath = flipud(startPath(2:end,:)); % omit the first point 
% end
% 
% % connect the end
% numPoints = numBeads-beads(end)+1;  
% endPath   = [];
% if beads(end)~=numBeads
%    
%     pos1     = chainPath(1,:);
%     for bIdx = 1:numPoints        
%         inDomain = false;
%         while ~inDomain 
%              pos2 = pos1 + randn(1,dp.dimension).*noiseSTD;
%              inDomain = domainClass.InDomain(pos2);            
%         end
%         endPath(bIdx,:) = pos1;
%         pos1 = pos2;
%         
%     end
%     endPath = (endPath(2:end,:)); % omit the first point 
% end
% 
% chainPath = [startPath; chainPath; endPath];
% toc
% figure;
% [sx, sy, sz]= sphere(15);
% sx = sx*dp.domainWidth;
% sy = sy*dp.domainWidth;
% sz = sz*dp.domainWidth;
% daspect([1 1 1]), cameratoolbar
% mesh(sx,sy,sz,'FaceColor','none','EdgeColor','m'), hold on
% plot3(chainPath(:,1),chainPath(:,2),chainPath(:,3),'Color',[rand rand rand],'LineWidth',4);
% 
% plot3(chainPath(beads,1), chainPath(beads,2),chainPath(beads,3),'o','MarkerSize',7,'MarkerFaceColor','r','LineStyle','none');
% % text(chainPath(beads,1)+0.001,chainPath(beads,2)-0.01,chainPath(beads,3)+0.002,num2str(beads'),'FontSize',15);


tic
pathsTotal = [];
for i = 1:length(beads)-1
    paths = [];
    % build a Brownian bridge between any two points on the boundary      
    B(i) = beads(i+1)-beads(i) ;
    paths(1,:)   = points(i,:);
    d            = 0;
    while min(d) == 0
    noise = [zeros(1,3);sqrt(2*dp.diffusionConst*dp.dt)...
            *randn(B(i)-1,3)];
    w     = cumsum(noise);%initialisation of a random process.
   
   x    = repmat(points(i,:),[B(i),1]);
   y    = repmat(points(i+1,:),[B(i),1]);
   t    = repmat([1:B(i)],[3,1])';
   T    = B(i);
   wT   = repmat(w(end,:),[B(i),1]);
   paths(2:B(i)+1,:) = x+w-(wT-y+x).*(t/(T));
   d                 = domainClass.InDomain(paths(2:B(i)+1,:));
    end
%   plot3(paths(:,1),paths(:,2),paths(:,3),'Color',[rand rand rand],'LineWidth',4);
%   hold on
  pathsTotal = [pathsTotal;paths(1:B(i),:)];
end
pathsReste1 = [];
pathsReste2 = [];
%constructor the rest of the chain
if beads(1)~= 1
     d1 = 0;
      while  min(d1) == 0
       
        pathsReste1 = cumsum([points(1,:);sqrt(2*dp.diffusionConst*dp.dt)...
            *randn(beads(1)-1,3)]);
        d1          =  domainClass.InDomain(pathsReste1);
      end
    plot3(pathsReste1(:,1),pathsReste1(:,2),pathsReste1(:,3),'Color',[rand rand rand],'LineWidth',4);
    hold on
end
if beads(end)== numBeads
    pathsTotal = [pathsTotal;points(end,:)];
else 
        d2 = 0;
         while min(d2) == 0
        pathsReste2 = cumsum([points(end,:);sqrt(2*dp.diffusionConst*dp.dt)...   
    *randn(numBeads-beads(end),3)]);
        d2          =  domainClass.InDomain(pathsReste2);
         end
       % plot3(pathsReste2(:,1),pathsReste2(:,2),pathsReste2(:,3),'Color',[rand rand rand],'LineWidth',4);
        %hold on
end


chainPath = [pathsTotal;pathsReste1(2:end,:);pathsReste2];
toc

% [sx, sy, sz]= sphere(20);
% sx = sx*dp.domainWidth;
% sy = sy*dp.domainWidth;
% sz = sz*dp.domainWidth;
% daspect([1 1 1]), cameratoolbar
% mesh(sx,sy,sz,'FaceColor','none','EdgeColor','m'), hold on

i = 1;

% while i <=length(beads)
% plot3(points(i,1), points(i,2),points(i,3),'o','MarkerSize',7,'MarkerFaceColor','r');
% text(points(i,1)+0.001,points(i,2)-0.01,points(i,3)+0.002,num2str(beads(i)),'FontSize',20);
% hold on
% i = i+1;
% hold on
% % pause(1.0)
% end 

end

    
