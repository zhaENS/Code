%Create a brownian bridge from beadStart to beadEnd
function chainPath = BrownianBridgeSim(initialPoint,domainClass,dp,beadsOnBoundary,numBeads)
% tic
beadsOnBoundary     = sort(beadsOnBoundary);
% obtain the coordinates of the points on the sphere.
domainCenter = domainClass.params.domainCenter;
points       = BeadsOnBoundary(initialPoint,dp.domainWidth,domainCenter,dp.dt, dp.diffusionConst,beadsOnBoundary);


pathsTotal = [];
for i = 1:length(beadsOnBoundary)-1
    %if ~ismember(i,beadsOnBoundary)
        
        % build a Brownian bridge between any two points on the boundary
        B(i)       = beadsOnBoundary(i+1)-beadsOnBoundary(i);
        paths      = zeros(B(i)+1,3);
        paths(1,:) = points(i,:);
        d          = 0;
        while min(d) == 0
            noise = [zeros(1,3);sqrt(2*dp.diffusionConst*dp.dt)...
                *randn(B(i)-1,3)];
            w     = cumsum(noise);%initialisation of a random process.
            
            x     = repmat(points(i,:),[B(i),1]);
            y     = repmat(points(i+1,:),[B(i),1]);
            t     = repmat([1:B(i)],[3,1])';
            T     = B(i);
            wT    = repmat(w(end,:),[B(i),1]);
            %the path of chain from current (beadsOnBoundary+1) to next
            %beadsOnBoundary;
            paths(2:end,:) = x+w-(wT-y+x).*(t/(T));
        
            %correct the paths(end,:) which is calculate by formula to the correct position
            %of beads on the boundary ;
            paths(end,:)   = points(i+1,:);
            [d,onDomain]   = domainClass.InDomain(paths(2:end-1,:));
        end
        
    
      
        %add the paths to the total paths;
        pathsTotal   = [pathsTotal;paths(2:end-1,:)];

   
        %   plot3(paths(:,1),paths(:,2),paths(:,3),'Color',[rand rand rand],'LineWidth',4);
        %   hold on
%         pathsTotal = [pathsTotal;paths(1:B(i),:)];

        % end
end
%add the beads on the boundary to the pathsTotal;
pathsTotal                    =[pathsTotal;points];

%change the position of beads to the right order;
pointsTempo                   =zeros(numel(beadsOnBoundary),3);
pointsTempo                   = pathsTotal(beadsOnBoundary,:);
pathsTotal(beadsOnBoundary,:) = points;
pathsTotal(numBeads-numel(beadsOnBoundary)+1:numBeads,:) = pointsTempo;
pathsReste1 = [];
pathsReste2 = [];
%constructor the rest of the chain
if beadsOnBoundary(1)~= 1
     d1 = 0;
      while  min(d1) == 0
       
        pathsReste1  = cumsum([points(1,:);sqrt(2*dp.diffusionConst*dp.dt)...
            *randn(beadsOnBoundary(1)-1,3)]);
       [d1,onDomain] =  domainClass.InDomain(pathsReste1(2:end,:));
      end
%     plot3(pathsReste1(:,1),pathsReste1(:,2),pathsReste1(:,3),'Color',[rand rand rand],'LineWidth',4);
%     hold on
end
if beadsOnBoundary(end)~= numBeads
   d2 = 0;
         while min(d2) == 0
        pathsReste2 = cumsum([points(end,:);sqrt(2*dp.diffusionConst*dp.dt)...   
    *randn(numBeads-beadsOnBoundary(end),3)]);
        [d2,onDomain]         =  domainClass.InDomain(pathsReste2(2:end,:));
         end
       % plot3(pathsReste2(:,1),pathsReste2(:,2),pathsReste2(:,3),'Color',[rand rand rand],'LineWidth',4);
        %hold on
end


chainPath = [pathsTotal;pathsReste1(2:end,:);pathsReste2(2:end,:)];
% toc

% [sx, sy, sz]= sphere(20);
% sx = sx*dp.domainWidth;
% sy = sy*dp.domainWidth;
% sz = sz*dp.domainWidth;
% daspect([1 1 1]), cameratoolbar
% mesh(sx,sy,sz,'FaceColor','none','EdgeColor','m'), hold on


% while i <=length(beads)
% plot3(points(i,1), points(i,2),points(i,3),'o','MarkerSize',7,'MarkerFaceColor','r');
% text(points(i,1)+0.001,points(i,2)-0.01,points(i,3)+0.002,num2str(beads(i)),'FontSize',20);
% hold on
% i = i+1;
% hold on
% % pause(1.0)
% end 

end

    
