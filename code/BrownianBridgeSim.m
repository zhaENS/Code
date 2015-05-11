%Create a brownian bridge from beadStart to beadEnd
function paths = BrownianBridgeSim(domainClass,dp,beads,numBeads)


points    = BeadsOnBoundary(dp.domainWidth, dp.dt, dp.diffusionConst,beads);



%initial for a random chain(wiener process)
% for j= 2:numSteps
% 
% w(j,:) = w(j-1,:)+sqrt(2*dp.diffusionConst+dp.dt)*randn(1,3);
% 
% end
% if beads(1)~=numBeads(1)
%     
%          noiseStart  = [sqrt(2*dp.diffusionConst*dp.dt)*randn(beads(1)-1,3)...
%                        ;points(1,:)]
%          path(:,:,1) = cumsum(noiseStart);
% end
% 
%     s = 1;
%     for k = 1:size(points,1)-1
%         path(1,:,s) = points(k,:);
%         for i = 2:beads(1)-1
%             d = 0;
%             while d==0
%     
%             path(i,:,s) = points(k,:)+w(i,:)-i/numSteps*(w(end,:)-points(k+1,:)+points(k,:));
%             d           =  domainClass.InDomain(path(i,:,s));
%    
%             w(1,:)    = randn(1,3);
%             
%                 
%           
%             end
%           
%         end
%     path(numSteps,:,s) = points(k+1,:);
%     s = s+1;
%     end
%     
        
%plot 
for i = 1:length(beads)-1
    B(i) = beads(i+1)-beads(i) ;
    noise = [zeros(1,3);sqrt(2*dp.diffusionConst*dp.dt)...
            *randn(B(i)-1,3)];
    w     = cumsum(noise);
   paths(1,:)   = points(i,:);
  
   paths(2:B(i)+1,:) = repmat(points(i,:),[])+w-...
                              (ones(B(i),1)*(w(end,:)-points(i+1,:)+points(i,:)))...
                              *(ones(B(i),1)*[1:B(i)])'/(B(i));
                              
   paths(B(i)+2,:)   = points(i+1,:); 
   pathsTotal(:,:,i) = paths;
end
[sx, sy, sz]= sphere(15);
sx = sx*dp.domainWidth;
sy = sy*dp.domainWidth;
sz = sz*dp.domainWidth;
daspect([1 1 1]), cameratoolbar
mesh(sx,sy,sz,'FaceColor','none','EdgeColor','m'), hold on
%plot3(noiseStart(:,1),noiseStart(:,2),noiseStart(:,3),'*')

i = 1;

while i <=length(beads)
plot3(points(i,1), points(i,2),points(i,3),'o','MarkerSize',7,'MarkerFaceColor','r');
text(points(i,1)+0.001,points(i,2)-0.01,points(i,3)+0.002,num2str(beads(i)),'FontSize',15);
hold on
% plot3(path(:,1,i),path(:,2,i),path(:,3,i),'Color',[rand rand rand]);
i = i+1;
hold on
% pause(1.0)
end 
c = rand(size(paths,3),3);
for pIdx = 1:size(pathsTotal,3)    
plot3(pathsTotal(:,1,pIdx),pathsTotal(:,2,pIdx),pathsTotal(:,3,pIdx),'Color',c(pIdx,:),'LineWidth',4)
end
end

    
