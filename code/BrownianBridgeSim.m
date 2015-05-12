%Create a brownian bridge from beadStart to beadEnd
function BrownianBridgeSim(domainClass,dp,beads,numBeads)

%obtain the coordinates of the points on the sphere.
points    = BeadsOnBoundary(dp.domainWidth, dp.dt, dp.diffusionConst,beads);

for i = 1:length(beads)-1
    B(i) = beads(i+1)-beads(i) ;
    paths(1,:)   = points(i,:);
    d            = 0;
    while min(d) == 0
    noise = [zeros(1,3);sqrt(2*dp.diffusionConst*dp.dt)...
            *randn(B(i)-1,3)];
    w     = cumsum(noise);%initialisation of a random process.
   
  
   paths(2:B(i)+1,:) = repmat(points(i,:),[B(i),1])+w-...
                              (repmat(w(end,:),[B(i),1])-repmat(points(i+1,:),[B(i),1])+repmat(points(i,:),[B(i),1]))...
                              .*(repmat([1:B(i)],[3,1]))'/(B(i));
   d                 = domainClass.InDomain(paths(2:B(i)+1,:));
    end
   paths(B(i)+2,:)   = points(i+1,:);
  plot3(paths(:,1),paths(:,2),paths(:,3),'Color',[rand rand rand],'LineWidth',4);
  hold on
end

%constructor the rest of the chain
if beads(1)~=1 
       d1 = 0;
       
      while  min(d1) == 0 
       
        pathsReste1 = cumsum([points(1,:);-sqrt(2*dp.diffusionConst*dp.dt)...
            *randn(beads(1)-1-1,3)])
        
        d1          =  domainClass.InDomain(pathsReste1);
      end
     plot3(pathsReste1(:,1),pathsReste1(:,2),pathsReste1(:,3),'Color',[rand rand rand],'LineWidth',4);
    points(end+1,:) = pathsReste1(end,:);
    beads(end+1)    = 1;
end
if beads(end)< numBeads 
       d2 = 0;
       
      while  min(d2) == 0 
       
         pathsReste2 = cumsum([points(end,:);sqrt(2*dp.diffusionConst*dp.dt)...   
    *randn(numBeads-beads(end),3)]);
        d2           =  domainClass.InDomain(pathsReste2);
      end
    plot3(pathsReste2(:,1),pathsReste2(:,2),pathsReste2(:,3),'Color',[rand rand rand],'LineWidth',4);
    points(end+2,:) = pathsReste2(end,:);
    beads(end+2)    = numBeads;
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
i = i+1;
hold on
% pause(1.0)
end 
% c = rand(size(paths,3),3);
% for pIdx = 1:size(pathsTotal,3)    
% plot3(pathsTotal(:,1,pIdx),pathsTotal(:,2,pIdx),pathsTotal(:,3,pIdx),'Color',c(pIdx,:),'LineWidth',4)
% end
end

    
