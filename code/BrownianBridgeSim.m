%Create a brownian bridge from beadStart to beadEnd
function path = BrownianBridgeSim(numSteps,domainClass,dp,beads)


points    = BeadsOnBoundary(dp.domainWidth, dp.dt, dp.diffusionConst,beads);

path      = zeros(numSteps,3,length(beads)-1);
w(1,:)    = randn(1,3);

%initial for a random chain(wiener process)
for j= 2:numSteps

w(j,:) = w(j-1,:)+sqrt(2*dp.diffusionConst+dp.dt)*randn(1,3);

end

    s = 1;
    for k = 1:size(points,1)-1
        path(1,:,s) = points(k,:);
        for i = 2:numSteps-1
            d = 0;
            while d==0
    
            path(i,:,s) = points(k,:)+w(i,:)-i/numSteps*(w(end,:)-points(k+1,:)+points(k,:));
            d           =  domainClass.InDomain(path(i,:,s));
   
   
            for j= 2:numSteps

                w(j,:) = w(j-1,:)+sqrt(2*dp.diffusionConst+dp.dt)*randn(1,3);

            end
            
                

            end
          
        end
    path(numSteps,:,s) = points(k+1,:);
    s = s+1;
    end
    
%plot    
daspect([1 1 1]), cameratoolbar
[sx, sy, sz]= sphere(20);
sx = sx*dp.domainWidth;
sy = sy*dp.domainWidth;
sz = sz*dp.domainWidth;

figure, mesh(sx,sy,sz,'FaceColor','none','EdgeColor','g'), hold on
plot3(points(1:2,1), points(1:2,2),points(1:2,3),'o-');
hold on
plot3(path(:,1,1),path(:,2,1),path(:,3,1),'m');
hold on
pause(1.0)
plot3(points(2:3,1), points(2:3,2),points(2:3,3),'o-');
hold on
plot3(path(:,1,2),path(:,2,2),path(:,3,2),'y');
hold on
pause(1.0)
plot3(points(3:4,1), points(3:4,2),points(3:4,3),'o-');
hold on
plot3(path(:,1,3),path(:,2,3),path(:,3,3),'r');
   
       
end

    
