  close all
  clear all
% define rotation matrix
Rx =@(roll) [1 0 0;0 cos(roll) -sin(roll);0 sin(roll) cos(roll)];
Ry =@(pitch) [cos(pitch) 0 sin(pitch) ;0 1 0 ;-sin(pitch) 0 cos(pitch)];
Rz =@(yaw) [cos(yaw) -sin(yaw) 0;sin(yaw) cos(yaw) 0;0 0 1];

rotationMatrix = @(yaw,pitch,roll) Rz(yaw)*Ry(pitch)*Rx(roll);



%initialPoint = [sqrt(50) 5 5];
initialPoint = [sqrt(50)+1 6 5];
dt           = 0.1;
diffConst    = 1;
numSteps     = 100;
domainCenter = [1 1 0];
radius       = 10;
w            = 0;
paths        = DiffusionOnSphere(initialPoint,dt,diffConst,numSteps,domainCenter,radius,w);
theta = atan(paths(:,2)./paths(:,1));
phi   = atan(sqrt(paths(:,1).^2+paths(:,2).^2)./paths(:,3));

% define translation and rotation 
yaw   = pi;
pitch = 0;
roll  = 0;

TransVector = [0 0 0];

[sx,sy,sz]    = sphere(20);
sx            = sx*radius;
sy            = sy*radius;
sz            = sz*radius;
cameratoolbar;
daspect([1 1 1]);
    plot3(paths(:,1),paths(:,2),paths(:,3),'LineWidth',3);
    hold on
    mesh(sx,sy,sz,'facecolor','none','EdgeColor','g')
   
%figure,
cameratoolbar;
daspect([1 1 1]);
%alpha = 45;
     %apply the velocity for each coordinate;
     sx = sx+TransVector(1);
     sy = sy+TransVector(2);
     sz = sz+TransVector(3);
     sh = mesh(sx,sy,sz,'facecolor','none','EdgeColor','g');
%      for rIdx = 1:size(sx,1);
%          for cIdx = 1:size(sx,2)
%              temp = rotationMatrix(phisx(rIdx,cIdx),       
%          end
%      end

  
      %rotate 
%      paths(:,1) = paths(:,1)+TransVector(1);
%      paths(:,2) = (radius.*sin(theta).*sin(phi)).*cos(alpha*pi/180)...
%                   -(radius.*cos(theta)).*sin(alpha*pi/180)+TransVector(2);
%      paths(:,3) =(radius.*sin(theta).*sin(phi)).*sin(alpha*pi/180)...
%                   +(radius.*cos(theta)).*cos(alpha*pi/180)+TransVector(3);     


pause(0.5);
for pIdx = 1:size(paths,1)
 paths(pIdx,1:3) = (rotationMatrix(yaw,pitch,roll)*(paths(pIdx,:)-domainCenter)')' +(domainCenter +TransVector);
end
plot3(paths(:,1),paths(:,2),paths(:,3),'LineWidth',3,'Color','r'); 
hold on
line('XData',[domainCenter(1) domainCenter(1)+radius],'YData',[0 0],'Color','r','LineWidth',4)
text(domainCenter(1)+radius,0,0,'X','FontSize',20);
line('XData',[0 0],'YData',[domainCenter(2) domainCenter(2)+radius],'Color','g','LineWidth',4)
text(0,domainCenter(2)+radius,0,'Y','FontSize',20);
line('XData',[0 0],'YData',[0 0],'ZData',[domainCenter(3) domainCenter(3)+radius],'Color','b','LineWidth',4)
text(0,0,domainCenter(3)+radius,'Z','FontSize',20);