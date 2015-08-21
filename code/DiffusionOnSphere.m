%This function is used to do the diffusion on sphere
%initialPoint : the initial point which is on the sphere
%dt :time step
%numSteps :number of steps that walks on the sphere
%domainCenter : the center of mass of sphere
%diffConst : a diffusion const for the processus 
%radius : the radius of the sphere
function paths = DiffusionOnSphere(initialPoint,dt,diffConst,numSteps,domainCenter,radius,w)
%check initialPoint is on the sphere or not
rho = sqrt(sum((initialPoint-domainCenter).^2,2));
if (rho-radius)^2>eps
    error('the initial point is not on the sphere');
end
%generate a path 2D with numStep;
if w==0
noise = [zeros(1,2);sqrt(2*diffConst*dt)*randn(numSteps-1,2)];
else
noise=[zeros(1,2);sqrt(2*diffConst*dt)*ones(numSteps-1,2)*cos(w*dt)];    
end
paths =  cumsum(noise);

%find the phi and theta on the sphere;
phi = atan(initialPoint(2)/initialPoint(1));
if phi<0
phi = phi+2*pi;
end
theta   = acos(initialPoint(3)./radius);
if initialPoint(1)<0
theta = -theta;
end
%rotate the paths to the initialPoint ;
paths = [paths zeros(numSteps,1)];

paths(:,3) = bsxfun(@plus, paths(:,3),radius);
% plot3([domainCenter,initialPoint(1,1)],[domainCenter initialPoint(1,2)],[domainCenter initialPoint(1,3)],'g','LineWidth',3);hold on
% plot3([paths(:,1)],[paths(:,2)],[paths(:,3)],'g','LineWidth',3);

%hold on
%defini the rotation matrix; 
% plot3(initialPoint(:,1),initialPoint(:,2),initialPoint(:,3),'o','MarkerSize',4);
% hold on;
Ry = [cos(theta) 0 sin(theta) ;0 1 0 ;-sin(theta) 0 cos(theta)];
Rz = [cos(phi) -sin(phi) 0;sin(phi) cos(phi) 0;0 0 1];

paths = (Rz*Ry*bsxfun(@minus,paths,domainCenter)')';
%paths = bsxfun(@plus, paths,initialPoint);
%paths(1,:) =initialPoint;
% c=[rand rand rand];
% [sx,sy,sz]  = sphere(20);
% sx          = sx*radius;
% sy          = sy*radius;
% sz          = sz*radius;
% mesh(sx,sy,sz,'FaceColor','none','EdgeColor','g'), hold on, 
% cameratoolbar
% daspect([1 1 1])
% plot3(paths(:,1), paths(:,2),paths(:,3),'LineWidth',2);
% 

%project the paths back on the sphere;
distToCenter = pdist2mex(paths(2:end,:)',domainCenter','euc',[],[],[]);
t            = radius./distToCenter;
paths(2:end,:) = bsxfun(@plus, domainCenter,bsxfun(@times,t,bsxfun(@minus,paths(2:end,:),domainCenter)));

% hold on ,
% plot3(paths(:,1),paths(:,2),paths(:,3),'LineWidth',2,'color',c);
% xlabel('x'), ylabel('y'), zlabel('z')
% plot3([domainCenter(1) radius],[domainCenter(2) domainCenter(2)],[domainCenter(3) domainCenter(3)])
% plot3([domainCenter(1) domainCenter(1)],[domainCenter(2) radius],[domainCenter(3) domainCenter(3)])
% plot3([domainCenter(1) domainCenter(1)],[domainCenter(2) domainCenter(2)],[domainCenter(3) radius])

end