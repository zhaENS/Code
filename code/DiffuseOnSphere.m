function points = DiffuseOnSphere(initialPoint,numSteps,radius, dt, diffusionConst)
% random walk on a sphere assuming the directions are indipendent
% initialPoint should be a point on the suface 
% numSteps is the number of steps to advance from initialPoint
% radius is the radius of the sphere
% dt is the time step 
% diffusionConst is the diffusion constant for the process

% close all
% numSteps       = 1000;
% radius         = 2;
% dt             = 0.001;
% diffusionConst = 0.1;


% translate initial position to spherical cooridinates 
rho = sqrt(sum(initialPoint.^2));
if (rho-radius)^2 >eps
    error('the initial point is not on the sphere')
end
theta = atan(initialPoint(2)/initialPoint(1));
phi   = atan(sqrt(initialPoint(1).^2 + initialPoint(2).^2)./initialPoint(3));

% % % random angle sampling for points on the sphere
% u     = rand(1);
% v     = rand(1);
% phi   = 2*pi*u;
% theta = acos(2*v -1);

w          = sqrt(2*diffusionConst*dt);%./(radius);
mu         = pi/2;
dAngles    = RandomWrappedNormalOnCircle(mu,(1/w),[numSteps-1,2]);
dTheta     = (dAngles(:,1)-mu);
dPhi       = (dAngles(:,2)-mu);

theta      = (cumsum([theta;dTheta]));
phi        = (cumsum([phi;(dPhi)]));


x           = radius*sin(phi).*cos(theta);
y           = radius*sin(phi).*sin(theta);
z           = radius*cos(phi);
points      = [x,y,z];

% --- Plot ---

% [sx,sy,sz]  = sphere(20);
% sx          = sx*radius;
% sy          = sy*radius;
% sz          = sz*radius;
% 
% figure, subplot(141),plot(phi),set(gca,'YLim',[-2*pi 2*pi]), title('phi'), 
%         subplot(142), plot(theta),set(gca,'YLim',[-2*pi 2*pi]), title('theta')
%         x = xcorr(theta,phi);
%         subplot(143), plot(x./sum(x))
%         dx = xcorr(dTheta,dPhi);
%         subplot(144), plot(dx./sum(dx)), title('xcorr dTheta, dPhi')
% 
% figure, plot(theta,phi), xlabel('theta'), ylabel('phi'), set(gca,'XLIm',2*[-pi pi],'YLim',2*[-pi pi])
% figure, mesh(sx,sy,sz,'FaceColor','none','EdgeColor','g'), hold on, 
% plot3(points(:,1), points(:,2),points(:,3))
% line('XData',[0 radius*sin(phi(1)).*cos(theta(1))],...
%      'YData',[0 radius*sin(phi(1)).*sin(theta(1))],...
%      'ZData',[0 radius*cos(phi(1))],...
%      'Color','r')
% 
% cameratoolbar
% daspect([1 1 1])

end