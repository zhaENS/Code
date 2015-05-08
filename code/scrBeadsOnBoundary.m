
addpath(genpath(fullfile(pwd,'..','..','PolymerChainDynamics')));
dp = DomainHandlerParams;
dp.domainShape    = 'sphere';
dp.domainWidth    = sqrt(14);
dp.domainCenter   = [0 0 0];
dp.showDomain     = true;
dp.dt             = 0.01;
dp.diffusionConst = 0.1;
initialPoint      = [2 1 3];


%points = DiffuseOnSphere(initialPoint,numSteps,radius, dt, diffusionConst);
domainClass = DomainHandler(dp);
beads       = [ 2 5 7 10];% beads on the boundary
%points     = DiffuseOnSphere(initialPoint,numSteps,radius, dt, diffusionConst)
points      = BeadsOnBoundary(dp.domainWidth, dp.dt, dp.diffusionConst,beads);
numSteps    = 20;
k           = 1;
for i = 1:size(points,1)-1
   
    path = BrownianBridgeSim(points(i,:),points(i+1,:),numSteps,domainClass);
       

pathTotal(:,:,k) = path;
k=k+1;
end

% plot 
[sx, sy, sz]= sphere(20);
sx = sx*dp.domainWidth;
sy = sy*dp.domainWidth;
sz = sz*dp.domainWidth;
c=rand(1,3);
figure, mesh(sx,sy,sz,'FaceColor','none'), hold on 
plot3(points(2:3,1), points(2:3,2),points(2:3,3),'o-')
hold on
plot3(pathTotal(:,1,2),pathTotal(:,2,2),pathTotal(:,3,2),'Color',c);
daspect([1 1 1]), cameratoolbar
