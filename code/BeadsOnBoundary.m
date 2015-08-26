function points = BeadsOnBoundary(initialPoint,radius,domainCenter, dt, diffusionConst,beadIndices)
% Choose the position of the beads on the boundary, starting from a random
% position and diffusing on the boundary to chose the remain

% make sure the beadIndices are sorted 
beadIndices = sort(beadIndices);
numSteps    = beadIndices(end)-beadIndices(1)+1;
% preallocation
% points  = zeros(length(beadIndices),3);

% choose a random position for the first bead
% phi           = 0+(pi-0)*rand;
% theta         = 0+(2*pi-0)*rand ; 
% points(1,:)   = radius.*[sin(phi)*cos(theta), sin(phi)*sin(theta),cos(phi)];
%pathOnSurface = DiffuseOnSphere(initialPoint,numSteps,radius,domainCenter,dt,diffusionConst);
w = 0;
mag =0.1;
pathOnSurface = DiffusionOnSphere(initialPoint,dt,diffusionConst,numSteps,domainCenter,radius,w,mag);
pointIndices  = beadIndices -beadIndices(1)+1; 
points        = pathOnSurface(pointIndices,:);

