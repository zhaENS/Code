function points = BeadsOnBoundary(radius, dt, diffusionConst,beadIndices)
% Choose the position of the beads on the boundary, strting from a random
% position and diffusing on the boundary to chose the rest

% make sure the beadIndices are sorted 
beadIndices = sort(beadIndices);
numSteps    = beadIndices(end)-beadIndices(1)+1;
% preallocation
points  = zeros(length(beadIndices),3);

% choose a random position for the first bead
phi           = 0+(pi-0)*rand;
theta         = 0+(2*pi-0)*rand ; 
points(1,:)   = radius.*[sin(phi)*cos(theta), sin(phi)*sin(theta),cos(phi)];
pathOnSurface = DiffuseOnSphere(points(1,:),numSteps,radius, dt, diffusionConst);
pointIndices  = beadIndices -beadIndices(1)+1; 
points        = pathOnSurface(pointIndices,:);

% %====
% for bIdx = 1:length(beadIndices)-1
%     numSteps         = beadIndices(bIdx+1)-beadIndices(bIdx);
%     pointsTotal      = DiffuseOnSphere(points(bIdx,:),numSteps,radius, dt, diffusionConst);
%     points(bIdx+1,:) = pointsTotal(end,:);% take the last point 
% end
