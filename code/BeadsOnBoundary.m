function points = BeadsOnBoundary(radius, dt, diffusionConst,beads)

points =zeros(length(beads),3);

% phi         = 0+(pi-0)*rand;
% theta       = 0+(2*pi-0)*rand ; 
%points(1,:) = [radius*sin(phi)*cos(theta) radius*sin(phi)*sin(theta)%radius*cos(phi)];*
points(1,:) = [2 1 3];
%radius      = sqrt(sum(points(1,:).^2));

for i = 1:length(beads)-1

pointsTotal   = DiffuseOnSphere(points(i,:),beads(i+1)-beads(i),radius, dt, diffusionConst);
points(i+1,:) = pointsTotal(end,:);
end
