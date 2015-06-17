%This function is used to do the diffusion on sphere
%initialPoint : the initial point which is on the sphere
%dt :time step
%numSteps :number of steps that walks on the sphere
%domainCenter : the center of mass of sphere
%diffConst : a diffusion const for the processus 
%radius : the radius of the sphere
function points = DiffusionOnSphere(initialPoint,dt,diffConst,numSteps,domainCenter,radius)
%check initialPoint is on the sphere or not
rho = sum(initialPoint.^2,2);
if (rho-radius)^2>eps
    error('the initial point is not on the sphere');
end

%Intersection between sphere and line
%t = radius/


end