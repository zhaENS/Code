%this a test for the function DiffusionOnSphere;


close all
initialPoint = randn(1,3);
dt           = 0.1;
diffConst    = 0.1;
numSteps     = 100;
domainCenter = [0 0 0];
radius       = sqrt(sum(initialPoint.^2));
w =0;
paths = DiffusionOnSphere(initialPoint,dt,diffConst,numSteps,domainCenter,radius,w);