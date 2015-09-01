%this a test for the function DiffusionOnSphere;


close all
initialPoint = randn(1,3);
dt           = 0.1;
diffConst    = 0.1;
numSteps     = 100;
domainCenter = [0 0 0];
radius       = sqrt(sum(initialPoint.^2));
w =1/2;
ang = linspace(0.1,2*pi,numSteps-1);
ang = ang';
paths = DiffusionOnSphere(initialPoint,dt,diffConst,numSteps,domainCenter,radius,ang,w);