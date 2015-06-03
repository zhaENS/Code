%This program is used to test the FirstTry class;


rp = RandomWalkParams('dimension',3,'numSteps',100,'diffusionConst',0.1,'dt',0.1,'numParticles',200);
                                

h = FirstTry(rp);

h.Calculate;
h.Plot;