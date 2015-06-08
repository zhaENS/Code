%This program is used to test the class IdealChain in order to verify the
%PDF is gaussian;


rp = RandomWalkParams('dimension',3,'numSteps',100,'simulation',1000);

iC = Idealchain(rp);

iC.Calculate
iC.Plot