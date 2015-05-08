%Create a brownian bridge from beadStart to beadEnd
function path = BrownianBridgeSim(beadStart,beadEnd,numSteps,domainClass)

path(1,:) = beadStart;

for i = 2:numSteps-1
    d = 0;
    while d==0
    path(i,:) = (1-1/numSteps*i)*beadStart+1/numSteps*i*beadEnd+rand(1,3);
    d         =  domainClass.InDomain(path(i,:));
    end
    disp('***')
end

path(numSteps,:) = beadEnd;

end