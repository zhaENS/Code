%Create a brownian bridge from beadStart to beadEnd
function path = BrownianBridge(beadStart,beadEnd,numSteps)

for i = 1:numSteps
    
    path(i,:) = (1-1/numSteps*i)*beadStart+1/numSteps*i*beadEnd+randn(1,3);


end


end