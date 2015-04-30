function R = RouseMatrix(numBeads, connectedBeads, fixedBeads)
% produce the numBeads X numBeas Rouse matrix.
% numBeads- a positive integer
% connectedBeads- an n by two array with pair wise bead numbers to connect
% fixedBeads - positive integer <= numBeads of beads that do not move

if ~exist('connectedBeads','var')
    connectedBeads = []; % non linear connections
end
if ~exist('fixedBeads','var')
    fixedBeads = []; % connectors that do not move
end
    
            R = zeros(numBeads);
            R = R+ diag(-1*ones(1,numBeads-1),-1) +...
                diag(-1*ones(1,numBeads-1),1);
            for bIdx = 1:size(connectedBeads,1)
                R(connectedBeads(bIdx,1),connectedBeads(bIdx,2)) = -1;
                R(connectedBeads(bIdx,2),connectedBeads(bIdx,1)) = -1;
            end
            R(fixedBeads,:) = 0;
            R = R+diag(sum(R==-1,2));
end