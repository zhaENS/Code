%Initialize the environnement

dbstop if error
curPath = pwd;   
addpath(genpath(pwd));
addpath(genpath(fullfile(curPath,'..','..','polymerChainDynamics','Code')));
addpath(genpath(fullfile(curPath,'..','..','Utils')));
addpath(genpath(fullfile(curPath,'..','..','Code')));
addpath(genpath(fullfile(curPath,'..','..','..','Utils')));
addpath(genpath(fullfile(curPath,'..','PolymerChainDynamicsResults')));
hName = system('HostName');
switch lower(hName)
    case('ofirDell-pc')
        
    case('PC-20150709IFBW')
end
