%Initialize the environnement

dbstop if error
curPath = pwd;   
addpath(genpath(pwd));
addpath(genpath(fullfile(curPath,'..','..','PolymerChainDynamics','Code')));
addpath(genpath(fullfile(curPath,'..','..','Utils')));
addpath(genpath(fullfile(curPath,'..','..','Code')));
addpath(genpath(fullfile(curPath,'..','..','..','Utils')));

hName = system('HostName');
switch lower(hName)
    case('ofirDell-pc')
        
    case('zhaHome')
end
