classdef RandomWalkParams<handle
    
properties

dimension      = 3;;%dimension=1,2 or 3
numParticles   = 100;%number of particles in polymers;
dt             =0.1;%pas de temps
numSteps       =100;% number of motion
diffusionConst =0.1; %constante diffusion
paths          =[];  ; %the paths of polymer;



end
methods

 function obj = RandomWalkParams(varargin)
                % parse name-value pair input parameters
                if ~isempty(varargin)
                    v = varargin;
                    if mod(numel(varargin{:}),2)~=0
                        error('name-value pair argument must be inserted')
                    end
                    
                    for vIdx = 1:(numel(v)/2)
                        obj.(v{2*vIdx-1})= v{2*vIdx};
                    end
                end
        end

 
end

end


