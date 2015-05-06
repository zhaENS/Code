%This program is used to simulate the first time such that 3 beads have met
%in each simulation
classdef RouseModelMetTime<handle
    properties
        dimension%dimension=1,2 or 3
        numParticles%number of particles in polymers;
        dt%pas de temps
        numSteps% number of motion
        numSimulations%number of simulations;
        diffusionConst %constante diffusion
        paths %the paths of polymer;
        frictionCoefficient;%the frictionCoefficient of a bead;
        metBeedNum; %the number of the beed which have met with the first and last beed;
        connectedBeads %an n by two array with pair wise bead numbers to connect
        fixedBeads %numBeads of beads that do not move
        encounterDistance%the distance between 3 beeds such that they have met
        encounterTime%a tableau contains the time they have met for every simulation
        b %length between 2 beeds;

    end
    
    methods
        
        
        function obj=RouseModelMetTime(dimension,numParticles,dt,diffusionConst,...
                numSteps,frictionCoefficient,connectedBeads, fixedBeads ,...
                numSimulations,metBeedNum,b,encounterDistance)
            
            obj.dimension = dimension;
            obj.numParticles = numParticles;
            obj.dt = dt;
            obj.numSteps = numSteps;
            obj.diffusionConst = diffusionConst;
            obj.paths = zeros(obj.numParticles,3,2);
            obj.frictionCoefficient = frictionCoefficient;
            obj.fixedBeads= fixedBeads;
            obj.connectedBeads = connectedBeads ;
            obj.numSimulations = numSimulations;
            obj.metBeedNum = metBeedNum;
            obj.b = b;
            obj.encounterDistance = encounterDistance;
            obj.encounterTime = zeros(obj.numSimulations,length(obj.b));
        
        end
        
        
        function Simulation(obj)
            % set initial position
       for i=1:length(obj.b)
  
            for s = 1:obj.numSimulations
                t1=clock;
                exitFlag = false;
                step     = 1;
                
                % set initial position
                for j=1
                    noise = [zeros(1,obj.dimension);...
                        sqrt(2*obj.diffusionConst*obj.dt)*randn(obj.numParticles-1,obj.dimension)];
                    obj.paths(:,1:obj.dimension,step)=cumsum(noise);
                    
                end
                
                a = obj.dimension*obj.diffusionConst*obj.frictionCoefficient/obj.b(i)^2;%the spring const
                R = RouseMatrix(obj.numParticles, obj.connectedBeads, obj.fixedBeads);
                step     = 2;
                while ~exitFlag && step <= obj.numSteps
                    
                 
                    
                    noiseSingle = sqrt(2*obj.diffusionConst*obj.dt)*randn(obj.numParticles,obj.dimension);
                   
                    % zero out noise for fixed particles
                    noiseSingle(obj.fixedBeads,:) = 0;
                    
                    obj.paths(:,:,2) = -a/obj.frictionCoefficient*R*obj.paths(:,:,1)*obj.dt+noiseSingle +obj.paths(:,:,1);
                    obj.paths(:,:,1) = obj.paths(:,:,2); 
                    centerMasse = 1/3*sum(obj.paths(obj.metBeedNum,:,2));
                    A           = zeros(numel(obj.metBeedNum)+1,3);
                    A(1:3,1:3)  = obj.paths(obj.metBeedNum,:,2);
                    A(end,1:3)  = centerMasse;
                  %  D = pdist2mex(obj.paths(obj.metBeedNum,:,2)',obj.paths(obj.metBeedNum,:,2)','euc',[],[],[]);%calculate the distance
                   D = pdist2mex(A',A','euc',[],[],[],[]);
 %                    D = Distance(obj.paths(1,:,1),obj.paths(obj.numParticles,:,1),obj.paths(obj.metBeedNum,:,1));
                    exitFlag = all(D(end,:)<obj.encounterDistance(i));
 %                    exitFlag = all([D(1),D(2),D(3)] < obj.encounterDistance(i));
                    
                    if exitFlag
                        obj.encounterTime(s,i)=step*obj.dt ;
                    
                    end
               
                   
                    step = step+1;
                    
                end
              t2 = clock;
                sprintf('%s%s%s','simulation ',num2str(s),' is done')
                sprintf('the time of simulation %s is %s',num2str(s),num2str(etime(t2,t1)))
            end
        end
    end
        
        
        function Plot(obj)
            
            for i=1:length(obj.b)
                figure(i+1)
                [h,bins]=hist(obj.encounterTime(:,i),50);
                h=h./trapz(bins,h);
                bar(bins,h)
            end
        end
       
    end
end
 




