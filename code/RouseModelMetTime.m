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
      %  metBeadNum; %the number of the beed which have met with the first and last beed;
        metBeadNum1;
        connectedBeads %an n by two array with pair wise bead numbers to connect
        fixedBeads %numBeads of beads that do not move
        encounterDistance%the distance between 3 beeds such that they have met
        encounterTime%a tableau contains the time they have met for every simulation
        encounterTime1;
        b %length between 2 beeds;

    end
    
    methods
        
        
        function obj=RouseModelMetTime(dimension,numParticles,dt,diffusionConst,...
                numSteps,frictionCoefficient,connectedBeads, fixedBeads ,...
                numSimulations,metBeadNum1,b,encounterDistance)
            
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
         %  obj.metBeadNum  = metBeadNum;
            obj.metBeadNum1 = metBeadNum1;
            obj.b = b;
            obj.encounterDistance = encounterDistance;
          %  obj.encounterTime = zeros(obj.numSimulations,size(obj.metBeadNum,1));
            obj.encounterTime1 = zeros(obj.numSimulations,size(obj.metBeadNum1,1));
        
        end
        
        
        function Simulation(obj)
            % set initial position
      
            for i = 1:size(obj.metBeadNum1,1)
            for s = 1:obj.numSimulations
                t1=clock;
              %  exitFlag  = false;
                exitFlag1 = false;
                step     = 1;
                
                % set initial position
                for j=1
                    noise = [zeros(1,obj.dimension);...
                        sqrt(2*obj.diffusionConst*obj.dt)*randn(obj.numParticles-1,obj.dimension)];
                    obj.paths(:,1:obj.dimension,step)=cumsum(noise);
                    
                end
                
                a = obj.dimension*obj.diffusionConst*obj.frictionCoefficient/obj.b^2;%the spring const
                R = RouseMatrix(obj.numParticles, obj.connectedBeads, obj.fixedBeads);
                step     = 2;
%                 while ~exitFlag && step <= obj.numSteps
%                     
%                  
%                     
%                     noiseSingle = sqrt(2*obj.diffusionConst*obj.dt)*randn(obj.numParticles,obj.dimension);
%                    
%                     % zero out noise for fixed particles
%                     noiseSingle(obj.fixedBeads,:) = 0;
%                     
%                     obj.paths(:,:,2) = -a/obj.frictionCoefficient*R*obj.paths(:,:,1)*obj.dt+noiseSingle +obj.paths(:,:,1);
%                     obj.paths(:,:,1) = obj.paths(:,:,2); 
%                  %   centerMasse = 1/3*sum(obj.paths(obj.metBeadNum1(i,:),:,2));
%                     A           = zeros(numel(obj.metBeadNum(i,:)),3);
%                %     A1          = zeros(numel(obj.metBeadNum1(i,:)+1),3);
%                     A(1:2,1:3)  = obj.paths(obj.metBeadNum(i,:),:,2);
%                 %    A1(1:3,:)   = obj.paths(obj.metBeadNum1(i,:),:,2);
%                  %   A1(end,1:3)  = centerMasse;
%                   %  D = pdist2mex(obj.paths(obj.metBeadNum,:,2)',obj.paths(obj.metBeadNum,:,2)','euc',[],[],[]);%calculate the distance
%                     D = pdist2mex(A',A','euc',[],[],[]);
%                %     D1 = pdist2mex(A1',A1','euc',[],[],[],[]);
%  %                    D = Distance(obj.paths(1,:,1),obj.paths(obj.numParticles,:,1),obj.paths(obj.metBeadNum,:,1));
%                     exitFlag = all(D(end,:)<obj.encounterDistance);
%  %                    exitFlag = all([D(1),D(2),D(3)] < obj.encounterDistance(i));
%                     
%                     if exitFlag
%                         obj.encounterTime(s,i)=step*obj.dt ;
%                     
%                     end
%                
%                    
%                     step = step+1;
%                     
%                 end
                while ~exitFlag1 && step <= obj.numSteps
                    
                 
                    
                    noiseSingle = sqrt(2*obj.diffusionConst*obj.dt)*randn(obj.numParticles,obj.dimension);
                   
                    % zero out noise for fixed particles
                    noiseSingle(obj.fixedBeads,:) = 0;
                    
                    obj.paths(:,:,2) = -a/obj.frictionCoefficient*R*obj.paths(:,:,1)*obj.dt+noiseSingle +obj.paths(:,:,1);
                    obj.paths(:,:,1) = obj.paths(:,:,2); 
                    centerMasse = 1/3*sum(obj.paths(obj.metBeadNum1(i,:),:,2));
                 %   A           = zeros(numel(obj.metBeadNum(i,:)),3);
                    A1          = zeros(numel(obj.metBeadNum1(i,:)+1),3);
                  %  A(1:2,1:3)  = obj.paths(obj.metBeadNum(i,:),:,2);
                    A1(1:3,:)   = obj.paths(obj.metBeadNum1(i,:),:,2);
                    A1(end,1:3) = centerMasse;
                  %  D = pdist2mex(obj.paths(obj.metBeadNum,:,2)',obj.paths(obj.metBeadNum,:,2)','euc',[],[],[]);%calculate the distance
                   % D = pdist2mex(A',A','euc',[],[],[]);
                    D1 = pdist2mex(A1',A1','euc',[],[],[],[]);
 %                    D = Distance(obj.paths(1,:,1),obj.paths(obj.numParticles,:,1),obj.paths(obj.metBeadNum,:,1));
                    exitFlag1 = all(D1(end,:)<obj.encounterDistance);
 %                    exitFlag = all([D(1),D(2),D(3)] < obj.encounterDistance(i));
                    
                    if exitFlag1
                        obj.encounterTime1(s,i)=step*obj.dt ;
                    
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
           for i = 1:size(obj.metBeadNum1,1)
                figure(i+1)
%             %   subplot(2,1,1)
%                  [h,bins]=hist(obj.encounterTime(:,i),50);
%                  h=h./trapz(bins,h);
%                   %bar(bins,h)
%                   plot(bins,h,'r')
%                   legend('the probability of encounter [1 16 32]')
%                title(['beads ',num2str(obj.metBeadNum(i,:))])
%                       % num2str(sum(obj.encounterTime(:,i))/obj.numSimulations)]);
%           hold on
                  [h,bins]=hist(obj.encounterTime1(:,i),50);
                 h=h./trapz(bins,h);
                  bar(bins,h)
                   %plot(bins,h)
                %  legend('the probability of encounter [1 32]');
               title(['beads ',num2str(obj.metBeadNum1(i,:))])
                    %   num2str(sum(obj.encounterTime1(:,i))/obj.numSimulations)]);
         
           end
        %figure
%         meanTime=sum(rouseModel.encounterTime)./rouseModel.numSimulations;
%         x=[1:1:14];plot(meanTime);set(gca,'XTick',[1:1:14]);legend('mean encounterTime');xlabel('metBeads from[1 2 32] to[1 15 32]');ylabel('mean encounterTime'); 
%         text(x,meanTime,num2cell(meanTime));
%         end
       
    end
    end
end
 




