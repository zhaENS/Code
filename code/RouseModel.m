%This program is used to simulate the spring-bead model
classdef RouseModel<handle
    properties
        dimension%dimension=1,2 or 3
        numParticles%number of particles in polymers;
        dt%pas de temps
        numSteps% number of motion
        numSimulations%number of simulations;
        diffusionConst %constante diffusion
        paths %the paths of polymer;
        pathsNormal %the paths of plymer with the normal coordinate;
        frictionCoefficient;%the frictionCoefficient of a bead;
        MetBeedNum; %the number of the beed which have met with the first and last beed;
        connectedBeads %an n by two array with pair wise bead numbers to connect
        fixedBeads %numBeads of beads that do not move
        b %length between 2 beads 
    end
    
    methods
        
        
        function obj=RouseModel(dimension,numParticles,dt,diffusionConst,...
                numSteps,frictionCoefficient,connectedBeads, fixedBeads,b)
            
            obj.dimension           = dimension;
            obj.numParticles        = numParticles;
            obj.dt                  = dt;
            obj.numSteps            = numSteps;
            obj.diffusionConst      = diffusionConst;
            obj.paths               = zeros(obj.numParticles,3,obj.numSteps);
            obj.pathsNormal         = zeros(obj.numParticles,3,obj.numSteps);
            obj.frictionCoefficient = frictionCoefficient;
            obj.fixedBeads          = fixedBeads;
            obj.connectedBeads      = connectedBeads ;
            obj.b                   =b;
        end
        
        
        function Simulation(obj)
     
                % set initial position
                for j=1
                    noise = [zeros(1,obj.dimension);...
                        sqrt(2*obj.diffusionConst*obj.dt)*randn(obj.numParticles-1,obj.dimension)];
                    obj.paths(:,1:obj.dimension,j)=cumsum(noise);
                    
                end
                
                a       = obj.dimension*obj.diffusionConst/obj.b^2;
                R       = RouseMatrix(obj.numParticles, obj.connectedBeads, obj.fixedBeads);
                [U,D,V] = eig(R); 
               %transform to the noraml coordinate;
               obj.pathsNormal(:,:,j) = obj.paths(:,:,j);
              
                for j=2:obj.numSteps
                    noiseSingle                   = U*sqrt(2*obj.diffusionConst*obj.dt)*randn(obj.numParticles,obj.dimension);
                    obj.pathsNormal(:,:,j)        = U*obj.pathsNormal(:,:,j);
                    %zero out noise for fixed beads
                    noiseSingle(obj.fixedBeads,:) = 0;
                    obj.pathsNormal(:,:,j)        = -a*D*obj.pathsNormal(:,:,j-1)*obj.dt+noiseSingle+obj.pathsNormal(:,:,j-1);
                 
                
                end
               
                for j=2:obj.numSteps
                    
                    noiseSingle = sqrt(2*obj.diffusionConst*obj.dt)*randn(obj.numParticles,obj.dimension);
                  
                    
                    % zero out noise for fixed particles
                    noiseSingle(obj.fixedBeads,:) = 0;
                    
                    obj.paths(:,1:3,j) = -a*R*obj.paths(:,1:3,j-1)*obj.dt+noiseSingle +obj.paths(:,1:3,j-1);
                    
                end
%               
            end
     
        
        function Plot(obj)
            f=figure;
            b=[-10 10];
            d= axes('Parent',f,'NextPlot','replaceChildren','XLim',b,'YLim',b,'ZLim',b);
            daspect([1 1 1]);
            cameratoolbar
           % c=rand(1,3);
            
            x  = obj.paths(:,1,1);
            y  = obj.paths(:,2,1);
            z  = obj.paths(:,3,1);
            x1 = obj.pathsNormal(:,1,1); 
            y1 = obj.pathsNormal(:,2,1);
            z1 = obj.pathsNormal(:,3,1);
            
            l  =line('XData',x,'YData',y,'ZData',z,'Color','g','linestyle','-.','Marker','o','markersize',10,'Parent',d);
            l1 =line('XData',x1,'YData',y1,'ZData',z1,'Color',[rand rand rand],'linestyle','-.','Marker','o','markersize',10,'Parent',d);
            for i=2:obj.numSteps
                
                set(l,'XData',obj.paths(:,1,i),'YData',obj.paths(:,2,i),'ZData',obj.paths(:,3,i));
                
                       pause(0.1)
                
                drawnow
               hold on
                set(l1,'XData',obj.pathsNormal(:,1,i),'YData',obj.pathsNormal(:,2,i),'ZData',obj.pathsNormal(:,3,i));
                pause(0.1)
                drawnow
            end
            
            
        end
  
    end
end




