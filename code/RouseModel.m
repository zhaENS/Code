classdef RouseModel<handle
    properties
        dimension%dimension=1,2 or 3
        numParticles%number of particles in polymers;
        dt%pas de temps
        numSteps% number of motion
        numSimulations%number of simulations;
        diffusionConst %constante diffusion
        paths %the paths of polymer;
        frictionCoefficient;%the frictionCoefficient of a bead;
        MetBeedNum; %the number of the beed which have met with the first and last beed;
        connectedBeads %an n by two array with pair wise bead numbers to connect
        fixedBeads %numBeads of beads that do not move
     
    end
    
    methods
        
        
        function obj=RouseModel(dimension,numParticles,dt,diffusionConst,...
                numSteps,frictionCoefficient,connectedBeads, fixedBeads)
            
            obj.dimension = dimension;
            obj.numParticles = numParticles;
            obj.dt = dt;
            obj.numSteps = numSteps;
            obj.diffusionConst = diffusionConst;
            obj.paths = zeros(obj.numParticles,3,obj.numSteps);
            obj.frictionCoefficient = frictionCoefficient;
            obj.fixedBeads= fixedBeads;
            obj.connectedBeads = connectedBeads ;
           
        end
        
        
        function Simulation(obj)
     
                % set initial position
                for j=1
                    noise = [zeros(1,obj.dimension);...
                        sqrt(2*obj.diffusionConst*obj.dt)*randn(obj.numParticles-1,obj.dimension)];
                    obj.paths(:,1:obj.dimension,j)=cumsum(noise);
                    
                end
                
                a = obj.diffusionConst/obj.frictionCoefficient;
                R = RouseMatrix(obj.numParticles, obj.connectedBeads, obj.fixedBeads);

                for j=2:obj.numSteps
                    
                    noiseSingle = sqrt(2*obj.diffusionConst*obj.dt)*randn(obj.numParticles,obj.dimension);
                    %                 obj.paths(:,:,j)=obj.paths(1,:,j-1)+obj.dt*(a*(obj.paths(2,:,j-1)-obj.paths(1,:,j-1)))+noiseSingle;
                    
                    % zero out noise for fixed particles
                    noiseSingle(obj.fixedBeads,:) = 0;
                    
                    obj.paths(:,1:3,j) = -a*R*obj.paths(:,1:3,j-1)*obj.dt+noiseSingle +obj.paths(:,1:3,j-1);
               
                
                    %               for i=2:obj.numParticles-1
                    %                 noiseSingle = sqrt(2*obj.diffusionConst*obj.dt)*randn(1,obj.dimension);
                    %                 obj.paths(i,:,j)=obj.paths(i,:,j-1)+obj.dt*(a*(obj.paths(i+1,:,j-1)...
                    %                     +obj.paths(i-1,:,j-1)-2*obj.paths(i,:,j-1)))+noiseSingle;
                    %
                    %               end
                    %               noiseSingle = sqrt(2*obj.diffusionConst*obj.dt)*randn(1,obj.dimension);
                    %               obj.paths(obj.numParticles,:,j)=obj.paths(obj.numParticles,:,j-1)+obj.dt*...
                    %                   (a*(obj.paths(obj.numParticles-1,:,j-1)-obj.paths(obj.numParticles,:,j-1)))+noiseSingle;
                    
                  
                    %             end
                    
                end
               
            end
     
        
        function Plot(obj)
            f=figure;
            b=[-10 10];
            d= axes('Parent',f,'NextPlot','replaceChildren','XLim',b,'YLim',b,'ZLim',b);
            daspect([1 1 1]);
            cameratoolbar
            c=rand(1,3);
            
            x = obj.paths(:,1,1);
            y = obj.paths(:,2,1);
            z = obj.paths(:,3,1);
            
            l=line('XData',x,'YData',y,'ZData',z,'Color',c,'linestyle','-.','Marker','o','markersize',10,'Parent',d);
            for i=2:obj.numSteps
                
                set(l,'XData',obj.paths(:,1,i),'YData',obj.paths(:,2,i),'ZData',obj.paths(:,3,i));
                
                        pause(0.5)
                
                drawnow
            end
            
            
        end
  
    end
end




