%This program is used to simulate the spring-bead model
classdef RouseModel<handle
    properties
        pathsNormal@double %the paths of plymer with the normal coordinate;
        frictionCoefficient@double;%the frictionCoefficient of a bead;
        connectedBeads@double %an n by two array with pair wise bead numbers to connect
        fixedBeads@double %numBeads of beads that do not move
        b@double %length between 2 beads 
        rp = RandomWalkParams;
  
    end
    
    methods
        
        
        function obj=RouseModel(varargin)
            
           
            obj.frictionCoefficient = 1;
            obj.fixedBeads          = [];
            obj.connectedBeads      = [];
            obj.b                   =  1;
            obj.pathsNormal         =zeros(obj.rp.numParticles,3,obj.rp.numSteps);
            
             if ~isempty(varargin)
                    v = varargin;
                    if mod(size(varargin,2),2)~=0
                        error('name-value pair argument must be inserted')
                    end
                    
                    for vIdx = 1:(numel(v)/2)
                        obj.(v{2*vIdx-1})= v{2*vIdx};
                    end
             end
        end
        
    
        function Simulation(obj)
     
                % set initial position
                for j=1
                    noise = [zeros(1,obj.rp.dimension);...
                        sqrt(2*obj.rp.diffusionConst*obj.rp.dt)*randn(obj.rp.numParticles-1,obj.rp.dimension)];
                    obj.rp.paths(:,1:obj.rp.dimension,j)=cumsum(noise);
                    
                end
                
                a       = obj.rp.dimension*obj.rp.diffusionConst/obj.b^2;
                R       = RouseMatrix(obj.rp.numParticles, obj.connectedBeads, obj.fixedBeads);
                [U,D,V] = eig(R); 
               %transform to the noraml coordinate;
               obj.pathsNormal(:,:,j) = obj.rp.paths(:,:,j);
              
                for j=2:obj.rp.numSteps
                    noiseSingle                   = U*sqrt(2*obj.rp.diffusionConst*obj.rp.dt)*randn(obj.rp.numParticles,obj.rp.dimension);
                    obj.pathsNormal(:,:,j)        = U*obj.pathsNormal(:,:,j);
                    %zero out noise for fixed beads
                    noiseSingle(obj.fixedBeads,:) = 0;
                    obj.pathsNormal(:,:,j)        = -a*D*obj.pathsNormal(:,:,j-1)*obj.rp.dt+noiseSingle+obj.pathsNormal(:,:,j-1);
                 
                
                end
               
                for j=2:obj.rp.numSteps
                    
                    noiseSingle = sqrt(2*obj.rp.diffusionConst*obj.rp.dt)*randn(obj.rp.numParticles,obj.rp.dimension);
                  
                    
                    % zero out noise for fixed particles
                    noiseSingle(obj.fixedBeads,:) = 0;
                    
                    obj.rp.paths(:,1:3,j) = -a*R*obj.rp.paths(:,1:3,j-1)*obj.rp.dt+noiseSingle +obj.rp.paths(:,1:3,j-1);
                    
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
            
            x  = obj.rp.paths(:,1,1);
            y  = obj.rp.paths(:,2,1);
            z  = obj.rp.paths(:,3,1);
%             x1 = obj.pathsNormal(:,1,1); 
%             y1 = obj.pathsNormal(:,2,1);
%             z1 = obj.pathsNormal(:,3,1);
            
            l  =line('XData',x,'YData',y,'ZData',z,'Color','g','linestyle','-.','Marker','o','markersize',10,'Parent',d);
%            l1 =line('XData',x1,'YData',y1,'ZData',z1,'Color',[rand rand rand],'linestyle','-.','Marker','o','markersize',10,'Parent',d);
            for i=2:obj.rp.numSteps
                
                set(l,'XData',obj.rp.paths(:,1,i),'YData',obj.rp.paths(:,2,i),'ZData',obj.rp.paths(:,3,i));
                
                       pause(0.1)
                
                drawnow
%                hold on
%                 set(l1,'XData',obj.pathsNormal(:,1,i),'YData',obj.pathsNormal(:,2,i),'ZData',obj.pathsNormal(:,3,i));
%                 pause(0.1)
%                 drawnow
             end
            
            
        end
  
    end
end




