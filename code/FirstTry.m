%This program is used to simulate the brownian motion
classdef FirstTry<handle

properties
dimension%dimension=1,2 or 3
numParticles%number of particles in polymers;
dt%pas de temps
numSteps% number of motion
diffusionConst %constante diffusion
paths %the paths of polymer;
end
  
methods
    
    %class constructor
    function obj=FirstTry(dimension,numParticles,dt,diffusionConst,numSteps)
       obj.dimension=dimension;
       obj.numParticles=numParticles;
       obj.dt=dt;
       obj.numSteps=numSteps;
       obj.diffusionConst=diffusionConst;
       obj.paths = zeros(obj.numParticles,3,obj.numSteps);
    end
     
    
    function Calculate(obj)
        
        for j=1
            noise = [zeros(1,obj.dimension);...
                     sqrt(2*obj.diffusionConst*obj.dt)*randn(obj.numParticles-1,obj.dimension)];
            obj.paths(:,1:obj.dimension,j)=cumsum(noise);
             
        end
        
        for j=2:obj.numSteps
            noise = [sqrt(2*obj.diffusionConst*obj.dt)*randn(obj.numParticles,obj.dimension)];
            obj.paths(:,1:obj.dimension,j)=obj.paths(:,1:obj.dimension,j-1)+noise;
             
        end
        
           
    end %'random walk 'simulation
    
    function Plot(obj)
       f=figure;
       b =[-20 20];
       a= axes('Parent',f,'NextPlot','replaceChildren','XLim',b,'YLim',b,'ZLim',b); 
       daspect([1 1 1]);
       cameratoolbar
       c=rand(1,3);
       
           x = obj.paths(:,1,1);
           y = obj.paths(:,2,1);
           z = obj.paths(:,3,1);
%        end
       l=line('XData',x,'YData',y,'ZData',z,'Color',c,'linestyle','-.','Marker','o','markersize',10,'Parent',a);
          for i=2:obj.numSteps
         
               set(l,'XData',obj.paths(:,1,i),'YData',obj.paths(:,2,i),'ZData',obj.paths(:,3,i));
               
              pause(3)
              
               drawnow
              
          end
    end
    
end
end
