%this program is used to simulate the probability distribution of Gaussian 
classdef Idealchain<handle

properties
dimension%dimension=3
numParticles%number of particles in polymers;
dt%pas de temps
numSteps% number of motion
diffusionConst %constante diffusion
paths %the paths of polymer;
R %end to end distance 
simulation %number of simulations
end
  
methods
    
    %class constructor
    function obj=Idealchain(dimension,numParticles,dt,diffusionConst,numSteps,simulation)
       obj.dimension=dimension;
       obj.numParticles=numParticles;
       obj.dt=dt;
       obj.numSteps=numSteps;
       obj.diffusionConst=diffusionConst;
       obj.simulation=simulation;
       obj.paths = zeros(obj.numParticles,3,obj.numSteps);
       obj.R=zeros(obj.simulation,3);
    end
     
    
    function Calculate(obj)
       for s=1:obj.simulation 
        for j=1
            noise = [zeros(1,obj.dimension);...
                     sqrt(2*obj.diffusionConst*obj.dt)*randn(obj.numParticles-1,obj.dimension)];
            obj.paths(:,1:obj.dimension,j)=cumsum(noise);
             
        end
        
        for j=2:obj.numSteps
            noise = [sqrt(2*obj.diffusionConst*obj.dt)*randn(obj.numParticles,obj.dimension)];
            obj.paths(:,1:obj.dimension,j)=obj.paths(:,1:obj.dimension,j-1)+noise;
            
            obj.R(s,:)=cumsum(obj.paths(obj.numParticles,:,j)-obj.paths(1,:,j-1));   
        end
        
      
        
       end
           
    end
    
    function Plot(obj)
%        f=figure;
%        b =[-4 4];
%        a= axes('Parent',f,'NextPlot','replaceChildren','XLim',b,'YLim',b,'ZLim',b); 
%        c=rand(1,3);
%        
%            x = obj.paths(:,1,1);
%            y = obj.paths(:,2,1);
%            z = obj.paths(:,3,1);
% %        end
%        l=line('XData',x,'YData',y,'ZData',z,'Color',c,'linestyle','-.','Marker','o','markersize',10,'Parent',a);
%           for i=2:obj.numSteps
%            
%                set(l,'XData',obj.paths(:,1,i),'YData',obj.paths(:,2,i),'ZData',obj.paths(:,3,i));
%                
% %                pause(0.1)
%               
%                drawnow
%               
%           end
          figure(2);
          
           s=[1:obj.simulation];
         
            plot(s,obj.R(:,1));
          
         
    end
    
end
            
    
 
    
    

end
    