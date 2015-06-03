%This program is used to simulate the brownian motion
classdef FirstTry<handle

properties


FT = RandomWalkParams; 


end
  
methods
 
    function obj = FirstTry(rp)
      
          obj.FT =  rp;
      
      end
   
    
    function Calculate(obj)
        %Initial position 
        for j=1
            noise = [zeros(1,obj.FT.dimension);...
                     sqrt(2*obj.FT.diffusionConst*obj.FT.dt)*randn(obj.FT.numParticles-1,obj.FT.dimension)];
            obj.FT.paths(:,1:obj.FT.dimension,j)=cumsum(noise);
             
        end
        
        for j=2:obj.FT.numSteps
            noise = [sqrt(2*obj.FT.diffusionConst*obj.FT.dt)*randn(obj.FT.numParticles,obj.FT.dimension)];
            obj.FT.paths(:,1:obj.FT.dimension,j)=obj.FT.paths(:,1:obj.FT.dimension,j-1)+noise;
             
        end
        
           
    end %'random walk 'simulation
    
    function Plot(obj)
       f=figure;
       b =[-20 20];
       a= axes('Parent',f,'NextPlot','replaceChildren','XLim',b,'YLim',b,'ZLim',b); 
       daspect([1 1 1]);
       cameratoolbar
       c=rand(1,3);
       
           x = obj.FT.paths(:,1,1);
           y = obj.FT.paths(:,2,1);
           z = obj.FT.paths(:,3,1);
%        end
       l=line('XData',x,'YData',y,'ZData',z,'Color',c,'linestyle','-.','Marker','o','markersize',10,'Parent',a);
          for i=2:obj.FT.numSteps
         
               set(l,'XData',obj.FT.paths(:,1,i),'YData',obj.FT.paths(:,2,i),'ZData',obj.FT.paths(:,3,i));
               
           %   pause(3)
              
               drawnow
              
          end
    end
    
end
end
