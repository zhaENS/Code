%this program is used to simulate the probability distribution of Gaussian 
classdef Idealchain<handle

properties


rp = RandomWalkParams;%the parameters of RandomWalk;
endToEndDist %end to end distance 

end
  
methods
    
    %class constructor
    function obj=Idealchain(randomWalkParams)

    obj.rp = randomWalkParams;
    end
     
    
    function Calculate(obj)
       for s=1:obj.rp.simulation 
        %step 1:connect the chain;
        obj.rp.paths(1,1:obj.rp.dimension,1)=[0 0 0];%position of first beed;
        obj.rp.paths(2,1:obj.rp.dimension,1)=obj.rp.paths(1,1:obj.rp.dimension,1);
        noise=sqrt(2*obj.rp.diffusionConst*obj.rp.dt)*randn(1,obj.rp.dimension);
        for i=1:obj.rp.numParticles
        obj.rp.paths(2,1:obj.rp.dimension,1)=obj.rp.paths(2,1:obj.rp.dimension,1)+noise;%position of last beed;
        end
        
        %step 2:end :the position of beedStart and beedEnd varity by time;
        for j=2:obj.rp.numSteps
                obj.rp.paths(1,1:obj.rp.dimension,2)=obj.rp.paths(1,1:obj.rp.dimension,1)+noise;
                obj.rp.paths(2,1:obj.rp.dimension,2)=obj.rp.paths(2,1:obj.rp.dimension,1)+noise;   
                obj.rp.paths(1,1:obj.rp.dimension,1)=obj.rp.paths(1,1:obj.rp.dimension,2);
                obj.rp.paths(2,1:obj.rp.dimension,1)=obj.rp.paths(2,1:obj.rp.dimension,2);  
        end
        obj.endToEndDist(s,:)=obj.rp.paths(2,:,2)-obj.rp.paths(1,:,2);
        
        
       end
           
    end
    
    function Plot(obj)
      
%plot the endToEndDist in each coordinate;
        f=@(N,R,b)(sqrt(1/(2*pi*N*b^2))*exp(-(sum(R.^2,2))./(2*N*b^2)));%PDF 
        subplot(2,2,1)
        [h, bins]=hist(obj.endToEndDist(:,1),50); h= h./sum(h); bar(bins,h),...
        hold on, 
       plot(bins,f(obj.rp.numParticles,[bins'],1),'r')
        title('PDF on x')
        
        subplot(2,2,2)
        [h, bins]=hist(obj.endToEndDist(:,2),50); h= h./sum(h);bar(bins,h),...
        hold on, plot(bins,f(obj.rp.numParticles,[bins'],1),'r')
        title('PDF on y')
        
        subplot(2,2,3)
        [h, bins]=hist(obj.endToEndDist(:,3),50); h= h./sum(h); bar(bins,h),...
        hold on, plot(bins,f(obj.rp.numParticles,[bins'],1),'r')
        title('PDF on z')
        %plot(3/(2*pi*obj.numParticles*1)^1.5*exp(3*obj.R.^2/(2*obj.numParticles*1)))
         
    end
    
end
end            
%    [h, bins]=hist(obj.R(:,2),50); h= h./sum(h); figure, bar(bins,h), hold on, plot(bins,f(obj.numParticles,[bins'],1,3),'r')
%     f=@(N,R,b,dim) ((1/(2*pi*N*(b/sqrt(dim))^2))^(dim/2)).*exp(-(sum(R.^2,2))./(2*N*(b/sqrt(dim))^2))
    
    


    