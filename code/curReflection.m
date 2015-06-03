
%This function is used to dectect the collision between the histones and
%reflect them;
%Parameters:
%numHistones :number of histones;
%dt:time step;
%Ite :number of iterations maximum for each step;
%curPos:the coordinates 3D of histones in each step;
%curPosVertex1(2):vertices left(right);
%chainPos:position of chain for each step;
%posHistone :position of histones in 1D;
%edgeLength:the length of edge of bond;
%edgeLengthTotal:the cumul somme of edgeLength;
%t:the number of times that 2 histones are colliding;
function [curPos,t] = curReflection(numHistones,dt,Ite,curPos,curPosVertex1, curPosVertex2,...
                            chainPos,posHistone,edgeLength,edgeLengthTotal,t,sIdx)   
                        
 
if sIdx >1
    
    %the velocity of each histone;
   v(sIdx,:) = (posHistone(sIdx,:)-posHistone(sIdx-1,:))/dt;  
  
           k    = zeros(numHistones,numHistones);
           ite  = 1;
           cols =zeros(1,numHistones);
      timeRem   = dt*ones(1,numHistones);
    
    while ite <=Ite && min(timeRem)>0
      
       for vIdx = 1:numHistones
      
            for kIdx = 1:numHistones
                 
           %stock all of time k that each 2 histones are colliding for each
           %row;
           k(vIdx,kIdx)         = (posHistone(sIdx-1,vIdx)-posHistone(sIdx-1,kIdx))/(v(sIdx,kIdx)-v(sIdx,vIdx));
         
            end
   
       
           column               = find(k(vIdx,:)>0 & k(vIdx,:)<timeRem(vIdx));
           
           %detect collision
          if ~isempty(column)
            sprintf('collision %s',int2str(t))
            t=t+1;
            
           %find the min time that 2 histones are colliding in each row;
          
            val   = min(k(vIdx,column(1)));
            cols(vIdx)  = column(1);
            for i =1:length(column)
                if k(vIdx,column(i))<val
                val=k(vIdx,column(i));
                cols(vIdx) = column(i);
                end
            end
           
          
            
            %collision point
            posHistoneCol1  = posHistone(sIdx-1,vIdx)+val*v(sIdx,vIdx);
            posHistoneCol2  = posHistone(sIdx-1,cols(vIdx))+val*v(sIdx,cols(vIdx));
         
            %the position of all the histones at the same time when 2 of them are colliding;     
             posHistone(sIdx-1,:)= posHistone(sIdx-1,:)+val*v(sIdx,:);
           
            %update the v after they are colliding;
            timeRem(vIdx)      = timeRem(vIdx)-val;
            v(sIdx,vIdx)       =  (posHistone(sIdx,vIdx)-posHistoneCol1)/timeRem(vIdx);
            v(sIdx,cols(vIdx)) = (posHistone(sIdx,cols(vIdx))-posHistoneCol2)/timeRem(vIdx);
            
       
                                                                
       
            %reflection
            %if 2 collision points have velocity opponent,change their direction ;
           
            if sign(v(sIdx,vIdx))~=sign(v(sIdx,cols(vIdx)))
              posHistone(sIdx-1,vIdx) =  posHistoneCol1-timeRem(vIdx)*v(sIdx,vIdx);
              posHistone(sIdx-1,cols(vIdx)) =  posHistoneCol2-timeRem(vIdx)*v(sIdx,cols(vIdx));            
             else if abs(v(sIdx,vIdx)) > abs(v(sIdx,cols(vIdx)))
             %if not,change the direction of velocity which is
             %faster than the other
                posHistone(sIdx-1,vIdx)       =  posHistoneCol1-timeRem(vIdx)*v(sIdx,vIdx);
                posHistone(sIdx-1,cols(vIdx)) =  posHistoneCol2+timeRem(vIdx)*v(sIdx,cols(vIdx));
                 else
                posHistone(sIdx-1,vIdx)       = posHistoneCol1+timeRem(vIdx)*v(sIdx,vIdx);
                posHistone(sIdx-1,cols(vIdx)) = posHistoneCol2-timeRem(vIdx)*v(sIdx,cols(vIdx));
                 end
            end           
            
          end
       
       end
  
 %transform the coordinate into 3D;
     curPos(:,:,sIdx-1) =  bsxfun(@times,((posHistone(sIdx-1,:)-edgeLengthTotal(sIdx-1,curPosVertex1(sIdx-1,:)))...
                                        ./edgeLength(sIdx,curPosVertex1(sIdx-1,:)))',...
                                        (chainPos(curPosVertex2(sIdx-1,:),:)-chainPos(curPosVertex1(sIdx-1,:),:)))...
                                        +chainPos(curPosVertex1(sIdx-1,:),:);
                                    
                                  
       ite = ite+1;
    end
    

  

   
end                                   

end


