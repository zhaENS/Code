classdef HistoneBis<handle
    % A class for a Histone object
    % histones are attached to chains by default and move on them by 1D
    % diffusion, the position in space is determined by the position of the
    % chain
    % the new chain position is passed to histone, the histone is then
    % projected onto the new position between the vertices it was on in the
    % previous step and then it is moved to its new position

    % TODO: omit the main loop to work on all the coordinates as a matrix 
    % TODO: expand MoveOnPolygon to MoveOnGraph
    
    properties
        curPos
        curPosVertex1
        curPosVertex2
        curPosSlope % position between vertices- used for transformation
        
        prevPos
        prevPosSlope % position between vertices- used for transformation
        prevPosVertex1
        prevPosVertex2
        params = HistoneParams;
                    
    end
    
    methods
        
        function obj = HistoneBis(histoneParams)
            % Class constructor 
            obj.params = histoneParams;                        
        end
        
        function Initialize(obj,initialChainPosition)
            % initialize histone position
            
            % choose a random starting point on the chain
            % choose a bond
            chainPosition = initialChainPosition;
            % choose a random bond
            for hIdx = 1:obj.params.numHistones
                bondNum  = randperm(size(chainPosition,1)-1);
                
                s1 = chainPosition(bondNum(1),:);
                s2 = chainPosition(bondNum(1)+1,:);
                r  = sort(rand(1,2));
                % linearily interpolate the bead position and choose a random
                % point on the bond
                
                obj.prevPos(hIdx,:)        = s1+r(1)*(s2-s1);
                obj.prevPosVertex1(hIdx,:) = bondNum(1);
                obj.prevPosVertex2(hIdx,:) = bondNum(1)+1;
                obj.prevPosSlope(hIdx,:)   = r(1);
                
                obj.curPos(hIdx,:)        = s1+r(2)*(s2-s1);
                obj.curPosVertex1(hIdx,:) = bondNum(1);
                obj.curPosVertex2(hIdx,:) = bondNum(1)+1;
                obj.curPosSlope(hIdx,:)   = r(2);
            end
        end
        
        function Step(obj,chainPosition,dt)
            % move by 1D diffusion
            % project onto new chain position
            chainPos          = chainPosition;
            
            % Update the histone position on the new chain position 
            obj.UpdateHistonePositionOnChain(chainPos);
            fp                 = obj.params.forceParams;
            diffusionForce     = ForceManager.GetDiffusionForce(fp.diffusionForce,obj.curPos,fp.diffusionConst,dt,[],[]);
            if fp.lennardJonesForce
            histoneDist        = ForceManager.GetParticleDistance(obj.curPos);
            lennardJonesForce  = ForceManager.GetLenardJonesForce(fp.lennardJonesForce,obj.curPos,histoneDist,fp.LJPotentialWidth,fp.LJPotentialDepth,[]);
            else
                lennardJonesForce = zeros(obj.params.numHistones,3);
            end
            
            % Advence to a tentative location
            newTempPos         = obj.curPos+ diffusionForce+lennardJonesForce*dt;

            obj.prevPos        = obj.curPos;
            obj.prevPosSlope   = obj.curPosSlope;
            
            for hIdx = 1:obj.params.numHistones
                         
                [obj.curPos(hIdx,:),vert1,vert2,obj.curPosSlope(hIdx)]= ...
                    obj.MoveOnChain(obj.curPos(hIdx,:),newTempPos(hIdx,:),obj.curPosVertex1(hIdx),obj.curPosVertex2(hIdx), chainPos,false);
                                    
                  % update previous
                  obj.prevPosVertex1(hIdx) = obj.curPosVertex1(hIdx);
                  obj.prevPosVertex2(hIdx) = obj.curPosVertex2(hIdx); 
                  
                  % update current
                  obj.curPosVertex1(hIdx) = vert1; 
                  obj.curPosVertex2(hIdx) = vert2; 

            end
        end
        
        function UpdateHistonePositionOnChain(obj,chainPos)
            % Update the position of the histones after the chain has moved

             % updqte the current position on the new chain position
             curDirVec   = chainPos(obj.curPosVertex2,:)-chainPos(obj.curPosVertex1,:);
             obj.curPos  = chainPos(obj.curPosVertex1,:)+bsxfun(@times,obj.curPosSlope,curDirVec);
             % update the previous chain position on the new chain position
             prevDirVec  = chainPos(obj.prevPosVertex2,:)-chainPos(obj.prevPosVertex1,:);
             obj.prevPos = chainPos(obj.prevPosVertex1,:)+bsxfun(@times,obj.prevPosSlope,prevDirVec);
        end
                
        function ParseInputParams(obj,varargin)
            % parse the input parameters
            p = varargin{:};
            if mod(numel(p),2)~=0
                error('name-value pair must be inserted')
            end
            
            for pIdx = 1:numel(p)/2
                obj.params.(p{2*pIdx-1}) = p{2*pIdx};
            end            
        end
    end
    
    methods (Static)
        
        function [vert1,vert2]= FindPointOnPolygon(point,vertices)
            % find a point on a polygon between vert1 and vert2
            flag    = false;
            vertNum = 1;
            numVert = size(vertices,1);
            if numVert <2
                error('polygon must contain at least 2 vertices')
            end
            
            while~flag
                t       = (point - vertices(vertNum,:))./(vertices(vertNum+1,:)-vertices(vertNum,:));
                flag    = all(abs(t-t(1))<1e-10);
                vert1   = vertNum;
                vert2   = vertNum+1;
                vertNum = vertNum+1;
                if vertNum>numVert
                    flag = true;
                    vert1 = [];
                    vert2 = [];
                end
            end
        end
        
        function [newPos,vert1,vert2,posRatio] = MoveOnChain(prevPos,curPos,curPosVertex1, curPosVertex2,vertices,isCircular)%TODO: expand to Move on graph            
            % start from the prevPos between prevPosVertex1 and
            % prevPosVertex2 on the polygon with vertices in vertices
            % find the direction according to curPos-prevPos and relocate
            % the particle position to newPos between vert1 and vert2
            
            % project the curPos-prevPos on the segment of prevPos 
            pPath      = curPos-prevPos;
            % find the direction of motion (toward or away from polygon
            % start)            
            numPolyVert = size(vertices,1);% number of polygon vertices
            motionDir   = sum(pPath.*(vertices(curPosVertex2,:)-vertices(curPosVertex1,:)));
            motionDir   = sign(motionDir);
            pathLength  = sqrt(sum(pPath.^2));
            
            % start from prevPos and subtract the path length from polygon
            % length
            flag = false;
             vert1 = curPosVertex1;
             vert2 = curPosVertex2;
            if motionDir>0
               vert0 = vert2;
            else
               vert0 = vert1;
            end
             cumDist = sqrt(sum((prevPos-vertices(vert0,:)).^2));
             
            while ~flag
               
               pathReminder = pathLength-cumDist;
               flag         = pathReminder<0;
               if ~flag %  move particle to the next segment 
                    
                    vert1 = vert1+motionDir;
                    vert2 = vert2+motionDir;
                    % reflection
                    if (vert2>numPolyVert) || (vert1<1)
                      motionDir = -motionDir;
                      vert1     = vert1+motionDir;
                      vert2     = vert2+motionDir;
                    end                                            
                    cumDist = cumDist+ sqrt(sum((vertices(vert1,:)-vertices(vert2,:)).^2));
               end
            end
                        
            edgeLength = sqrt(sum((vertices(vert2,:)-vertices(vert1,:)).^2));
            if motionDir>0
               posRatio   = (pathReminder+edgeLength)/edgeLength;% the position ratio between vertices
            else
                posRatio  = -pathReminder/edgeLength;
            end
            
            % update the new position 
             newPos     = vertices(vert1,:)+ posRatio*(vertices(vert2,:)-vertices(vert1,:));
        end

    end
    
end