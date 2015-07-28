classdef RouseSimulatorFramework<handle      
    %TODO: add presets for simulations,(load/save)
    %TODO: move the initialization of the chains to allow initialization on
    %      the fly
    %TODO: add a listener that draws emails from the server, to allow
    %      controlling the simulations from a distance
    %TODO: save noise terms seeds to allow re-simulations
    %TODO: save parameters and a short explanation with the results   
    %TODO: make some sense of the parameters passing through the framework.
    %TODO: remove the option to plot online. plotting is always offline.
    %      remove parentAxes property from domainHandlerParams
    %TODO: chain parameters should be initialized before the chains according to the number specified by the simulator params 
    %TODO: build a simple GUI
    %TODO: insert the correct indices for fixed particles in Step 
    %TODO: fix reflection for general meshes shapes
    %TODO: consider moving Domain to the ObjectManager
    % TODO: find a uniform framework for domains and chains
    % TODO: integrate the domain with the object manager.
    % TODO: allow connection of sticky beads
    % TODO: expand ChainParams.m to include parameters for attaching and detaching either from domain or from other chains
    % TODO: complete diffusion constrained on domain
    % TODO: complete diffusion constrained on chains
    % TODO: add documentation
    % TODO: check graph representation for object mapper
    % TODO: allow some bonds to be stiff, find a way to simulate a stiff connector chain efficiently
    % TODO: add reflect outside the domain
    % TODO: add constrained diffusing objects to exert force on the objects
    % TODO: consider joining Domain to objectManager
    properties
        handles     
        objectManager      % data and object manipulation
        dataRecorder       % simulation data module
        simulationGraphics % graphical module
        runSimulation      % flag 
        params             % object parameters
        simulationData  = struct('step',[],'stepTime',[]); % simulation data (should move to dataRecorder)
        batchRound      = 0;
        simulationRound = 0; 
    end
    
    properties (Access=private)
        chainsCenterOfMass
        recipe = struct('PreSimulationBatchActions','',...
                        'PreRunActions','',...
                        'PostRunActions','',...
                        'PostSimulationBatchActions','');
                    
       SetRecipeParamsFlag            = false;
       PreSimulationBatchActionsFlag  = false;
       PreRunActionsFlag              = false;
       PostRunActionsFlag             = false;
       PostSimulationBatchActionsFlag = false;
       
    end    
    
    methods 
        function obj = RouseSimulatorFramework(frameworkParams)
                                
           % Set initial parameters
           obj.OrganizeParams(frameworkParams)

           % Initizlie classes
           obj.InitializeClasses
           
           % set recipe files 
           obj.ReadRecipeFile           
        end
        
        function OrganizeParams(obj,frameworkParams)% obsolete-move to params class
            obj.params = frameworkParams;           
        end
        
        function InitializeClasses(obj)
            
            obj.InitializeDomain;                        
            
            obj.InitializeObjectManager;
                        
            obj.InitializeDataRecorder;
            
            obj.InitializeGraphics;
        end                              
        
        function InitializeDataRecorder(obj)
            obj.dataRecorder = SimulationDataRecorder(obj.params.dataRecorder);
        end
        
        function InitializeGraphics(obj)
          % Initialize framework's graphics 
          obj.simulationGraphics = SimulationFrameworkGraphics(obj);   
          obj.simulationGraphics.CreateControls;
        end
        
        function InitializeDomain(obj)
             % Parameters for domains are set in the organizeParams
            
            obj.handles.classes.domain  = DomainHandler(obj.params.domain);
            
        end   
        
        function InitializeObjectManager(obj)
            % Initialize objectManager
            obj.objectManager = ObjectManager(obj.params.chain);
            obj.objectManager.InitializeObjects(obj.handles.classes.domain);
                                                                               
        end 
        
        function ReadRecipeFile(obj)
            
         t = fileread([obj.params.simulator.recipeFileName,'.rcp']);
         % Search for the function marker 
         [funcStartPos1,funcStartPos2] = regexp(t,'<func>');
         [funcEndPos1,funcEndPos2]     = regexp(t,'</func>');
         for fIdx = 1:numel(funcStartPos1)-1
             % Sort functions into categories
             functionName = strtrim(t(funcStartPos2(fIdx)+1:funcEndPos1(fIdx)-1));
             obj.recipe.(functionName) = t(funcEndPos2(fIdx)+1:funcStartPos1(fIdx+1)-1);
         end
             functionName = strtrim(t(funcStartPos2(end)+1:funcEndPos1(end)-1));
             obj.recipe.(functionName) = t(funcEndPos2(end)+1:end);     
        end               
        
        function Run(obj)
            % Run the simulation according to the specified simulationType
            % parameter
           
            for bIdx = 1:obj.params.simulator.numSimulationBatches
                % perform actions before each simulation batch
                obj.PreSimulationBatchActions
                for sIdx = 1:obj.params.simulator.numSimulations
                    % perform actions before each simulation
                    obj.PreRunActions                    
                    while obj.runSimulation
                        % perform actions before the current step
                        obj.PreStepActions
                        % advance one step of the polymer chain    
                        obj.Step;   
                       % perform action post the current step 
                        obj.PostStepActions                       
                    end
                    % perform actions post simulation
                    obj.PostRunActions
                end
                % perform actions post simulation batch
                obj.PostSimulationBatchActions
            end           
        end         
        
        function PreSimulationBatchActions(obj)
            % Actios performed before a simulation batch      
            obj.batchRound      = obj.batchRound+1;
            obj.simulationRound = 0;
            obj.dataRecorder.CreateNewSimulationBatch
            
            % evaluate the recipe file at the PreSimulationBatchActions
            eval(obj.recipe.PreSimulationBatchActions)
        end
        
        function PreRunActions(obj)
            % Actions called before running the simulation        
            obj.simulationRound = obj.simulationRound+1;
            obj.runSimulation   = true;                        
            eval(obj.recipe.PreRunActions);            
            obj.dataRecorder.NewSimulation(obj.objectManager.handles.chain,obj.params);%(obj.handles.classes.rouse,obj.params);                        
            obj.dataRecorder.SetSimulationStartTime;
            obj.simulationData(obj.batchRound,obj.simulationRound).step     = 1;
            obj.simulationData(obj.batchRound,obj.simulationRound).stepTime = 0;
        end
        
        function PreStepActions(obj)
            eval(obj.recipe.PreStepActions);
        end
        
        function Step(obj,varargin)
             % Next simulation step 
            objList        = 1:obj.objectManager.numObjects;
            dt             = obj.params.simulator.dt;
            % Advance one step and apply object forces                        
           % obj.objectManager.Step(objList,dt)% update current object position
            
            % Update object list
           % objList        = 1:obj.objectManager.numObjects;
            domainInd     = obj.objectManager.GetDomainIndices;% the boundary number the beads belong to 
          
            % Apply domain (global) forces on all objects in the domain 
            %prevParticlePosition     = obj.objectManager.prevPos;      % prev position 
            %curParticlePosition      = obj.objectManager.curPos;       % new pos after internal forces
           % particleDist             = obj.objectManager.particleDist; % distance before applying internal forces
           % fixedParticleNum         = [obj.objectManager.fixedParticles{:}];             
            stickyBeads              = [obj.objectManager.stickyParticles{:}];
           % particlesOnBoundaryAll   = [obj.objectManager.particlesOnBoundary{:}];
           cNb = 1;
           numBeads = [obj.objectManager.objParams.numBeads];
           for oIdx = 1:obj.objectManager.numObjects
            obj.objectManager.Step(oIdx,dt);
            [prevpos,curpos]   = obj.objectManager.GetMembersPosition(oIdx);
            members = obj.objectManager.map.GetObjectMembers(oIdx);
            curParticlePosition = cell(numel(members),1);
            prevParticlePosition = cell(numel(members),1);
            for mIdx = 1:numel(members)
            particlesOnBoundaryAll = [obj.objectManager.objParams(members(mIdx)).beadsOnBoundary];
           
            prevParticlePosition{mIdx} = prevpos{mIdx};
            curParticlePosition{mIdx} = curpos{mIdx};
            particleDist = obj.objectManager.GetMemberDistance(members(mIdx));
            fixedParticleNum = [obj.objectManager.objParams(members(mIdx)).fixedBeadNum];
             
            % Apply external forces from all domains and reflect
            curParticlePosition{mIdx} = obj.handles.classes.domain.Step(prevParticlePosition{mIdx},...
                                    curParticlePosition{mIdx},particleDist,fixedParticleNum,...
                                    particlesOnBoundaryAll,domainInd(cNb:(cNb+numBeads(members(mIdx))-1),1),dt);            
             
             cNb = cNb+numBeads(members(mIdx));
            % diffuse particles on the boundary(currently works only on
            % spheres)
            for dIdx =1:numel(obj.params.domain)
             p    = domainInd(particlesOnBoundaryAll);
             pInd = p==dIdx;
             particleInitPos     = curParticlePosition{mIdx}(particlesOnBoundaryAll(pInd),:); 
             domainRad           = obj.handles.classes.domain.params(dIdx).domainWidth;
             dc                  = obj.handles.classes.domain.params(dIdx).domainCenter;
             %diffusionConst      = obj.handles.classes.domain.params(dIdx).forceParams.diffusionConst;
             diffusionConst      = obj.params.domain(dIdx).forceParams.diffusionConst;
             pIdx                = particlesOnBoundaryAll(pInd);
             pSize               = numel(particlesOnBoundaryAll(pInd));
             for nIdx = 1:pSize
                 s = randperm(50,10);
                 s = max(s);
                 poscurTempo         = DiffusionOnSphere(particleInitPos(nIdx,:),dt,diffusionConst,s,dc,domainRad); 
                 curParticlePosition{mIdx}(pIdx(nIdx),:) = poscurTempo(s,:);
             end
           
           
           %if the fixed beads beads are also on the boundary ,find the
           %positions and correct the position such that all the chains
           %have the same position for the fixed beads;
           fIdxCumul = [];
           fSize = numel(fixedParticleNum);
            for bIdx = 1:fSize
                fIdx      = find(particlesOnBoundaryAll==fixedParticleNum(bIdx));
                fIdxCumul = [fIdxCumul fIdx];
            end
         
            if ~isempty(fIdxCumul)
            curParticlePosition{mIdx}(particlesOnBoundaryAll(fIdxCumul),:) = repmat(curParticlePosition{mIdx}(particlesOnBoundaryAll(fIdxCumul(1)),:),[numel(fIdxCumul),1]);
            obj.objectManager.DealCurrentPosition(members(mIdx),curParticlePosition{mIdx});
            end
            end
            end
           obj.objectManager.DealCurrentPosition(oIdx,reshape([curParticlePosition{:}],numel(members)*64,3)); 
          % obj.objectManager.DealPreviousPosition(oIdx,reshape([prevParticlePosition{:}],numel(members)*64,3));
           end           
             s    = domainInd(stickyBeads);
             sInd = s==dIdx;
             %check if the beads should be sticked
             if ~isempty(stickyBeads(sInd))
                 stickyDistance = obj.params.simulator.encounterDist;
                 %stick the beads
                 obj.objectManager.ConnectStickyParticles(stickyDistance);
                 %update the objList
                %objList = 1:obj.objectManager.numObjects;
             end
            
             %test if the beads in domain can be attached to the boundary
             encounterDistance = obj.params.simulator.encounterDist;
             obj.objectManager.AttachToBoundary(encounterDistance,obj.handles.classes.domain);
            
            % Show simulation
            obj.simulationGraphics.ShowSimulation
            
            % check for interaction between objects   
            if obj.params.simulator.objectInteraction %TODO: move to objectinteractio parameters
             obj.objectManager.ObjectInteraction;
            end
             
            % Update simulation data
            obj.simulationData(obj.batchRound,obj.simulationRound).step = ...
            obj.simulationData(obj.batchRound,obj.simulationRound).step+1;
        end                
                            
                
        function PostStepActions(obj)
            eval(obj.recipe.PostStepActions);
              obj.Record % record data related to the Rouse polymer                        
              % raise the Stop flag if the number of steps is more than the allowed
              obj.runSimulation = obj.runSimulation && ...
             (obj.simulationData(obj.batchRound,obj.simulationRound).step<obj.params.simulator.numSteps);
        end
     
        function PostRunActions(obj)
            % Calculate distributions related to the Rouse polymer
            % chain
            
            eval(obj.recipe.PostRunActions);
            
            obj.dataRecorder.SetSimulationEndTime;
            obj.dataRecorder.SaveResults;% save results
            
            if obj.params.simulator.notifyByEmail
                if mod(obj.simulationRound,obj.params.simulator.notifyCycleLength)==0
                    msg = sprintf('%s%s%s\n%s%s%s\n%s%s\n%s','simulation batch',num2str(obj.batchRound),', ',...
                        'simulation ', num2str(obj.simulationRound),' is done. ',...
                        'Simulation Time: ',num2str(obj.handles.classes.dataRecorder.simulationData(obj.simulationRound).simulationTime),...
                        'Simulating the distance between monomenrs at relexation time, each time bead 1 is connected to bead j, where j=3..64');
                    try
                        SendMail('ofirshukron','Simulation Status',msg);
                    catch
                    end
                end
            end
        end
        
        function PostSimulationBatchActions(obj)
            % actions performed post simulation batch
            eval(obj.recipe.PostSimulationBatchActions);            
            obj.dataRecorder.ClearAllSimulationData;            
%             obj.simulationRound = 0;               
        end                                        
        
        function Record(obj)%TODO: verify data passes correctly 
            % Add samples to the recorded distributions in
            % distributionHandler
            if obj.params.simulator.recordData
               obj.dataRecorder.Add;
            end
        end                                                
        
    end
       
end