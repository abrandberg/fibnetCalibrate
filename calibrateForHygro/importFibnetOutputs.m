function [output,extraOutput] = importFibnetOutputs(PaperDir,netgen,fibnet,simulation,ctrl)

% extraOutput = 0; % Initialize so we always have something to return.

numGrammage = numel(netgen.grammage);
numThickness = numel(netgen.thickness);
numAngle = numel(netgen.angleStd);

    
maxTests = max([numGrammage numThickness numAngle]);
minTests = min([numGrammage numThickness numAngle]);

if numGrammage < maxTests && numGrammage == 1
    netgen.grammage = netgen.grammage*ones(maxTests,1);
end

if numThickness < maxTests && numThickness == 1
    netgen.thickness = netgen.thickness*ones(maxTests,1);
end


counter = 0;

for aLoop = 1:numel(netgen.angleStd)

    for bLoop = 1:numel(simulation)
        counter = counter + 1;

        Project = [simulation(bLoop).type '_AngStd' num2str(netgen.angleStd(aLoop)) simulation(bLoop).dir];
        disp(Project)
        
        while (exist([ctrl.runDir Project '.rea'],'file')~=2)
            pause(10)
        end
%         if not(2==exist([ctrl.runDir Project '.rea'],'file')),continue,end
%         try
            M = importdata([ctrl.runDir Project '.rea'],' ',1);
            M = M.data;
%         catch
%             M = nan(5,5);
%         end
        
        
        
        output(counter).name = Project;
        output(counter).angStd = netgen.angleStd(aLoop);
        if strcmp(simulation(bLoop).type,'Hyg')
            
            output(counter).type =  simulation(bLoop).type;
            output(counter).dir =  '_x';
            output(counter).data = M(end,1);

            counter = counter + 1;
            
            output(counter).name = Project;
            output(counter).angStd = netgen.angleStd(aLoop);
            output(counter).type =  simulation(bLoop).type;
            output(counter).dir =  '_y';
            output(counter).data = M(end,2);
            

        elseif strcmp(simulation(bLoop).type,'E')

            output(counter).type =  simulation(bLoop).type;
            output(counter).dir =  simulation(bLoop).dir;
            output(counter).data = 1e-6*(M(end,3)+M(end,4))*0.5 / (M(end,1)*1e-2 * fibnet.width*1e-6*netgen.thickness(aLoop)*1e-6); % Force divided by strain
            
            
            
            % Here we also put a new structure which should capture the FEM
            % output network data (aka the fiberBuckling library).
            %
            % We mainly use this to check how many bonds we have.
            if Project(end) == 'x' && exist([ctrl.runDir Project '_bondData.csv'],'file')==2
                extraOutput(counter) = importFEMFiberNetwork(Project,ctrl,fibnet,netgen,aLoop);
            end
            
            
            
            
        
        elseif strcmp(simulation(bLoop).type,'FreeCU')
            output(counter).type =  simulation(bLoop).type;
            output(counter).dir =  simulation(bLoop).dir;
            
            % Works on the network submitted to the solver
                            referenceNodes = importPointCloud([ctrl.runDir Project '.xyz']);
                            currentNodes = importPointCloud([ctrl.runDir Project '.xyzEnd']);
    
                            % Derive some representative delta_z values
                            maxX = max(referenceNodes(:,2));
                            selIdx = referenceNodes(:,2) == maxX;
                            deltaZ = mean(currentNodes(selIdx,4));

                            % Estimate the curvatures of the sheet
                            output(counter).data = deltaZ;%[ estimateCurvatureLS(referenceNodes+currentNodes,currentNodes,'2ndOrder-reduced')];
       elseif strcmp(simulation(bLoop).type,'Curl')
            output(counter).type =  simulation(bLoop).type;
            output(counter).dir =  simulation(bLoop).dir;
            
            M = importFODI([ctrl.runDir Project '%nax%.fodi']);
            output(counter).data = M(end,12);
        end % if
        
    end % bLoop
end % aLoop


if exist('extraOutput','var') == 0
    extraOutput = 0;
end



