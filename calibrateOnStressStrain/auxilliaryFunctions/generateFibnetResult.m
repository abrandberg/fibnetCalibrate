function output = generateFibnetResult(x,ctrl,netgen,fibnetIn,paramName)

% Prevent negative values
x(x<0.1) = 0.1;


% Now wrap a loop around it!
for tLoop = 1:numel(fibnetIn)
    
    % Extract this loop's fibnet data
    fibnet = fibnetIn(tLoop);   
    
    % Generate network/check that network exists
    [PaperDir,FName] = generateNetworks(netgen,ctrl);
    
    fibnet.INPF = ["'" PaperDir{:} '/' FName "'"];
    
    % Construct input file
    constructSolverInput(x,fibnet,ctrl,paramName)

    
    % Solver call
    % Specify ANSYS call according to your installation
    
    % Run input file
    system(['SET KMP_STACKSIZE=4096k & "C:\Program Files\ANSYS Inc\v150\ANSYS\bin\winx64\ansys150.exe"  -p aa_t_a -np 4 -dir "' ctrl.runDir '" -j "' ctrl.simName '" -s read -l en-us -b -i "' ctrl.runDir ctrl.simName '.dat" -o "' ctrl.runDir ctrl.simName '.out" -custom "' ctrl.customExecutable '" ']);

    % Import output
    importedREA = importfile_REA([ctrl.runDir filesep ctrl.simName '.rea']);

    if strcmp(ctrl.optiVar,'Stiffness')
        % Calculate the elastic modulus
        selIdxTS = [1 size(importedREA,1)];
        output(tLoop)   = 1e-6*diff(0.5*(importedREA(selIdxTS,2)+importedREA(selIdxTS,3)))./ ...
                         (diff(importedREA(selIdxTS,1))/100*fibnet.widthIn*1e-6*netgen.thickness*1e-6);
                     
    elseif strcmp(ctrl.optiVar,'Strength')
        
        if max(0.5*(importedREA(:,2)+importedREA(:,3))) == max(0.5*(importedREA(end,2)+importedREA(end,3)))
            disp('No failure!')
            output(tLoop) = 1e9;
        else
            output(tLoop) = 1e-6*max(0.5*(importedREA(:,2)+importedREA(:,3)))/(fibnet.widthIn*1e-6*netgen.thickness*1e-6);
        end
    end
end





