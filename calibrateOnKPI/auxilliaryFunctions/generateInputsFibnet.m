function [costFcn] = generateInputsFibnet(x,fixedParams)

% clear all; close all; clc
ctrl.workingDir = cd;
% addpath(fullfile(ctrl.workingDir,'auxilliaryFunctions'))

runDir = 'optiDir_v2';
% mkdir('optiDir')
% disp(['Current factor ',num2str(Efiber)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bondOptions = fixedParams.debonding;%[0];
% Options: 0 - No delamination, 1 - Delamination

plastOptions = fixedParams.plasticity;%[0];                                                         
% Options: 0 - No plasticity,   1 - Plasticity

keepString = 'file1_L8.0_W8.0_g60.0';

for zLoop = 1:1
    inputOptions{zLoop} = strrep(keepString,'file1',horzcat('file',num2str(zLoop)));
end
% disp(['Size of inputOptions: ' num2str(numel(inputOptions))])
            
% Network - no file extension
% CAVLIN AND FELLERS SET AKA NETWORK ALPHA: 'file1_L30.0_W6.0_g300.0'
% STOCKMANN SET AKA NETWORK BETA:           

% Network width, in [um]
widthOptions = fixedParams.width;%2e3;
% Network length, in [um]
lengthOptions = fixedParams.length;%2e3;

% Options 'C' - Compression, 'T' - Tension
loadOptions = {'T'};

hingeOptions = [0];%0;%[0:0.02:0.1 0.15 0.2:0.1:1.0 1.5 2:1:8];
% Hinge: Number of hinges per mm fiber times 0.025.

clampOptions = [0];                %[0 -0.1 -0.2 -0.3 -0.4]                                          
% Clamping: E_z strain to be applied at clamps

shearOptions = [0];            % [0 1]                                              
% Options: 0 - No correction, 
%          1 - Worst case (90% reduction)
%          2 - 50% reduction
%          3 - Realistic 

materialOptions = [99];            % [1:5]                                             
% Options: 1 - Normal (CTMP),        2 - Better plastic properties
%          3 - Higher fiber modulus, 4 - Stronger bonds
%          5 - Stiffer bonds         6 - Custom
materialMultiplier = [1];

loadDirectionOptions = [1];
% Options: 1 - MD (along x)
%          2 - CD (along y)


bondStrengthScaling = [1];
% Options: 0 - Bond strength indenpendent of underlying geometry
% Options: 1 - Bond strength \propto area in contact



% Figure out what we are keeping constant.
if strcmp(fixedParams.variable,'Efiber')
    Ebond = fixedParams.Ebond;
    Efiber = x;
    PfiberS = fixedParams.Pfiber;
    PfiberH = fixedParams.Pfiber;
    Sbond = fixedParams.Sbond;
    
elseif strcmp(fixedParams.variable,'Ebond')
    Ebond = x;
    Efiber = fixedParams.Efiber;
    PfiberS = fixedParams.Pfiber;
    PfiberH = fixedParams.Pfiber;
    Sbond = fixedParams.Sbond;
    
elseif strcmp(fixedParams.variable,'Pfiber')
    Ebond = fixedParams.Ebond;
    Efiber = fixedParams.Efiber;
    PfiberS = x;
    PfiberH = x;
    Sbond = fixedParams.Sbond;
    
elseif strcmp(fixedParams.variable,'Sbond')
    Ebond = fixedParams.Ebond;
    Efiber = fixedParams.Efiber;
    PfiberS = fixedParams.Pfiber;
    PfiberH = fixedParams.Pfiber;
    Sbond = x;

elseif strcmp(fixedParams.variable,'Demo')
    Ebond = fixedParams.Ebond;
    Efiber = fixedParams.Efiber;
    PfiberS = fixedParams.PfiberS;
    PfiberH = fixedParams.PfiberH;
    Sbond = fixedParams.Sbond;
    
elseif strcmp(fixedParams.variable,'Nonlin')
    Ebond = fixedParams.Ebond;
    Efiber = fixedParams.Efiber;
    
    if x(1) < 0
        disp('Warning: Yield stress bottoming out')
    end
    if x(2) < 0
        disp('Warning: Hardening modulus bottoming out')
    end
    if x(3) < 0
        disp('Warning: Bond strength bottoming out')
    end
    PfiberS = max(0,x(1));
    PfiberH = max(0,x(2));
    Sbond = max(0,x(3));

    elseif strcmp(fixedParams.variable,'Nonlin2')
%     Ebond = fixedParams.Ebond;
    Efiber = fixedParams.Efiber;
    
    if x(1) < 0
        disp('Warning: Yield stress bottoming out')
    end
    if x(2) < 0
        disp('Warning: Hardening modulus bottoming out')
    end
    if x(3) < 0
        disp('Warning: Bond strength bottoming out')
    end
    Ebond = max(0,x(4));
    PfiberS = max(0,x(1));
    PfiberH = max(0,x(2));
    Sbond = max(0,x(3));
    
else
%     Efiber = fixedParams.Efiber;
    disp(stop)
end




% Execution environment
execEnvir = 'Fred'; % Tensor, Bertil, Kebnekaise

if strcmp(execEnvir,'Kebnekaise')
    shBaseFile = 'sbatchFile2.sh';
    pDirectory = '/psearch,/pfs/nobackup/home/a/augustbr/April2018/MACS';
elseif strcmp(execEnvir,'Tensor')
    shBaseFile = 'sbatchTensor.sh';
    pDirectory = '/psearch,/scratch/USERS/AUGUST/fibnetRuns_J2018_2/MACS';
elseif strcmp(execEnvir,'Bertil')
    shBaseFile = 'sbatchBertil.sh';
    pDirectory = '/psearch,/scratch/users/august/fibnet/MACS';
else
    shBaseFile = 'sbatchBertil.sh'; % Dummy
    pDirectory = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            FILE GENERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for aLoop = 1:length(bondOptions)
for bLoop = 1:length(plastOptions)
for cLoop = 1:length(inputOptions)
for dLoop = 1:length(widthOptions)
for eLoop = 1:length(loadOptions)
for fLoop = 1:length(lengthOptions)
for gLoop = 1:length(hingeOptions)
for hLoop = 1:length(clampOptions)
for iLoop = 1:length(shearOptions)
for jLoop = 1:length(materialOptions)
for kLoop = 1:length(loadDirectionOptions)
    for lLoop = 1:length(materialMultiplier)
   
    debondingSwitch     = bondOptions(aLoop);
    plasticitySwitch    = plastOptions(bLoop);
    inputFile           = inputOptions{cLoop};
    realWidth           = widthOptions(dLoop);
    loadDirection       = loadOptions{eLoop};         
    realLength          = lengthOptions(fLoop);
    nHin                = hingeOptions(gLoop);%*(realLength*realWidth)*1e-6; % hinges per square mm
    clam                = clampOptions(hLoop);
    shearCorr           = shearOptions(iLoop);
    matChoice           = materialOptions(jLoop);
    loadMDCD            = loadDirectionOptions(kLoop);
    matMult           = materialMultiplier(lLoop);
    

    % Import the submission outline
    shSkeleton = shBaseFile;%'sbatchFile2.sh';
    fileID = fopen(shSkeleton,'r');
    importedText = fread(fileID,'*char');
    fclose(fileID);

%     if hingeOptions(gLoop) == 0
%         hingeNameTemp = 0;
%     else 
%         hingeNameTemp = log10(hingeOptions(gLoop));
%     end
    
    if clampOptions(hLoop) > 0 || clampOptions(hLoop) < 0
        lengthTemp = lengthOptions(fLoop)+350;
    else
        lengthTemp = lengthOptions(fLoop);
    end

 simName = horzcat('delmeOpt');
 %   simName = horzcat('M',num2str(matChoice),sprintf('%3.2f',matMult),  ...
 %                     '_L',num2str(realLength), ...
 %                     '_W',num2str(realWidth),  ...
 %                     '_',strtok(inputFile(5:end),'_'),loadDirection,num2str(loadMDCD),        ...
 %                     '_DP',num2str(debondingSwitch),num2str(plasticitySwitch), ...
 %                     '_HS',num2str(nHin), ...
 %                            num2str(shearOptions(iLoop)), ...
 %                     '_C',sprintf('%03.2f',clampOptions(hLoop)));
%     if strcmp(execEnvir,'Kebnekaise')
%         finalLine = horzcat('srun -n 1 -c 14 /pfs/nobackup/home/a/artkula/August/ansyscust.e150 -i ',simName,'.dat -p aa_r -j ',simName,' -o ',simName,'.out -s read -b -np 14');
%     elseif strcmp(execEnvir,'Tensor')
%         finalLine = horzcat('/usr/ansys_inc/v150/ansys/bin/ansys150 -p aa_r -np 8  -j ',simName,' -s read -l en-us -b -i ',simName,'.dat -o ',simName,'.out -custom "/home/august/FibNet/Heavy/ansyscust.e150"');
%     elseif strcmp(execEnvir,'Bertil')
%         finalLine = horzcat('/usr/ansys_inc/v150/ansys/bin/ansys150 -p aa_r -np 8  -j ',simName,' -s read -l en-us -b -i ',simName,'.dat -o ',simName,'.out -custom "/scratch/users/artem/FibreNet/FibNet20190227/ansyscust.e150"');
%     end
    if strcmp(execEnvir,'Kebnekaise')
        finalLine = horzcat('srun -n 1 -c 8 /pfs/nobackup/home/a/artkula/src/FibNet20190227/ansyscust.e150 -i ',simName,'.dat -p aa_r -j ',simName,' -o ',simName,'.out -s read -b -np 2');
    elseif strcmp(execEnvir,'Tensor')
        %/usr/ansys_inc/v150/ansys/bin/ansys150 -p aa_r -np 8  -j M1_L700_W6000_1C_DP11_HS00_C0.00 -s read -l en-us -b -i M1_L700_W6000_1C_DP11_HS00_C0.00.dat -o M1_L700_W6000_1C_DP11_HS00_C0.00.out -custom "/home/artem/src/FibNet20190227/ansyscust.e150"
        finalLine = horzcat('/usr/ansys_inc/v150/ansys/bin/ansys150 -p aa_r -np 2  -j ',simName,' -s read -l en-us -b -i ',simName,'.dat -o ',simName,'.out -custom "/home/artem/src/FibNet20190227/ansyscust.e150"');
        %finalLine = horzcat('/usr/ansys_inc/v150/ansys/bin/ansys150 -p aa_t_a -np 2  -j ',simName,' -s read -l en-us -b -i ',simName,'.dat -o ',simName,'.out -custom "/home/artem/src/FibNet20190227/ansyscust.e150"');
    elseif strcmp(execEnvir,'Bertil')
%         finalLine = horzcat('/usr/ansys_inc/v150/ansys/bin/ansys150 -p aa_t_a -np 2  -j ',simName,' -s read -l en-us -b -i ',simName,'.dat -o ',simName,'.out -custom "/scratch/users/artem/FibreNet/FibNet20190227/ansyscust.e150" > ansysOut.txt');
        finalLine = horzcat('/usr/ansys_inc/v150/ansys/bin/ansys150 -p aa_r -np 1  -j ',simName,' -s read -l en-us -b -i ',simName,'.dat -o ',simName,'.out -custom "/scratch/users/artem/FibreNet/FibNet20200319/ansyscust.e150" > ansysOut.txt');
    else
        finalLine = horzcat('SET KMP_STACKSIZE=4096k &  "C:\Program Files\ANSYS Inc\v150\ANSYS\bin\winx64\ansys150.exe" -p aa_r -np 1  -j ',simName,' -s read -l en-us -b -i ',simName,'.dat -o ',simName,'.out -custom "c:\Program Files\ANSYS Inc\v150\ansys\custom\user\Mossab_drying\ansys.exe"');
    end
    % Write real submission file
    fileID = fopen(horzcat(simName,'.sh'),'w');
    fprintf(fileID,'%s\n',importedText,finalLine);
    fclose(fileID);

    % Write the simulation file
    % Import the submission outline
    part1SimFile = 'part1Sim.txt';
    fileID = fopen(part1SimFile,'r');
    part1Sim = fread(fileID,'*char');
    fclose(fileID);
    
%     part2SimFile = 'part2Sim_v4.txt';
    part2SimFile = 'part2Sim_v5.txt';
    fileID = fopen(part2SimFile,'r');
    part2Sim = fread(fileID,'*char');
    fclose(fileID);

    inputF          = horzcat('INPF = ''',inputFile,'''');
    debo            = horzcat('debo = ',num2str(debondingSwitch),'                 ! 0: no debonding');
    plast           = horzcat('plast = ',num2str(plasticitySwitch),'                ! 0: no plasticity ');
    realWidth       = horzcat('RWID = ',num2str(realWidth));
    realLength      = horzcat('RLID = ',num2str(lengthTemp)); % OBS adapted for clamps
    numHinges       = horzcat('hinge = ',num2str(nHin));
    percentClamped  = horzcat('clam = ',num2str(clam));
    shearCorrMode   = horzcat('scor = ',num2str(shearCorr));
    materialMode    = horzcat('msch = ',num2str(matChoice));
    loadMode        = horzcat('LDIR = ',num2str(loadMDCD));
    multiplierMode  = horzcat('matM = ',num2str(matMult));
    fiberModulus = horzcat('efib = ',num2str(Efiber));
    bondModulus = horzcat('ebon = ',num2str(Ebond));
    bondStrength = horzcat('sbon = ',num2str(Sbond));
    pModulusS = horzcat('pfib1 = ',num2str(PfiberS));
    pModulusH = horzcat('pfib2 = ',num2str(PfiberH));
    
    appliedStrain = horzcat('fs = ',num2str(fixedParams.appliedStrain));
    numSteps = horzcat('nstep = ',num2str(fixedParams.numSteps));
    
    bssSetting = horzcat('BSS = ',num2str(bondStrengthScaling));

    if strcmp(loadDirection,'T')
        signOfLoad = horzcat('SOL = 1');
    elseif strcmp(loadDirection,'C')
        signOfLoad = horzcat('SOL = -1');
    end

    % Write the real simulation file
    fileID = fopen(horzcat(simName,'.dat'),'w');
    fprintf(fileID,'%s\n',part1Sim,inputF,debo,plast,realWidth,realLength,     ...
                          numHinges,percentClamped,shearCorrMode,materialMode,multiplierMode, ...
                          loadMode,signOfLoad,fiberModulus,bondModulus,bondStrength,pModulusS,pModulusH,appliedStrain,numSteps,bssSetting,pDirectory,part2Sim);
    fclose(fileID);
    
    % CD to directory and submit
simFile = horzcat(simName,'.dat');
subFile = horzcat(simName,'.sh');
movefile(simFile,horzcat(runDir,filesep,simFile))
movefile(subFile,horzcat(runDir,filesep,subFile))

cd(runDir)
% runString = horzcat('cd ',runDir,'; qsub ',subFile,'; cd ..');
%system(runString);
% disp('Submitting the job')
% system(horzcat('cd ',runDir,'; ',finalLine,'; cd ..'));
% disp('Exited job')
system(finalLine);

% switchVar = 1;
% tic
% while switchVar
%    pause(60)
%    if isfile(horzcat(simName,'.rea'))
%         switchVar = 0;
%    else
%        toc
% %        disp('still waiting')
%    end
% end
cd(ctrl.workingDir)

    end
end
end
end
end
end
end
end
end
end
end
end

save('efib.mat','Efiber')



% Import results and evaluate the stiffness
% disp('Reading results')
resuArray = importfile_REA(horzcat(runDir,filesep,simName,'.rea'));
pause(10)



modulusToFit = 1610.2; % Nm/g
targetTEA = 60*0.04685;
targetStrength = 9.60;
targetFailureStrain = 0.8171;

if strcmp(fixedParams.variable(1),'E') % Finding stiffness
    % selIdxTS = [2 min(8,size(resuArray,1))];
    selIdxTS = [1 size(resuArray,1)];

    tensileStiffness = 1e-6*diff(0.5*(resuArray(selIdxTS,2)+resuArray(selIdxTS,3)))./(diff(resuArray(selIdxTS,1))/100)/(widthOptions*1e-6);

    costFcn = sqrt((tensileStiffness - modulusToFit*60)^2)./(modulusToFit*60);
    
    subplot(1,2,1)
    plot(resuArray(:,1),resuArray(:,2),'m')
    hold on
%     xline(targetFailureStrain);
%     yline(1e6*targetTensileStrength*widthOptions*1e-6);
    plot([0 0.2], widthOptions*modulusToFit*60.*[0 0.2]/100,'k-.')

elseif strcmp(fixedParams.variable(1),'P')
    
    
    breakStrain = targetFailureStrain;
    selIdxTEA = resuArray(:,1) <= breakStrain;
    tensileEnergyAbsorption = 1e-6*1e-2*trapz(resuArray(selIdxTEA,1),resuArray(selIdxTEA,2))./(widthOptions*1e-6);
%     targetTEA = 60*0.05166;
    costFcn = sqrt((tensileEnergyAbsorption - targetTEA)^2)/(targetTEA);
    
    subplot(1,2,1)
    plot(resuArray(:,1),resuArray(:,2),'b')
    xline(targetFailureStrain);
    hold on
    subplot(1,2,2)
    plot(resuArray(:,1),1e-6*1e-2*cumtrapz(resuArray(:,1),resuArray(:,2))./(widthOptions*1e-6),'r')
    xline(targetFailureStrain);
    yline(targetTEA);
    hold on
    pause(0.5)
elseif strcmp(fixedParams.variable(1),'S') % Finding strength
    [maxF,maxIdx] = max(resuArray(:,2));
    if maxIdx == size(resuArray,1)
        disp(stop)
    end
    tensileStrength =  1e-6*maxF/(widthOptions*1e-6);
    breakStrain = resuArray(maxIdx,1);
    costFcn1 = sqrt((tensileStrength - targetStrength*60)^2)/(targetStrength*60);
    costFcn2 = sqrt((breakStrain - targetFailureStrain)^2)/(targetFailureStrain);
    
    costFcn = costFcn1 + costFcn2;
    
    fprintf('%20s %20s \n','','x(1)');
    fprintf('%20s %20.4f\n','',x(1));
    fprintf('%20s %20s %20s %20s %20s\n','','Tensile error','Strain error','costFcn');
    fprintf('%20s %20.4f %20.4f %20.4f %20.4f\n','',costFcn1,costFcn2,costFcn);
    
    
    targetTensileStrength = targetStrength*60;
     
    subplot(1,3,1)
    plot(resuArray(:,1),resuArray(:,2),'b')
    hold on
    xline(targetFailureStrain);
    yline(1e6*targetTensileStrength*widthOptions*1e-6);
    plot([0 1], widthOptions*modulusToFit*60.*[0 0.01],'k-.')
    
    subplot(1,3,2)
    plot(resuArray(:,1),1e-6*1e-2*cumtrapz(resuArray(:,1),resuArray(:,2))./(widthOptions*1e-6),'r')
    hold on
    xline(targetFailureStrain);
%     targetTEA = targetTEA;
    yline(targetTEA);
    plot(resuArray(:,1),1e-2*1e-2*cumtrapz(resuArray(:,1),resuArray(:,1).*modulusToFit*60),'k-.')
    
    subplot(1,3,3)
    plot(x,costFcn1,'ob')
    hold on
    plot(x,costFcn2,'sr')
    plot(x,costFcn,'+k')
    
elseif strcmp(fixedParams.variable(1),'N')
    
    [maxF,maxIdx] = max(resuArray(:,2));
    
    breakStrain = resuArray(maxIdx,1);
    selIdxTEA = resuArray(:,1) <= breakStrain;
    
    
    
    tensileStrength =  1e-6*maxF/(widthOptions*1e-6);
    tensileEnergyAbsorption = 1e-6*1e-2*trapz(resuArray(selIdxTEA,1),resuArray(selIdxTEA,2))./(widthOptions*1e-6);
%     tensileStrain = breakStrain;
%     targetTEA = 60*0.05166;
%     targetTensileStrength = targetStrength*60;
    
%     costFcn(1) = sqrt((tensileStrength - targetStrength*60)^2)/(targetStrength*60);
%     costFcn(2) = sqrt((breakStrain - targetFailureStrain)^2)/(targetFailureStrain);
%     costFcn(3) = sqrt((tensileEnergyAbsorption - targetTEA)^2)/(targetTEA);
    
%     costFcn(1) = ((tensileStrength - targetStrength*60));%/(targetStrength*60);
%     costFcn(2) = ((breakStrain - targetFailureStrain));%/(targetFailureStrain);
%     costFcn(3) = ((tensileEnergyAbsorption - targetTEA));%/(targetTEA);
    
    selIdxTS = [1 4];
    tensileStiffness = 1e-6*diff(0.5*(resuArray(selIdxTS,2)+resuArray(selIdxTS,3)))./(diff(resuArray(selIdxTS,1))/100)/(widthOptions*1e-6);
    % 1770.8*60
    term4 = sqrt((tensileStiffness - modulusToFit*60)^2)./(modulusToFit*60);

    term1 =  ((tensileStrength - targetStrength*60)^1)/(targetStrength*60);
    term2 = ((breakStrain - targetFailureStrain)^1)/(targetFailureStrain);
    term3 = ((tensileEnergyAbsorption - targetTEA)^1)/(targetTEA);
%     costFcn = term1 + term2 + term3;
    costFcn = term1 + term2 + term3;
    
    costFcn(1) = term1;
    costFcn(2) = term2;
    costFcn(3) = 0;%term3;
    
    fprintf('%20s %20s %20s %20s\n','','x(1)','x(2)','x(3)');
    fprintf('%20s %20.4f %20.4f %20.4f\n','',x(1),x(2),x(3));
    fprintf('%20s %20s %20s %20s %20s\n','','Tensile error','Strain error','Energy error','costFcn');
    fprintf('%20s %20.4f %20.4f %20.4f %20.4f\n','',term1,term2,term3,sum(costFcn));
%     fprintf('%20s %20s %20s %20s %20s\n','','x(1)','x(2)','x(3)','x(4)');
%     fprintf('%20s %20.4f %20.4f %20.4f %20.4f\n','',x(1),x(2),x(3),x(4));
%     fprintf('%20s %20s %20s %20s %20s %20s\n','','Tensile error','Strain error','Energy error','Stiff. error','costFcn');
%     fprintf('%20s %20.4f %20.4f %20.4f %20.4f %20.4f\n','',term1,term2,term3,term4,costFcn);

    subplot(2,2,1)
    plot(resuArray(:,1),resuArray(:,2),'b')
    hold on
    xline(targetFailureStrain);
    yline(1e6*targetStrength*60*widthOptions*1e-6);
    plot([0 1], widthOptions*modulusToFit*60.*[0 0.01],'k-.')
    
    subplot(2,2,2)
    plot(resuArray(:,1),1e-6*1e-2*cumtrapz(resuArray(:,1),resuArray(:,2))./(widthOptions*1e-6),'r')
    hold on
    xline(targetFailureStrain);
    yline(targetTEA);
    plot(resuArray(:,1),1e-2*1e-2*cumtrapz(resuArray(:,1),resuArray(:,1).*modulusToFit*60),'k-.')
    
    subplot(2,2,3)
    plot([1 2 3],x,'ob')
    hold on
    
    subplot(2,2,4)
    plot([1 2 3],[term1 term2 term3],'sr')
    hold on
    
    
    
    pause(0.5)
    
    if maxIdx == size(resuArray,1)
%         dbstop if error
        disp('Warning, failure is endpoint')
    end
    
elseif strcmp(fixedParams.variable(1),'D') % Demo
    
    targetTensileStrength = targetStrength*60;
    targetTEA = 60*0.05166;
    
    figure;
    subplot(1,2,1)
    plot(resuArray(:,1),resuArray(:,2),'b')
    hold on
    xline(targetFailureStrain);
    yline(1e6*targetTensileStrength*widthOptions*1e-6);
    plot([0 1], widthOptions*modulusToFit*60.*[0 0.01],'k-.')
    
    subplot(1,2,2)
    plot(resuArray(:,1),1e-6*1e-2*cumtrapz(resuArray(:,1),resuArray(:,2))./(widthOptions*1e-6),'r')
    hold on
    xline(targetFailureStrain);
    yline(targetTEA);
    plot(resuArray(:,1),1e-2*1e-2*cumtrapz(resuArray(:,1),resuArray(:,1).*modulusToFit*60),'k-.')    
    pause(0.5)
    
end

% disp('Deleting outputs')
delete(horzcat(runDir,filesep,simName,'.*'))
delete(horzcat(runDir,filesep,'*.png'))
pause(10)










% fprintf('%20s %20s %20s \n','Current','Target','costFcn');
% fprintf('%20.1f %20.1f %20.5f \n',tensileStiffness,1770.8*60,costFcn);

% disp('End of optiFcn')
