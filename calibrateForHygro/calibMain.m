%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calibration script for the fibnet hygroexpansion
%
%
% created by: August Brandberg
% date: 2021-05-10

% Meta-instructions
clear; close all; clc;
format compact

% Set up ctrl structure
ctrl.workDir = cd;
% ctrl.runDir = [ctrl.workDir filesep 'runs1\'];
ctrl.saveDir = [ctrl.workDir filesep 'saveme\'];


ctrl.interpreter = 'latex';
ctrl.histogramInstructions = {'normalization','probability'};

ctrl.colors = [255 94 105 ;
               0 188 213]./255;
ctrl.lineInstructions = {'markersize',8,'linewidth',1.25};

ctrl.parallelNetGetRuns = 12;
ctrl.parallelFibNetRuns = 6;
ctrl.fibnetNP = 4;                  % Number of cores per FibNet run

ctrl.executable = 'C:\Program Files\ANSYS Inc\v150\ansys\custom\user\Hygro3D_Kurosh\ansys.exe';
ctrl.pulpDataFile = 'Euca2.txt'%'Euca2.txt';

% assert(strcmp(ctrl.runDir(end),filesep),'ctrl.runDir is missing trailing file separator.')

% Generate necessary folders
% if not(exist(ctrl.runDir,'dir'))
%     mkdir(ctrl.runDir)
% end
if not(exist(ctrl.saveDir,'dir'))
    mkdir(ctrl.saveDir)
end

pulpCharacterizationRaw = characterizeInputOFA(ctrl);


netgen.length = 8;                                          % Length along MD, cm
netgen.width = 8;                                           % Length along CD, cm
netgen.grammage = 80;                                       % Surface weight per area, g/m^2
netgen.thickness = 120;%180;                                     % Sheet thickness, um
netgen.angleStd = [linspace(30,80,12)];
netgen.interfaceAngle = 9;
fibnet.length = 4000;
fibnet.width = 4000;

materialProps.EX = 29e3; % 21e3 good
materialProps.EZ = 17e3; %12e3;
materialProps.Kct = 8000;
materialProps.searchRadiusMultiplier = 1.00;

% Here we could generate a master curve of fiber orientation anisotropies
densitiesToTry = [400 600 800]; % kg/m^3
resultingThicknessToTry = 1e6*(1e-3*netgen.grammage./densitiesToTry);


if exist('surrogateModelFiberAnisotropy_v2.mat','file')
    load('surrogateModelFiberAnisotropy_v2.mat')
    ctrl.runDir = [ctrl.workDir filesep 'fiberOrientationCalib_1' filesep];
else    
    figure;
    for aLoop = 1:length(resultingThicknessToTry)

        ctrl.runDir = [ctrl.workDir filesep 'fiberOrientationCalib_' num2str(aLoop) filesep];

        netgen.thickness = resultingThicknessToTry(aLoop);
        [XangleStd,YTSI_Ratio] = fiberOrientationAnisotropySurrogateModel(netgen,fibnet,ctrl,materialProps);

        plot(XangleStd,YTSI_Ratio,'-o')
        hold on

        fiberOrientationTable(aLoop).density = densitiesToTry(aLoop);
        fiberOrientationTable(aLoop).thickness = netgen.thickness;
        fiberOrientationTable(aLoop).XangleStd = XangleStd;
        fiberOrientationTable(aLoop).YTSI_Ratio = YTSI_Ratio;
    end
    xlabel('netgen.angleStd [deg]')
    ylabel('Ratio of elastic moduli')

    save('surrogateModelFiberAnisotropy_v2.mat','fiberOrientationTable','-v7.3') % ,'netgenRef','fibnetRef','ctrlRef'
end
% plot(XangleStd,[pulpCharacterizationGenerated.FORatio],'-s')


% Check the pulp characterization from the files
% 1. Check element data
% 2. Check what was submitted to the solver
pulpCharacterizationGenerated = characterizeFEMInputs(ctrl,fibnet);



figure;
tiledlayout(1,4,'padding','none')
nexttile;
xVec = linspace(0,10e-3,50);
histogram([pulpCharacterizationRaw.lc]./[pulpCharacterizationRaw.curl + 1]*1e-3,xVec,'normalization','probability')
hold on
histogram(pulpCharacterizationGenerated(1).lp,xVec,'normalization','probability')

nexttile;
xVec = linspace(0,100e-6,50);
histogram(1.3*pulpCharacterizationRaw.width*1e-6,xVec,'normalization','probability')
hold on
histogram(pulpCharacterizationGenerated(1).width,xVec,'normalization','probability')

nexttile;
xVec = linspace(0,50e-6,50);
histogram(pulpCharacterizationRaw.wallTkn*1e-6,xVec,'normalization','probability')
hold on
histogram(pulpCharacterizationGenerated(1).wallTkn,xVec,'normalization','probability')

nexttile;
xVec = linspace(0,1.5,50);
histogram(pulpCharacterizationRaw.curl,xVec,'normalization','probability')
hold on
histogram(pulpCharacterizationGenerated(1).curl,xVec,'normalization','probability')


% plot(XangleStd,[pulpCharacterizationGenerated.a]./[pulpCharacterizationGenerated.b],'-d')
% plot([fiberOrientationTable(1).XangleStd],[pulpCharacterizationGenerated.a]./[pulpCharacterizationGenerated.b],'-d')
% xlabel('netgen.angleStd [deg]')
% ylabel('Ratio of elastic moduli')
% legend('Mechanical simulation','Fiber seg. orientation','location','northeast')
% pause(0.5)




%%%
X = [reshape(repmat([fiberOrientationTable.density],numel(netgen.angleStd),1),1,[])]';
Y = [fiberOrientationTable.YTSI_Ratio]';
V = [fiberOrientationTable.XangleStd]';
VFcn = scatteredInterpolant(X,Y,V);

% figure;
% plot3(X,Y,V,'.')

[xq,yq] = meshgrid([400:50:800],[1:0.1:3]);
vq1 = VFcn(xq,yq);
% mesh(xq,yq,vq1)

% Vq = VFcn(76.66*1e-3/(179.33*1e-6),1.44)
%,
%%%

ctrl.runDir = [ctrl.workDir filesep 'calib_EAndBeta\'];
netgen.grammage = [76.66 72.51 84.11];
netgen.thickness = [179.33 169.83 194.46];
TSI_RatioToMatch = [1.44 1.63 1.78];

averageDensity = 1e-3*[netgen.grammage]./(1e-6*[netgen.thickness]);



% netgen.angleStd = interp1(YTSI_Ratio,XangleStd,TSI_RatioToMatch);
netgen.angleStd = VFcn(averageDensity,TSI_RatioToMatch);

% Initializing Run variables
simulation(1).type = 'Hyg';
simulation(1).dir  = '';

simulation(2).type = 'E';
simulation(2).dir = '_x';

simulation(3).type = 'E';
simulation(3).dir = '_y';

% Generate networks
PaperDir = generateNetworks(netgen,ctrl)

fibnet.length = 2000;
fibnet.width = 2000;
% fibnet runs
executeFibnetRuns(PaperDir,netgen,fibnet,simulation,ctrl,materialProps)


% Import results
[results,extraResults] = importFibnetOutputs(PaperDir,netgen,fibnet,simulation,ctrl)



% Plot the results
selE = strcmp({results.type},'E');
selH = strcmp({results.type},'Hyg');

coeffs = determineRHtoMCRatio(ctrl,'handsheetLowDens');
   



% DMA data taken from:
% C:\Users\augus\Documents\softwareProjects\appliedNetworkCurl\data\sheetData
% Unbeaten
experimentalData = [1494.914    0.167705694
                    1361.01     0.184982605
                    1684.616    0.168184264
                    900.4948    0.227689482
                    911.8092    0.213632986
                    976.9585    0.215860242 ]; % DMA data, taken from
experimentalData(:,1) = experimentalData(:,1)*1e6;
experimentalData(:,2) = experimentalData(:,2)./(diff(polyval(coeffs,[33 66])));
experimentalSTD = [ 105.131  0.00739    
                     61.122  0.02168    
                    141.214  0.00829  
                     31.380  0.00651   
                     65.639  0.00456   	
                     26.774  0.01822]; % DMA data, taken from
experimentalSTD(:,1) = experimentalSTD(:,1)*1e6;
experimentalSTD(:,2) = experimentalSTD(:,2)./(diff(polyval(coeffs,[33 66])));

linFitExp = polyfit(experimentalData(:,1).*1e-9,experimentalData(:,2),1)
linFitModel = polyfit([results(selE).data].*1e-9,[results(selH).data],1)

fitXvec = [min([min(experimentalData(:,1).*1e-9),min([results(selE).data].*1e-9)]) max([max(experimentalData(:,1).*1e-9),max([results(selE).data].*1e-9)]) ];

A=figure('color','w','units','centimeters','OuterPosition',[10 10 16 16]);

plot(fitXvec, polyval(linFitExp,fitXvec),'color',ctrl.colors(1,:))
hold on   
plot(fitXvec, polyval(linFitModel,fitXvec),'color',ctrl.colors(2,:))
% plot(experimentalData(:,1).*1e-9,experimentalData(:,2),'o','color',ctrl.colors(1,:),ctrl.lineInstructions{:})
hold on   
% plot(experimentalSTD(:,1).*1e-9,experimentalSTD(:,2),'o','color',ctrl.colors(1,:),ctrl.lineInstructions{:},'HandleVisibility','off')
errorbar(experimentalData(:,1).*1e-9,experimentalData(:,2), ...
         experimentalSTD(:,2), ...
         experimentalSTD(:,2), ...
         experimentalSTD(:,1).*1e-9, ...
         experimentalSTD(:,1).*1e-9, ...
        'o','color',ctrl.colors(1,:),'markersize',8)
hold on 

plot([results(selE).data].*1e-9,[results(selH).data],'s','color',ctrl.colors(2,:),ctrl.lineInstructions{:});


xlabel('Sheet modulus [GPa]','interpreter',ctrl.interpreter)
ylabel('Hygroexpansion $\beta$ [\%/\%mc]','interpreter',ctrl.interpreter)
xlim([0.8 1.8])
legend('Data (TH)','Model','location','northeast','interpreter',ctrl.interpreter)
set(gca,'TickLabelInterpreter',ctrl.interpreter,'fontsize',14)
print('cdl07_3','-dpng','-r0')
pause(0.5)


%%% Control data set: Change the density to higher value, see what comes
%%% out
ctrl.runDir = [ctrl.workDir filesep 'calib_highDens\'];
TSI_RatioToMatch = [1.49 1.70 2.02];

netgen.grammage = [86.10 76.0 86.2];%82.24;
netgen.thickness = [143.4 129.9 143.16];%139;
% netgen.angleStd = interp1(YTSI_Ratio,XangleStd,TSI_RatioToMatch);
averageDensity = 1e-3*[netgen.grammage]./(1e-6*[netgen.thickness]);
netgen.angleStd = VFcn(averageDensity,TSI_RatioToMatch);

ctrlHandsheetHighDens(ctrl,netgen,fibnet,simulation,materialProps)


%%% Control data set #2: Change the density to higher value
%%% This corresponds to the handsheet data
% Initializing Run variables
simulation(1).type = 'Hyg';
simulation(1).dir  = '';

simulation(2).type = 'E';
simulation(2).dir = '_x';

simulation(3).type = 'E';
simulation(3).dir = '_y';

thicknessToUse = [153 122 103];
grammageToUse = [81 81 81];
counter = 0;

isothermToImport = {'handsheetLowDens','handsheetMidDens','handsheetHighDens'};

clear res3
for aLoop = 1:3
    ctrl.runDir = [ctrl.workDir filesep 'calib_isoSheets' num2str(aLoop) filesep];
    netgen.grammage = grammageToUse(aLoop);
    netgen.thickness = thicknessToUse(aLoop);
    netgen.angleStd = 0;

    % Generate networks
    PaperDir = generateNetworks(netgen,ctrl)

    % fibnet runs
    executeFibnetRuns(PaperDir,netgen,fibnet,simulation,ctrl,materialProps)

    % Import results
    [res3Temp,~] = importFibnetOutputs(PaperDir,netgen,fibnet,simulation,ctrl)

    % Plot the results
%     selE = strcmp({results.type},'E');
%     selH = strcmp({results.type},'Hyg');
% 
    coeffs(aLoop,:) = determineRHtoMCRatio(ctrl,isothermToImport{aLoop});
    for bLoop = 1:numel(res3Temp)
        counter = counter + 1;
       res3(counter) = res3Temp(bLoop); 
    end
    dividerMX(aLoop) = diff(polyval(coeffs(aLoop,:),[33 66]));
end

selE = strcmp({res3.type},'E');
selH = strcmp({res3.type},'Hyg');
experimentalData = [1818    0.1604
                    3309    0.1957
                    4064    0.2040]; % DMA data, taken from
experimentalData(:,1) = experimentalData(:,1)*1e6;
experimentalData(:,2) = experimentalData(:,2)./dividerMX';
experimentalSTD = [126.984  0.005
                   332.274  0.018
                   653.568  0.039]; % DMA data, taken from
experimentalSTD(:,1) = experimentalSTD(:,1)*1e6;
experimentalSTD(:,2) = experimentalSTD(:,2)./dividerMX';

linFitExp = polyfit(experimentalData(:,1).*1e-9,experimentalData(:,2),1)
linFitModel = polyfit([res3(selE).data].*1e-9,[res3(selH).data],1)

fitXvec = [min([min(experimentalData(:,1).*1e-9),min([res3(selE).data].*1e-9)]) max([max(experimentalData(:,1).*1e-9),max([res3(selE).data].*1e-9)]) ];

%figure('color','w','units','centimeters','OuterPosition',[10 10 16 16]);
figure(A)
plot(fitXvec, polyval(linFitExp,fitXvec),'color',ctrl.colors(1,:))
hold on   
plot(fitXvec, polyval(linFitModel,fitXvec),'color',ctrl.colors(2,:))
% plot(experimentalData(:,1).*1e-9,experimentalData(:,2),'o','color',ctrl.colors(1,:),ctrl.lineInstructions{:})
hold on   
errorbar(experimentalData(:,1).*1e-9,experimentalData(:,2), ...
         experimentalSTD(:,2), ...
         experimentalSTD(:,2), ...
         experimentalSTD(:,1).*1e-9, ...
         experimentalSTD(:,1).*1e-9, ...
        'o','color',ctrl.colors(1,:),'markersize',8)
hold on 

plot([res3(selE).data].*1e-9,[res3(selH).data],'s','color',ctrl.colors(2,:),ctrl.lineInstructions{:});

xlabel('Sheet modulus [GPa]','interpreter',ctrl.interpreter)
ylabel('Hygroexpansion $\beta$ [\%/\%mc]','interpreter',ctrl.interpreter)
xlim([0.8 1.8])
legend('Data (TH)','Model','location','northeast','interpreter',ctrl.interpreter)
set(gca,'TickLabelInterpreter',ctrl.interpreter,'fontsize',14)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CURL SIMULATIONS
%
netgen.thickness = 180; 
netgen.grammage = 80;

TSI_RatioToMatch = [1.44 1.63 1.78];

% netgen.angleStd = interp1(YTSI_Ratio,XangleStd,TSI_RatioToMatch)
averageDensity = 1e-3*[netgen.grammage]./(1e-6*[netgen.thickness]) .*ones(1,3);
netgen.angleStd = VFcn(averageDensity,TSI_RatioToMatch);

ctrl.runDir = [ctrl.workDir filesep 'runs2\'];

clear simulation
simulation(1).type = 'FreeCU';
simulation(1).dir  = '_x';

simulation(2).type = 'FreeCU';
simulation(2).dir  = '_y';

% Generate networks
PaperDir = generateNetworks(netgen,ctrl)


% fibnet runs
executeFibnetRuns(PaperDir,netgen,fibnet,simulation,ctrl,materialProps)

% Import results
curlResults = importFibnetOutputs(PaperDir,netgen,fibnet,simulation,ctrl)

figure
plot([results(selE).data].*1e-9,[curlResults.data],'s','color',ctrl.colors(2,:),ctrl.lineInstructions{:});
pause(0.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CURL SIMULATIONS
%
ctrl.runDir = [ctrl.workDir filesep 'runs3\'];

clear simulation
simulation(1).type = 'Curl';
simulation(1).dir  = '_x';

simulation(2).type = 'Curl';
simulation(2).dir  = '_y';



% Generate networks
PaperDir = generateNetworks(netgen,ctrl)


% fibnet runs
executeFibnetRuns(PaperDir,netgen,fibnet,simulation,ctrl,materialProps)

% Import results
curlResults2 = importFibnetOutputs(PaperDir,netgen,fibnet,simulation,ctrl)

figure
plot([results(selE).data].*1e-9,[curlResults2.data],'s','color',ctrl.colors(2,:),ctrl.lineInstructions{:});




