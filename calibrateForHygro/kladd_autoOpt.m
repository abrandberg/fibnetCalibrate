clear; close all; clc


% Set up ctrl structure
ctrl.workDir = cd;
ctrl.runDir = [ctrl.workDir filesep 'autoOpt\'];
ctrl.saveDir = [ctrl.workDir filesep 'autoOptSave\'];


ctrl.interpreter = 'latex';
ctrl.histogramInstructions = {'normalization','probability'};

ctrl.colors = [255 94 105 ;
               0 188 213]./255;
ctrl.lineInstructions = {'markersize',8,'linewidth',1.25};

ctrl.parallelNetGetRuns = 12;
ctrl.parallelFibNetRuns = 6;
ctrl.fibnetNP = 2;

ctrl.executable = 'C:\Program Files\ANSYS Inc\v150\ansys\custom\user\Hygro3D_Kurosh\ansys.exe';
ctrl.pulpDataFile = 'Euca2.txt'

ctrl.forceRewrite = 1;



weightFactor = [ 0.5 0.5 0.5 0.5 0.5 0.5]'; % 1 1 1
% ExperimentalDataToMatch = [... 1818             % Handsheet
%                            ...3309
%                            ...4064
%                            1494.914         % Low density oriented sheet
%                            1361.01 
%                            1684.616
%                            900.4948
%                            911.8092
%                            976.9585
%                            ...4706.746        % High density oriented sheet
%                            ...4057.886
%                            ...4559.684
%                            ...2703.32
%                            ...2764.12
%                            ].*1e6;...2609.578].*1e6; 
  ExperimentalDataToMatch = [                           1494.914         % Low density oriented sheet
                           1361.01 
                           1684.616
                           900.4948
                           911.8092
                           976.9585].*1e6;...2609.578].*1e6;                          
                       
                       
% experimentalSTD = [ ...126.984
%                    ...332.274 
%                    ...653.568 
%                     105.131 
%                      61.122 
%                     141.214 
%                      31.380 
%                      65.639 
%                      26.774 ]; % DMA data, taken from
experimentalSTD = [ 
                    105.131 
                     61.122 
                    141.214 
                     31.380 
                     65.639 
                     26.774 ]; % DMA data, taken from
experimentalSTD(:,1) = experimentalSTD(:,1)*1e6;


                        
materialProps.EX                        = [1];
materialProps.Kct                       = [10e3];
materialProps.interfaceAngle            = [9];
materialProps.searchRadiusMultiplier    = [1.0];


x0 = [materialProps.EX ]; % materialProps.searchRadiusMultiplier 
ctrl.fitParam = 1;
unloadFitFun = @(x,ctrl) fitSheetStiffness(x,ctrl) - ExperimentalDataToMatch ;
            
unloadFitMinFun = @(x) sqrt(  ...
                            1/length(ExperimentalDataToMatch)* ... 
                            sum( ...
                                    ( ...
                                        weightFactor .* unloadFitFun(x,ctrl) ...
                                         ./ExperimentalDataToMatch... 
                                        ).^2 ...
                                     ) ...
                            );                    
                       

% unloadFitMinFun(x)




options = optimset('Display','iter','TolFun',0.01,'TolX',0.01);
tic
[xOut,fval] = fminsearch(unloadFitMinFun,x0,options)
toc
tic
ap = fitSheetStiffness(xOut,ctrl);
toc


figure;
% plot(ExperimentalDataToMatch,'o-')
errorbar(ExperimentalDataToMatch(:,1), ...
         experimentalSTD(:,1), ...
        'o','color',ctrl.colors(1,:),'markersize',8)
hold on
plot(ap,'s-b')

% 0.95
% 1.18 bara  low dens


% Kör bond stiffness mult * Kct * dens/dens_ref
weightFactor = [ 1 1 1]';
ExperimentalDataToMatch = [1818             % Handsheet
                           3309
                           4064].*1e6;

                       
% weightFactor = [1 1 1 1 1 1]';                     
% ExperimentalDataToMatch = [4706.746        % High density oriented sheet
%                            4057.886
%                            4559.684
%                            2703.32
%                            2764.12
%                            2609.578].*1e6;
%                        
ctrl.fitParam = 2;
ctrl.EFiberMult = xOut;
x0 = 1.0;

unloadFitFun = @(x,ctrl) fitSheetStiffness(x,ctrl) - ExperimentalDataToMatch ;
            
unloadFitMinFun = @(x) sqrt(  ...
                            1/length(ExperimentalDataToMatch)* ... 
                            sum( ...
                                    ( ...
                                        weightFactor .* unloadFitFun(x,ctrl) ...
                                         ./ExperimentalDataToMatch... 
                                        ).^2 ...
                                     ) ...
                            );      
                        
                        
options = optimset('Display','iter','TolFun',0.01,'TolX',0.01);
tic
[xOut,fval] = fminsearch(unloadFitMinFun,x0,options)
toc
tic
ap = fitSheetStiffness(xOut,ctrl);
toc
%     2.1563

figure;
plot(ExperimentalDataToMatch,'o-')
hold on
plot(ap,'s-b')


ctrl.fitParam = 3;
ctrl.EBondMult = xOut;
ExperimentalDataToMatch = [ 1818             % Handsheet
                           3309
                           4064
                           1494.914         % Low density oriented sheet
                           1361.01 
                           1684.616
                           900.4948
                           911.8092
                           976.9585
                           4706.746        % High density oriented sheet
                           4057.886
                           4559.684
                           2703.32
                           2764.12
                           2609.578].*1e6;


cp = fitSheetStiffness(1,ctrl) 


figure;
plot(ExperimentalDataToMatch,'o-')
hold on
plot(cp,'s-b')




% fiberMult = 1.1938
% bondMult = 1.2250


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Time for the hygroexpansion coefficients
%
%
coeffs = determineRHtoMCRatio(ctrl,'handsheetLowDens');
weightFactor = [ 1 1 1 1 1 1]';
ExperimentalDataToMatch = [ 0.167705694
                            0.227689482
                            0.184982605
                            0.213632986
                            0.168184264
                            0.215860242 ]./(diff(polyval(coeffs,[33 66])));

                        
ctrl.fitParam = 4;
% ctrl.EFiberMult = xOut;
% ctrl.EBondMult = xOut;
x0 = 1.0;                        
                        
                        
                        
unloadFitFun = @(x,ctrl) fitSheetStiffness(x,ctrl) - ExperimentalDataToMatch ;
            
unloadFitMinFun = @(x) sqrt(  ...
                            1/length(ExperimentalDataToMatch)* ... 
                            sum( ...
                                    ( ...
                                        weightFactor .* unloadFitFun(x,ctrl) ...
                                         ./ExperimentalDataToMatch... 
                                        ).^2 ...
                                     ) ...
                            );      


                        

ctrl.fibnetNP = 4;
%xz = fitSheetStiffness(1,ctrl)
x0 = 1;
options = optimset('Display','iter','TolFun',0.02,'TolX',0.02);
tic
[xOut,fval] = fminsearch(unloadFitMinFun,x0,options)
toc
tic
ap4 = fitSheetStiffness(xOut,ctrl);
toc

figure;
plot(ExperimentalDataToMatch,'o-')
hold on
plot(ap4,'s-b')

% xOut = 1.3750
