%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FIBNET CALIBRATE
% 
% fitTensileTest.m is the main function of the fibnetCalibrate repository, hosted at:
%
%   https://github.com/abrandberg/fibnetCalibrate
%
% The essential task of this script is to calibrate the fibnet micromechanical model for paper
% mechanics. Currently, it works for calibration in uniaxial tension, but it can easily be extended
% to other situations.
%
% This script works as follows:
%
%   TO DO
%
%
% created by: August Brandberg augustbr at kth dot se
% date: 2021-06-01
%

% Meta instructions
clear; close all; clc
format compact
addpath('auxilliaryFunctions')



% CTRL structure
ctrl.workDir    = cd;
ctrl.runDir     = [ctrl.workDir filesep 'runs4' filesep];
ctrl.saveDir    = [ctrl.workDir filesep 'runs1Saves' filesep];

ctrl.myPackingPointer       = ['femInputs' filesep 'MyPacking.exe'];
ctrl.modelingDataPointer    = ['femInputs' filesep 'ModelingData.txt'];
ctrl.pulpDataFile           = ['femInputs' filesep 'Euca2.txt'];
ctrl.apdlHeader             = ['femInputs' filesep 'ssCalib_Part1.txt'];
ctrl.apdlInput              = ['femInputs' filesep 'ssCalib_Part2.txt'];
ctrl.simName                = 'delme';
ctrl.customExecutable       = 'c:\Program Files\ANSYS Inc\v150\ansys\custom\user\Mossab_drying\ansys.exe';

% Generate folders that will be needed
if not(exist(ctrl.runDir,'dir'))
    mkdir(ctrl.runDir)
end
if not(exist(ctrl.saveDir,'dir'))
    mkdir(ctrl.saveDir)
end


tic
% Set network generation options
netgen.length           = 2;                                      % Length along MD, mm
netgen.width            = 2;                                      % Length along CD, mm
netgen.grammage         = 40;                                     % Surface weight per area, g/m^2
netgen.thickness        = 40;                                    % Sheet thickness, um
netgen.angleStd         = [0];                                    % Orientation spread, deg
netgen.interfaceAngle   = 5;


fibnet(1).widthIn        = 1000;                                   % FEM model length, um
fibnet(1).lengthIn       = 1000;                                   % FEM model width, um
fibnet(1).dirIn            = 0;                                      % Prerotation angle, e.g. 0, 45, 90

ctrl.optiVar = 'Stiffness';

opts = optimset('Display','iter','TolFun',0.05,'TolX',0.05);
% Solver settings



% Experimental value to match:
Eexp = [2e9 ];

assert(numel(Eexp) == numel(fibnet))

% Always specify pairs of paramName, x
% Always specify x as a multiplier of magnitude circa 1, to help solver.
paramName = {'EFiberMult'};
xIn         = [          1 ];


lossFcn = @(x) generateFibnetResult(x,ctrl,netgen,fibnet,paramName) - Eexp;

costFcn = @(x) sqrt(  ...
                    1/length(Eexp)* ... 
                    sum( ...
                            ( ...
                                lossFcn(x) ...
                                 ./Eexp ... 
                                ).^2 ...
                             ) ...
                    );


[xOut,fval] = fminsearch(costFcn,xIn,opts);


generateFibnetResult(xOut,ctrl,netgen,fibnet,paramName)
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now that we have the Fiber modulus, we can fit the KOD factor
tic
fibnet(1).widthIn        = 1000;                                   % FEM model length, um
fibnet(1).lengthIn       = 1000;                                   % FEM model width, um
fibnet(1).dirIn            = 0;                                      % Prerotation angle, e.g. 0, 45, 90
fibnet(1).EFiberMult     = xOut;


fibnet(2).widthIn        = 1000;                                   % FEM model length, um
fibnet(2).lengthIn       = 1000;                                   % FEM model width, um
fibnet(2).dirIn            = 45;                                      % Prerotation angle, e.g. 0, 45, 90
fibnet(2).EFiberMult     = xOut;

fibnet(3).widthIn        = 1000;                                   % FEM model length, um
fibnet(3).lengthIn       = 1000;                                   % FEM model width, um
fibnet(3).dirIn            = 90;                                      % Prerotation angle, e.g. 0, 45, 90
fibnet(3).EFiberMult     = xOut;

Eexp = [2e9 1.5e9 1e9];
assert(numel(Eexp) == numel(fibnet))


ctrl.optiVar = 'Stiffness';


paramName = {'KODMult'};
xIn       = [       1 ];

lossFcn = @(x) generateFibnetResult(x,ctrl,netgen,fibnet,paramName) - Eexp;

costFcn = @(x) sqrt(  ...
                    1/length(Eexp)* ... 
                    sum( ...
                            ( ...
                                lossFcn(x) ...
                                 ./Eexp ... 
                                ).^2 ...
                             ) ...
                    );


[xOut,fval] = fminsearch(costFcn,xIn,opts);


generateFibnetResult(xOut,ctrl,netgen,fibnet,paramName) 

toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main steps for next fitting parameter
%
% 1. Add another placeholder variable to the end of "ssCalib_Part1.txt"
%
% 2. Append KODMult to the fibnet structure, as we did for EFiberMult
%
% 3. Define a set of experimental results (Exp above)
%
% 4. Define a new fitting/calibration parameter in "paramName". Supply it
%    with a starting guess, and make sure the starting guess is around 1.
%    
%    That means, if you have e.g. the tangent modulus, then the starting
%    guess needs to be a multiplier. It CANNOT be the real value, because
%    FMINSEARCH will not understand the relevance of the magnitude.
%
% 5. Call the loss and the cost function so that their implicit variables
%    are updated (i.e., ctrl, fibnet, netgen).
%
% 6. Ensure that the correct postprocessing function is used. This part is not
%    changed between the two examples above, but can be done via e.g.
%    passing a function rather than using a hardcoded function at the end
%    of generateFibnetResult.m.
%
%

fibnet = fibnet(1);
fibnet(1).widthIn        = 1000;                                   % FEM model length, um
fibnet(1).lengthIn       = 1000;                                   % FEM model width, um
fibnet(1).dirIn          = 0;                                      % Prerotation angle, e.g. 0, 45, 90
fibnet(1).EFiberMult     = fibnet(1).EFiberMult;
fibnet(1).KODMult        = xOut;


ctrl.optiVar = 'Strength';
Eexp = [20e6];

fibnet(1).deboIn        = 1;
paramName = {'SBMult'};
xIn       = [      0.2 ];

lossFcn = @(x) generateFibnetResult(x,ctrl,netgen,fibnet,paramName) - Eexp;

costFcn = @(x) sqrt(  ...
                    1/length(Eexp)* ... 
                    sum( ...
                            ( ...
                                lossFcn(x) ...
                                 ./Eexp ... 
                                ).^2 ...
                             ) ...
                    );

[xOut,fval] = fminsearch(costFcn,xIn,opts);

generateFibnetResult(xOut,ctrl,netgen,fibnet,paramName) 










