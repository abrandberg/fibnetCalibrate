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
% date: 2021-05-06
%

% Meta instructions
clear; close all; clc
addpath('auxilliaryFunctions')




% fixedParams is a struct that provides build information for the finite element model.
fixedParams.variable        = 'Efiber';         % {String} Variable to fit
fixedParams.debonding       = 0;                % {Boolean} Debonding active (1) or not (0)
fixedParams.plasticity      = 0;                % {Boolean} Fiber plasticity active (1) or not (0)
fixedParams.width           = 4e3;              % {Double} [um] Network width
fixedParams.length          = 4e3;              % {Double} [um] Network length
fixedParams.Ebond           = 1;                % {Double} [-] Multiplier for bond stiffness
fixedParams.Pfiber          = 1;                % {Double} [-] Multiplier for fiber plasticity
fixedParams.Sbond           = 1;                % {Double} [-] Multiplier for bond strength
fixedParams.appliedStrain   = 0.2;              % {Double} [%] Applied strain in per cent
fixedParams.numSteps        = 2;                % {Integer} [-] Number of solution steps

% fittingThresholds is a struct that provides acceptance criteria for the fits. For example, some sheets
% may have very compliant bonds (e.g. sheets made from TMP pulp where bonding is poor and fibers 
% remain quite oval/open inside the sheet). You can specify acceptance criteria which will force the 
% script to vary "secondary parameters" such as the bond stiffness when fitting the elastic modulus.
fittingThresholds.EfiberLowerBound = 0.4;       % {Double} [-] If multiplier goes below, fix it
                                                % at this value and vary bond stiffness instead.
fittingThresholds.EfiberUpperBound = 2.0;       % {Double} [-] If multiplier goes above, fix it
                                                % at this value and vary bond stiffness instead.

% Optimizer settings
% The fitting is done by recasting the calibration problem as an optimization problem and then
% minimizing the cost function of that optimization problem via the standard MATLAB optimization
% functions.
%
% There are XX types of solvers to consider that each comes with advantages and disadvantages.
%
%   FMINSEARCH  - Unconstrained derivative free method.
%
%   FMINUNC     - Unconstrained derivative (analytical or numerical) method.
%
%   FMINCON     - Constrained derivative (analytical or numerical) method.



% Step 1:   Fit the elastic modulus of the fibers and the elastic stiffness of the fiber-to-fiber 
%           bonds. 
%
%           This process is done in series: First the modulus of the fibers is fitted for a 
%           reasonable bond stiffness. If the fitted fiber modulus is very low, a second loop adjusts
%           the bond stiffness downwards. This is considerably easier than calibrating both values
%           at once, as they are essentially co-linear when considering only the small strain tensile
%           response.
opts = optimset('Display','iter','TolFun',0.01,'TolX',0.02);
% Solver settings

Efac = 1;
% Starting guess, multiplier for the fiber elastic modulus

optiFcnCaller = @(x) generateInputsFibnet(x,fixedParams);
% Cost function

[x,fval] = fminsearch(optiFcnCaller,Efac,opts);
% Solver call. x returns the "optimal" multiplier factor for the fiber elastic modulus

fprintf('%s\n','Fitting of fiber modulus complete:');
fprintf('%20s %20s %20s %20s\n','','x0','xEnd','Error for xEnd [%]');
fprintf('%20s %20.4f %20.4f %20.4f\n','',Efac,x,fval);



if ( x < fittingThresholds.EfiberLowerBound ) 
    fprintf('%s\n','Modulus outside accepted bounds.');
    fprintf('%s %2.2f\n','Modulus multiplier will be fixed at ',fittingThresholds.EfiberLowerBound);
    fprintf('%s\n','Bond stiffness will be calibrated now.');
    fixedParams.Efiber = fittingThresholds.EfiberLowerBound;

elseif ( x > fittingThresholds.EfiberUpperBound )
    fprintf('%s\n','Modulus outside accepted bounds.');
    fprintf('%s %2.2f\n','Modulus multiplier will be fixed at ',fittingThresholds.EfiberUpperBound);
    fprintf('%s\n','Bond stiffness will be calibrated now.');
    fixedParams.Efiber = fittingThresholds.EfiberUpperBound;

end

if ( x < fittingThresholds.EfiberLowerBound ) || ( x > fittingThresholds.EfiberUpperBound )
    fixedParams = rmfield(fixedParams,'Ebond');
    fixedParams.variable = 'Ebond';

    optiFcnCaller = @(x) generateInputsFibnet(x,fixedParams);

    [x,fval] = fminsearch(optiFcnCaller,Bfac,opts);
    fixedParams.Ebond = x;

else
    fprintf('%s\n','Modulus accepted.');
    fixedParams.Efiber = x;
    fixedParams.Ebond = 1;

end

fprintf('%s\n','Fitting of bond modulus complete:');
fprintf('%20s %20s %20s %20s\n','','x0','xEnd','Error [%]');
fprintf('%20s %20.4f %20.4f %20.4f\n','',Bfac,fixedParams.Ebond,fval);





% Step 2:   Fit the inelastic constants: Bond damage and fiber plasticity
%
%           This process can be done in some different ways: The plastic properties can be simulated
%           first, followed by the damage parameters (sequential estimation) or all three parameters
%           can be investigated simultaneously. 
%
%           Sometimes it is necessary to make a judgement call about which of these is best in a 
%           given application.

%solverSetting = 'sequential';
solverSetting = 'parallel';


if strcmp(solverSetting,'sequential')

    % New third step is to fit the plastic properties
    fixedParams = rmfield(fixedParams,'Pfiber');
    fixedParams.debonding = 0;
    fixedParams.plasticity = 1;
    fixedParams.appliedStrain = 1;
    fixedParams.numSteps = 100;
    fixedParams.variable = 'Pfiber';
    Pguess = 180;
    Pguess2 = 6e3;
    Pfac = 1;
    optiFcnCaller = @(x) generateInputsFibnet(x,fixedParams);
    [x,fval] = fminsearch(optiFcnCaller,Pfac,opts);

    fixedParams.Pfiber = x;

    fprintf('%s\n','Fitting of fiber plasticity complete:');
    fprintf('%20s %20s %20s %20s %20s %20s\n','','Starting Sy [MPa]','Final Sy [MPa]','Starting Et [MPa]','Final Et [MPa]','Error [%]');
    fprintf('%20s %20.4f %20.4f %20.4f %20.4f %20.4f\n','',Pguess*Pfac,Pguess*x,Pguess2*Pfac,Pguess2*x,fval);



    % Third step is to fit the failure point, starting with strength and
    % working on the bond strength
    %fprintf('%s\n'
    fixedParams = rmfield(fixedParams,'Sbond');
    fixedParams.debonding = 1;
    fixedParams.plasticity = 0;
    fixedParams.appliedStrain = 2;
    fixedParams.numSteps = 100;
    fixedParams.variable = 'Sbond';
    Sguess = 10000;
    Sfac = 0.6;
    optiFcnCaller = @(x) generateInputsFibnet(x,fixedParams);
    [x,fval] = fminsearch(optiFcnCaller,Sfac,opts);

    fixedParams.Sbond = x;

    fprintf('%s\n','Fitting of bond strength complete:');
    fprintf('%20s %20s %20s %20s\n','','Starting S [uN]','Final S [uN]','Error [%]');
    fprintf('%20s %20.4f %20.4f %20.4f\n','',Sguess*Sfac,Sguess*x,fval);


elseif strcmp(solverSetting,'parallel')
    
    options = optimoptions('fsolve','Display','iter','FiniteDifferenceStepSize',0.005);
    options = optimoptions('fsolve','Display','iter','FiniteDifferenceStepSize',0.02);
    
    
    fixedParams.variable = 'Nonlin';
    fixedParams.debonding = 1;
    fixedParams.plasticity = 1;
    fixedParams.appliedStrain = 1.5;
    fixedParams.numSteps = 100;
%     fixedParams.Ebond = 1;
    Nguess = [180 6e3 10000]; % PfiberS PfiberH SBond
    Nfac = [0.5 2.0 0.2 ];
    Nfac = [3 0.4 0.3 ];
%     Nfac = [1 1 1];
    Nfac = [2.1778               0.4193               0.4994]
    Nfac = [0.51778               1.4193               0.4994]
    Nfac = [2.0 1.4 0.8]
    optiFcnCaller = @(x) generateInputsFibnet(x,fixedParams);
%     [x,fval] = fminsearch(optiFcnCaller,Nfac,opts);
%     
%     
%     optiFcnCaller = @(x) generateInputsFibnet(x,fixedParams);
%     [x,fval] = fsolve(optiFcnCaller,Nfac,options);%,opts);
    
%     options2 = optimoptions('fgoalattain','Display','iter','FiniteDifferenceStepSize',0.005);
%     [x,~,attainfactor] = ...
%         fgoalattain(optiFcnCaller,Nfac,[0 0 0],[1 1 1],[],[],[],[],[0.1 0.1 0.1],[],[],options2);
    Nfac = [1.6433               1.0027               1];%0.7958] % SAVE SAVE SAVE works for 4x2
    Nfac = [1.6419 1.0127 0.9830] % Save best values for 4x2
    Nfac = [2.0 1.15 0.875] % worked for 266 um
%     Nfac = [1.6 2.40 1] 
    Nfac = [0.5 0.6 0.5] % worked for 277 um, 4 deg
    Nfac = [2.0 1.0 0.75]
    Nfac = [1.6 1.0 0.28]
    Nfac = [0.5 1.0 0.40]
    Nfac = [0.4715               0.8897               0.3707]
    Nfac = [0.3515               0.8897               0.2007]
    opt3 = optimoptions('lsqnonlin','Display','iter','FiniteDifferenceStepSize',0.02);
    [x,fval] = lsqnonlin(optiFcnCaller,Nfac,[0 0 0],[],opt3);%,opts);
    % Nfac = [0.3634               0.9917               0.2133] yields good
    % results

    % Fit everything
%     fixedParams.variable = 'Nonlin2';
%     fixedParams.debonding = 1;
%     fixedParams.plasticity = 1;
%     fixedParams.appliedStrain = 4;
%     fixedParams.numSteps = 100;
% %     fixedParams.Ebond = 1;
%     Nguess = [180 6e3 10000 10032*0.8]; % PfiberS PfiberH SBond
% %     Nfac = [0.5 2.0 0.2 ];
%     Nfac = [3 0.4 0.3 1.5];
% %     Nfac = [1 1 1];
% %     Nfac = [2.1778               0.4193               0.4994]
%     optiFcnCaller = @(x) generateInputsFibnet(x,fixedParams);
%     [x,fval] = fminsearch(optiFcnCaller,Nfac,opts);



    % Map outputs
    fixedParams.PfiberS = x(1);
    fixedParams.PfiberH = x(2);
    fixedParams.Sbond = x(3);
    
end



% Final demo
fixedParams.debonding = 1;
fixedParams.plasticity = 1;
fixedParams.appliedStrain = 2;
fixedParams.numSteps = 100;
fixedParams.variable = 'Demo';
generateInputsFibnet(0,fixedParams);


toc

