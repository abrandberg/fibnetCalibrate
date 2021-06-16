function [result] = fitSheetStiffness(x,ctrl)

% The results to compare against are:
% Isotropic sheets, three different densities


% First step is to make an array of the different combinations that should
% be tried
% netgenAll.grammage         = [  83   82   81   76.66     72.51     84.11    76.66     72.51     84.11 ];%  86.10      76.0      86.2    86.10    76.0    86.2]';
% netgenAll.thickness        = [ 153  122  103  179.33    169.83    194.46   179.33    169.83    194.46 ];%  143.4     129.9    143.16    143.4   129.9  143.16]';
% netgenAll.width            = [   8    8    8       8         8         8        8         8         8 ];%      8         8         8        8       8       8]';
% netgenAll.length           = [   8    8    8       8         8         8        8         8         8 ];%      8         8         8        8       8       8]';
% netgenAll.angleStd         = [  79   80   81 52.0296   47.7956   45.2986  52.0296   47.7956   45.2986 ];%50.3038   46.0547   41.2320  50.3038 46.0547 41.2320]';
% simDir                     = {'_x','_x','_x',   '_x',     '_x',     '_x',    '_y',     '_y',      '_y'};%   '_x',     '_x',     '_x',    '_y',   '_y',   '_y'};

if ctrl.fitParam == 1 
netgenAll.grammage         = [   76.66     72.51     84.11    76.66     72.51     84.11 ];%  86.10      76.0      86.2    86.10    76.0    86.2]';
netgenAll.thickness        = [  179.33    169.83    194.46   179.33    169.83    194.46 ];%  143.4     129.9    143.16    143.4   129.9  143.16]';
netgenAll.width            = [       8         8         8        8         8         8 ];%      8         8         8        8       8       8]';
netgenAll.length           = [      8         8         8        8         8         8 ];%      8         8         8        8       8       8]';
netgenAll.angleStd         = [52.8346 48.6838 46.2105  52.8346 48.6838 46.2105 ];%50.3038   46.0547   41.2320  50.3038 46.0547 41.2320]';
simDir                     = {   '_x',     '_x',     '_x',    '_y',     '_y',      '_y'};%   '_x',     '_x',     '_x',    '_y',   '_y',   '_y'};
%
elseif ctrl.fitParam == 2
    netgenAll.grammage         = [  83   82   81];%[   86.10      76.0      86.2    86.10    76.0    86.2]';%[  83   82   81];
    netgenAll.thickness        = [ 153  122  103];%[   143.4     129.9    143.16    143.4   129.9  143.16]';%[ 153  122  103];
    netgenAll.width            = [   8    8    8];%[       8         8         8        8       8       8]';%[   8    8    8];
    netgenAll.length           = [   8    8    8];%[       8         8         8        8       8       8]';%[   8    8    8];
    netgenAll.angleStd         = [  79   80   81];%[51.2450   47.0398   42.3813  51.2450 47.0398 42.3813]';%[  79   80   81];
    simDir                     = {'_x','_x','_x'};%{   '_x',     '_x',     '_x',    '_y',   '_y',   '_y'};%{'_x','_x','_x'};
    
    
elseif ctrl.fitParam == 3
    netgenAll.grammage         = [  83   82   81   76.66     72.51     84.11    76.66     72.51     84.11   86.10      76.0      86.2    86.10    76.0    86.2]';
    netgenAll.thickness        = [ 153  122  103  179.33    169.83    194.46   179.33    169.83    194.46   143.4     129.9    143.16    143.4   129.9  143.16]';
    netgenAll.width            = [   8    8    8       8         8         8        8         8         8       8         8         8        8       8       8]';
    netgenAll.length           = [   8    8    8       8         8         8        8         8         8       8         8         8        8       8       8]';
    netgenAll.angleStd         = [  79   80   81 52.8346   48.6838   46.2105   52.8346  48.6838   46.2105 51.2450   47.0398   42.3813  51.2450 47.0398 42.3813]';
    simDir                     = {'_x','_x','_x',   '_x',     '_x',     '_x',    '_y',     '_y',      '_y'   '_x',     '_x',     '_x',    '_y',   '_y',   '_y'};
    
elseif ctrl.fitParam == 4
    netgenAll.grammage     = [   76.66     72.51     84.11 ];%  86.10      76.0      86.2    86.10    76.0    86.2]';
netgenAll.thickness        = [  179.33    169.83    194.46 ];%  143.4     129.9    143.16    143.4   129.9  143.16]';
netgenAll.width            = [       8         8         8 ];%      8         8         8        8       8       8]';
netgenAll.length           = [      8         8          8 ];%      8         8         8        8       8       8]';
netgenAll.angleStd         = [52.8346 48.6838 46.2105  ];%50.3038   46.0547   41.2320  50.3038 46.0547 41.2320]';
simDir                     = {   '',     '',     ''};%   '_x',     '_x',     '_x',    '_y',   '_y',   '_y'};
%
end
%
fibnet.width            = 4000;%[3500 ];
fibnet.length           = 4000;%[3500 ];

% simulation(1).type = 'E';
% simulation(1).dir = '_x';
simType = 'E';

if ctrl.fitParam == 4
   simType = 'Hyg'; 
%    simDir = {'','','','','',''};
end


% Then we should make an informal list of all the things that we will allow
% to change.
% materialProps.EX                        = [30e3];
% materialProps.Kct                       = [10e3];
% materialProps.interfaceAngle            = [9];
% materialProps.searchRadiusMultiplier    = [1.0];
materialProps.EZ = 17e3;

counter = 0;
for aLoop = 1:numel(netgenAll.grammage)
    
    netgen.grammage = netgenAll.grammage(aLoop);
    netgen.thickness = netgenAll.thickness(aLoop);
    netgen.width = netgenAll.width(aLoop);
    netgen.length = netgenAll.length(aLoop);
    netgen.angleStd = netgenAll.angleStd(aLoop);
    
    
    simulation(1).type = simType;
    simulation(1).dir = simDir{aLoop};
    simulation(1).idx = aLoop;
    
    
    if ctrl.fitParam == 1
        materialProps.EX = x(1)*25e3;
        materialProps.Kct = 8000*((netgen.grammage/netgen.thickness)/(0.4275) );%x(2);
        netgen.interfaceAngle = 9.0;%x(3);
        materialProps.searchRadiusMultiplier = 1.00;%x(2);
    elseif ctrl.fitParam == 2
        materialProps.EX = ctrl.EFiberMult*25e3;
%         materialProps.Kct = x(1)*8000 * (netgen.grammage/netgen.thickness)/(0.4275);%x(2);
        materialProps.Kct = 8000 * ((netgen.grammage/netgen.thickness)/(0.4275) ) .*x(1);%x(2);

        netgen.interfaceAngle = 9.0;%x(3);
        materialProps.searchRadiusMultiplier = 1.00;%x(2);
        
    elseif ctrl.fitParam == 3
        materialProps.EX = ctrl.EFiberMult*25e3;
        materialProps.Kct = 8000 * ((netgen.grammage/netgen.thickness)/(0.4275) ).*ctrl.EBondMult;
        netgen.interfaceAngle = 9.0;%x(3);
        materialProps.searchRadiusMultiplier = 1.00;%x(2);
        
    elseif ctrl.fitParam == 4
        materialProps.EX = ctrl.EFiberMult*25e3;
        materialProps.Kct = 8000 * ((netgen.grammage/netgen.thickness)/(0.4275) ).*ctrl.EBondMult;
        netgen.interfaceAngle = 9.0;%x(3);
        materialProps.searchRadiusMultiplier = 1.00;%x(2);
        
        materialProps.EZ = x(1)*17e3; 
    end
    
    
    % Now we should build networks and solve them
    PaperDir = generateNetworks(netgen,ctrl);

    % fibnet runs
    executeFibnetRuns(PaperDir,netgen,fibnet,simulation,ctrl,materialProps);


end

for bLoop = 1:numel(netgenAll.grammage)
    
    simulation(1).type = simType;
    simulation(1).dir = simDir{bLoop};
    simulation(1).idx = bLoop;
    
    netgen.grammage = netgenAll.grammage(bLoop);
    netgen.thickness = netgenAll.thickness(bLoop);
    netgen.width = netgenAll.width(bLoop);
    netgen.length = netgenAll.length(bLoop);
    netgen.angleStd = netgenAll.angleStd(bLoop);
    % Then we should import a specific piece of data (i.e., the force so we can
    % calculate the modulus).
    % Import results
    % Then we calculate the modulus
    [res3Temp,~] = importFibnetOutputs(PaperDir,netgen,fibnet,simulation,ctrl);

    for vLoop = 1:numel(res3Temp)
        counter = counter + 1;
        res3(counter) = res3Temp(vLoop);
    end

% When we have the modulus, we compare it against the modulus of the
% corresponding experimental observation.
%
% This comparison is essentially the loss.
%


% Once we have the loss, we can return it. We are done!!


end


result = [res3.data]';



