function [XangleStd,YTSI_Ratio] = fiberOrientationAnisotropySurrogateModel(netgen,fibnet,ctrl,materialProps)
% Takes the geometry you want to use, generates a lot of networks, submits
% them to tensile load cases

if exist('surrogateModelFiberAnisotropy.mat','file')
    load('surrogateModelFiberAnisotropy.mat','netgenRef','fibnetRef','ctrlRef','XangleStd','YTSI_Ratio')
else
    netgenRef = nan;
    fibnetRef = nan;
    ctrlRef = nan;
end

% Check that fields match
conditionOne = isequaln(netgen,netgenRef) && isequaln(fibnet,fibnetRef) && isequaln(ctrl,ctrlRef);

if conditionOne
    
else
    % Method A: Tensile load case, TSI-ratio via ratio of elastic moduli
    simulation(1).type = 'E';
    simulation(1).dir = '_x';

    simulation(2).type = 'E';
    simulation(2).dir = '_y';

    % Generate networks
    PaperDir = generateNetworks(netgen,ctrl)


    % fibnet runs
    executeFibnetRuns(PaperDir,netgen,fibnet,simulation,ctrl,materialProps);


    % Import results
    results = importFibnetOutputs(PaperDir,netgen,fibnet,simulation,ctrl);



    % Plot the results
    selE = strcmp({results.type},'E');
    selE_MD = selE & strcmp({results.dir},'_x');
    selE_CD = selE & strcmp({results.dir},'_y');
    
    

    XangleStd = netgen.angleStd;
    YTSI_Ratio = [results(selE_MD).data]./[results(selE_CD).data];
    
    nanIdx = isnan(YTSI_Ratio);
    XangleStd(nanIdx) = [];
    YTSI_Ratio(nanIdx) = [];
    
    netgenRef = netgen;
    fibnetRef = fibnet;
    ctrlRef = ctrl;
    save('surrogateModelFiberAnisotropy.mat','netgenRef','fibnetRef','ctrlRef','XangleStd','YTSI_Ratio','-v7.3')
    
%     rmdir(ctrl.runDir)
end
% Method B: Fiber geometry postprocessing, TSI-ratio via fiber orientation



