function [lc,width,wallTkn,curl,fibrillation] = generalizedPulpDataImport_v2(targetFile,breakFlag)
if nargin < 2
    breakFlag = 0;
end

if contains(targetFile,'outOG.txt') 
    fprintf('%20s %s\n','','importOutOG.m is used to import the data.');
    networkData = importOutOG(targetFile);
    
    networkData(networkData(:,2)<eps,:) = [];
    networkData(networkData(:,3)<eps,:) = [];
    networkData(networkData(:,4)<eps,:) = [];
    networkData(networkData(:,5)<eps,:) = [];
    networkData(networkData(:,6)<eps,:) = [];
    networkData(networkData(:,16)<eps,:) = [];
    
    lc = networkData(:,3);
    width = networkData(:,4);
    wallTkn = networkData(:,5);
    curl = networkData(:,6)/100;
    
elseif contains(targetFile,'Stora32')
    fprintf('%20s %s\n','','importStora32.m is used to import the data.');
    networkData = importStora32(targetFile);
    networkData(networkData(:,3)<eps,:) = [];
    networkData(networkData(:,4)<eps,:) = [];
    networkData(networkData(:,5)<eps,:) = [];
    networkData((networkData(:,3)-networkData(:,2))<eps,:) = [];
    
    lc = networkData(:,3);
    width = networkData(:,4);
    wallTkn = networkData(:,5);
    curl = networkData(:,3)./networkData(:,2)-1;
    
elseif contains(targetFile,'Merge')
    % Number length width form area perimeter
    fprintf('%20s %s\n','','importMossab.m is used to import the data.');
    fprintf('%20s %s\n','','Warning: No wall thickness available.');
    networkData = importMossab(targetFile);
    
    networkData(networkData(:,2)<eps,:) = [];
    networkData(networkData(:,3)<eps,:) = [];
    networkData(networkData(:,4)<eps,:) = [];
    networkData(networkData(:,5)<eps,:) = [];
    networkData(networkData(:,6)<eps,:) = [];
    
    
    
    lc = networkData(:,2)/1000;
    width = networkData(:,3);
    wallTkn = nan(size(lc));
    curl = [100./networkData(:,4) - 1]';

elseif contains(targetFile,'collectedPulps') || breakFlag
    fprintf('%20s %s\n','','importOutOG.m is used to import the data.');
    networkData = importOutOG(targetFile);
    
    networkData(networkData(:,2)<eps,:) = [];
    networkData(networkData(:,3)<eps,:) = [];
    networkData(networkData(:,4)<eps,:) = [];
    networkData(networkData(:,5)<eps,:) = [];
    networkData(networkData(:,6)<eps,:) = [];
    %networkData(networkData(:,16)<eps,:) = [];
    
    lc = networkData(:,3);
    width = networkData(:,4);
    wallTkn = networkData(:,5);
    curl = networkData(:,6)/100;    
    fibrillation = zeros(size(networkData(:,16),1));
else
    fprintf('%20s %s\n','','importSofiaPulp.m is used to import the data.');
    networkData = importSofiaPulp(targetFile);
    
    % Extract the four parameters.
    lc = networkData(:,4);
    width = networkData(:,11);
    wallTkn = networkData(:,12);
    curl = networkData(:,4)./networkData(:,5)-1;
    
end