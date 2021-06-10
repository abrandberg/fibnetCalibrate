function networkStructure =  characterizeFEMInputs(ctrl,fibnet)

% Simulations to investigate are in [ctrl.runDir PaperDir{i}]


% Add path
addpath(['auxilliaryFunctions' filesep 'femInputCharacterization'])


spaceToSearch = [-0.5*fibnet.length 0.5*fibnet.length ; -0.5*fibnet.width 0.5*fibnet.width];
spaceToSearch = inf.*[-1 1 ; -1 1];

availableFolders = subdirImport(ctrl.runDir);

for aLoop = 1:numel(availableFolders) % Go through each folder, importing results and then tabulating them somehow
    
    if exist([ctrl.runDir filesep availableFolders{aLoop} filesep 'networkStructure.mat'],'file')
        load([ctrl.runDir filesep availableFolders{aLoop} filesep 'networkStructure.mat'],'thisNetwork')
    else
    
        filesInDir = subdirImport([ctrl.runDir filesep availableFolders{aLoop}],'regex','file');
        [nodalData,realSetData,elementData] = importNetworks([ctrl.runDir filesep availableFolders{aLoop}],filesInDir{1});
        tic
        thisNetwork = tabulateNetworksAB_v2(nodalData,elementData,realSetData, ...
                                                        spaceToSearch, ...
                                                        filesInDir{1},[ctrl.runDir filesep availableFolders{aLoop}]);
        toc
        save([ctrl.runDir filesep availableFolders{aLoop} filesep 'networkStructure.mat'],'thisNetwork','-v7.3')
    end
    
    
    
    
    
    networkStructure(aLoop) = thisNetwork;
    
end





% Remove path
rmpath(['auxilliaryFunctions' filesep 'femInputCharacterization'])