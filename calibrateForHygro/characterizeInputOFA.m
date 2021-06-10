function pulpCharacterization = characterizeInputOFA(ctrl)

addpath('auxilliaryFunctions\pulpCharacterization')

[lc,width,wallTkn,curl,~] = generalizedPulpDataImport_v2(ctrl.pulpDataFile,1);

% plotPulp(ctrl,lc,width,wallTkn,curl,fibrillation)


pulpCharacterization(1).name = 'Reference file';
pulpCharacterization(1).obs = length(lc);
pulpCharacterization(1).lc = lc;
pulpCharacterization(1).width = width;
pulpCharacterization(1).wallTkn = wallTkn;
pulpCharacterization(1).curl = curl;



rmpath('auxilliaryFunctions\pulpCharacterization')













