function [nodalData,realSetData,elementData] = importNetworks(targetDir,networkName)

targetPointer = horzcat(targetDir,'\',networkName);


disp('-> Importing nodal data.')
tic
nodalData = importNodalData(horzcat(targetPointer,'.xyz'));
toc

disp('-> Importing type data')
tic 
realSetData = importTypeData(horzcat(targetPointer,'.typ'));
toc

disp('-> Importing element data')
tic 
elementData = importElementData(horzcat(targetPointer,'.nod'));
toc
