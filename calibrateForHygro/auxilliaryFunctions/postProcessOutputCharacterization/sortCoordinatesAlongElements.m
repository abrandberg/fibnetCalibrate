function sortedListOfContacts = sortCoordinatesAlongElements(nodalData,elementData,cordTemp1,plotFlag)
% 
%
% The problem is as follows:
% - We need to calculate the free spans between bonds in the network.
% - To do that we need to create an ordered list, where the bonds are
%   sequential. Otherwise, we will calculate the wrong free span.
% - To create an ordered list between elements is simple, because we know
%   that the beam elements themselves are along a line.
% - However, within that line we do not know which order they have.
%
% The solution is as follows:
% 1. For each element, we isolate the contacts on that element. Each
%    contact is labeled (C_i).
%
% 2. We form the vector from node 1 of the element to node 2 (end to end).
%    Call this vector E1.
%
% 3. We project the vector from N1 to C1 onto the vector E1.
%
%                 x C1
%                /|
%               / |
%              /  |
%             /   |
%            /    |
%        N1 o--------------------------------------o N2
%
%
%

% For each beam element that has at least one contact element attached to
% it.

uniqueElements = unique(cordTemp1(:,9)); % Find all the unique elements in the model.

for xLoop = 1:length(uniqueElements)
    evalElement = uniqueElements(xLoop);
    
    % Form the vector
    N1 = nodalData(elementData(evalElement,2),2:4);
    N2 = nodalData(elementData(evalElement,3),2:4);
    
    elementVector = N2-N1;
    elementVectorNormalized = elementVector/norm(N2-N1);
    
    % Form the vector to each contact point
    d1 = cordTemp1(cordTemp1(:,9)==evalElement,4:6) - N1;
    
    % Projecting d1 onto elementVector
    localCoordinate = dot(d1,elementVectorNormalized.*ones(size(d1)),2)/norm(N2-N1);
    cordTemp1(cordTemp1(:,9)==evalElement,13) = localCoordinate;
    
    
    if plotFlag == 1
        figure();
        plot3([N1(1) N2(1)],[N1(2) N2(2)],[N1(3) N2(3)],'-o')
        axis equal
        hold on
        plot3(cordTemp1(cordTemp1(:,9)==evalElement,4),cordTemp1(cordTemp1(:,9)==evalElement,5),cordTemp1(cordTemp1(:,9)==evalElement,6),'xr')
        xlabel x; ylabel y; zlabel z;
        quiver3(N1(1)*ones(size(d1,1),1),N1(2)*ones(size(d1,1),1),N1(3)*ones(size(d1,1),1), ...
                d1(:,1),d1(:,2),d1(:,3))
        plot3(N1(1)+localCoordinate.*elementVector(1),N1(2)+localCoordinate.*elementVector(2),N1(3)+localCoordinate.*elementVector(3),'ms')    
    end
end
 
% Sort the result, first on element number, then on position inside that
% element.
[~,sortIdx1] = sort(cordTemp1(:,13));
[~,sortIdx2] = sort(cordTemp1(sortIdx1,9));
sortIdx = sortIdx1(sortIdx2);
sortedListOfContacts = cordTemp1(sortIdx,:);




