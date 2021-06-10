function freeSpanTemp = calculateElementSegment(elementCoordinates,contactCoordinates)

% Form an element vector
elementVector = elementCoordinates(2,:)-elementCoordinates(1,:);
elementVectorNormalized = elementVector/norm(elementVector);

% Form the vector to each contact point
d1 = contactCoordinates - elementCoordinates(1,:);

% Projecting d1 onto elementVector
localCoordinate = dot(d1,elementVectorNormalized.*ones(size(d1)),2)/norm(elementVector);

if length(localCoordinate) == 2
    freeSpanTemp = abs(diff(localCoordinate))*norm(elementVector);
elseif length(localCoordinate) == 1
    freeSpanTemp = (1-localCoordinate)*norm(elementVector);
end

if 0% plotFlag == 1
    figure();
    plot3([elementCoordinates(1,1) elementCoordinates(2,1)],[elementCoordinates(1,2) elementCoordinates(2,2)],[elementCoordinates(1,3) elementCoordinates(2,3)],'-o')
    axis equal
    hold on
    plot3(contactCoordinates(:,1),contactCoordinates(:,2),contactCoordinates(:,3),'xr')
    xlabel x; ylabel y; zlabel z;
%     quiver3(N1(1)*ones(size(d1,1),1),N1(2)*ones(size(d1,1),1),N1(3)*ones(size(d1,1),1), ...
%             d1(:,1),d1(:,2),d1(:,3))
    plot3(elementCoordinates(1,1)+localCoordinate.*elementVector(1), ...
          elementCoordinates(1,2)+localCoordinate.*elementVector(2), ...
          elementCoordinates(1,3)+localCoordinate.*elementVector(3),'ms')    
end

% a = elementCoordinates(2,:) - elementCoordinates(1,:);
% b = contactCoordinates - elementCoordinates(1,:);
% alpha = dot(a,b)/norm(a)/norm(a);
% 
% % Decide which part should be added and which part should be
% % discarded. This is done by comparing to the next element.
% 
% distToCorner1 = norm(elementCoordinates(1,:)-contactCoordinates);
% distToCorner2 = norm(elementCoordinates(2,:)-contactCoordinates);
% 
% if distToCorner1 < distToCorner2 % If distance from contact to corner node 1 is less than to corner node 2
%     freeSpanTemp = norm(a)*alpha;
% else
%     freeSpanTemp = norm(a)*(1-alpha);
% end