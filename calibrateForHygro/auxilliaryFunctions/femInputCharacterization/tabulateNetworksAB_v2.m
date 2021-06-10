function networkStructure = tabulateNetworksAB_v2(nodalData,elementData,realSetData,probeSpace,networkName,networkDir)
%function
%tabulateNetworks(nodalData,elementData,realSetData,probeSpace,networkName,networkDir)
%uses the fibnet data to calculate various things. In this version, the
%characterization of interest is the fiber orientation.
%
% INPUT:    nodalData   : All nodal coordinates {IDX X Y Z}
%           elementData : All element connectivity {IDX START_NODE END_NODE MID_NODE REAL_IDX MAT_IDX}
%           realSetData : Real values {IDX TYPE WIDTH HEIGHT WALL_TKN}
%                         TYPE = 1 : Solid rectangle
%                         TYPE = 2 : Hollow rectangle
%           probeSpace  : Rectangle specifying region of interest ([xMin xMax ; yMin yMax])
%           networkDir  : The directory to import from
%           networkName : A string matching the network to import
%
% OUTPUT:   networkStructure : Calculated data, in this case fiber segment
%                              orientation. Also contains some fields for
%                              filtering purposes.
%
% ABOUT:
%
% TO DO: 
% - This script does not bother itself with partial fiber segments (where 1
%   part is inside and 1 part is outside the region of interest). I have
%   not checked but I do not think it matters. /AB, 09-07-2019
%
% created: 09-07-2019
% author: August Brandberg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xDim = probeSpace(1,:); % Boundaries for the region of interest
yDim = probeSpace(2,:);

% Pre-allocation
newFiberMarker = diff(elementData(:,5));
newFiberMarker = [1 ; newFiberMarker];
newFiberIndex = find(newFiberMarker);
angYX = nan(size(elementData(:,1)));%[];
% angXZ = nan(size(elementData(:,1)));%[];
cMap = parula(128);
% Cartesian axes
yAx = [0 1 0]';
% zAx = [0 0 1]';
% 
angYXTrunc = [];
fiberPropTrunc = [];

counter = 1;
for xLoop = 1:sum(newFiberMarker)-1   % Goes through the data fiber by fiber (useful for measuring fiber length)
   startPos = newFiberIndex(xLoop);
   endPos = newFiberIndex(xLoop+1)-1;
   truncMarker = startPos;
   truncSelector = 2;
   truncSelStart = 3;
   yLoop = -1;
   
   LpTemp = sqrt(sum([ nodalData(elementData(endPos,3),2) - nodalData(elementData(startPos,2),2) , ...
                       nodalData(elementData(endPos,3),3) - nodalData(elementData(startPos,2),3) , ...
                       nodalData(elementData(endPos,3),4) - nodalData(elementData(startPos,2),4) ].^2));
   
   LcTemp = sum(sqrt(sum([ nodalData(elementData(startPos:endPos,3),2) - nodalData(elementData(startPos:endPos,2),2) , ...
                           nodalData(elementData(startPos:endPos,3),3) - nodalData(elementData(startPos:endPos,2),3) , ...
                           nodalData(elementData(startPos:endPos,3),4) - nodalData(elementData(startPos:endPos,2),4) ].^2,2)));
%        LcTemp = sum(sqrt(sum([ nodalData(elementData(startPos:endPos,3),2) - nodalData(elementData(startPos:endPos,2),2)].^2,2)));
%                        
%                        
%                        
%     LcTemp = sum(sqrt(sum([ nodalData(elementData(endPos,3),2) - nodalData(elementData(endPos-1,2),2) , ...
%                            nodalData(elementData(endPos,3),3) - nodalData(elementData(endPos-1,2),3) , ...
%                            nodalData(elementData(endPos,3),4) - nodalData(elementData(endPos-1,2),4) ].^2,2)));             
                 
                 
                 
    fiberProp(xLoop,:) = [LcTemp LpTemp realSetData(elementData(startPos,5),3:5)];
    
    lc(xLoop,1) = LcTemp;
    lp(xLoop,1) = LpTemp;
    width(xLoop,1) = realSetData(elementData(startPos,5),3);
    
    if realSetData(elementData(startPos,5),5) == 0 
       wallTkn(xLoop,1) = realSetData(elementData(startPos,5),4)*0.5;
    else
       wallTkn(xLoop,1) = realSetData(elementData(startPos,5),5) ;
    end
    
    curl(xLoop,1) = LcTemp/LpTemp - 1; 
    

        
                 
%               LX = 0;   
%          for hLoop = startPos:endPos-1
%             LX = LX + sqrt(sum([ nodalData(elementData(hLoop+1,3),2) - nodalData(elementData(hLoop,2),2) ;
%                      nodalData(elementData(hLoop+1,3),3) - nodalData(elementData(hLoop,2),3) ;
%                      nodalData(elementData(hLoop+1,3),4) - nodalData(elementData(hLoop,2),4) ].^2));
%          end
                 
                 
                 %sqrt(sum(sum(diff(realCenterline,                   1,1))).^2)
   while yLoop <=(endPos-startPos)-1 % Goes through each fiber, element by element
       
       yLoop = yLoop + 1; 
       ia = startPos+yLoop;

       % This conditional checks that the start and end node of the element
       % in question is inside the region of interest so that only that
       % material is added. Since many fibers are partially inside, a while
       % loop is used to incrementally "step through" the fiber elements.
       %
       % This implementation is by no means the fastest, but it is easy
       % enough to understand and debug.
       conditionalOne = nodalData(elementData(ia,2),2) > xDim(1) & nodalData(elementData(ia,2),2) < xDim(2) & ...
                        nodalData(elementData(ia,2),3) > yDim(1) & nodalData(elementData(ia,2),3) < yDim(2) ; 
       conditionalTwo = nodalData(elementData(ia,3),2) > xDim(1) & nodalData(elementData(ia,3),2) < xDim(2) & ...
                        nodalData(elementData(ia,3),3) > yDim(1) & nodalData(elementData(ia,3),3) < yDim(2) ; 
       
       
	   
       % If the element is entirely inside the region of interest, add it
       % to the calculation of fiber properties.
       if conditionalOne && conditionalTwo           
          % Structural anisotropy module:
          % Goal: Measure the direction of each element to determine the
          % average in-plane orientation
          %
          % Elements are simplified to linear pieces, which they are in
          % fact not (3-noded Timoshenko beam elements are quadratic)
          
          tempRand = rand(1);
          if tempRand > 0.5
             sIdx = 3; 
             eIdx = 2;
          else
             sIdx = 2;
             eIdx = 3;
          end
          
          elementVector = [ nodalData(elementData(ia,eIdx),2) - nodalData(elementData(ia,sIdx),2) ;
                            nodalData(elementData(ia,eIdx),3) - nodalData(elementData(ia,sIdx),3) ;
                            nodalData(elementData(ia,eIdx),4) - nodalData(elementData(ia,sIdx),4) ]; 
          elementVectorNorm = elementVector./norm(elementVector);
          elementVectorRed = elementVector(1:2)./norm(elementVector(1:2)); % In plane vector

          angYX(counter) = findAngleInPlane(elementVectorRed);

%           fiberProp(counter,:) = [LcTemp LpTemp realSetData(elementData(ia,5),3:5)];
          
          truncVector = [ nodalData(elementData(ia,truncSelStart),2) - nodalData(elementData(truncMarker,truncSelector),2) ;
                            nodalData(elementData(ia,truncSelStart),3) - nodalData(elementData(truncMarker,truncSelector),3) ;
                            nodalData(elementData(ia,truncSelStart),4) - nodalData(elementData(truncMarker,truncSelector),4) ]; 
          if tempRand > 0.5
              truncVector = - truncVector; 
          end
                        
                        
          if norm(truncVector)> 120e-6
              elementVectorTruncNorm = truncVector./norm(truncVector);
              elementVectorTruncRed = elementVectorTruncNorm(1:2)./norm(elementVectorTruncNorm(1:2));
              angYXTrunc(end+1) = findAngleInPlane(elementVectorTruncRed);
              fiberPropTrunc(end+1,:) = realSetData(elementData(ia,5),3:5);
              truncMarker = ia;
              truncSelector = 3;
              truncSelStart = 2;
          end
          
          
          
          counter = counter + 1;
          
%           plot(counter,rad2deg(angYX(counter)),'s','color',cMap(elementData(ia,5),:))
%           hold on
%           pause(0.5)
            
          
          
            if 0%1
               subplot(1,3,1)
               plot([nodalData(elementData(ia,2),2) nodalData(elementData(ia,3),2)],[nodalData(elementData(ia,2),3) nodalData(elementData(ia,3),3)],'s--')
               hold on
               plot([-inf inf],[0 0 ],'k--')
               plot([-0 0],[-inf inf],'k--')
               axis equal
               xlim(xDim); ylim(yDim);
               xlabel x; ylabel y;
                
                
               subplot(1,3,2) 
               plot([0 elementVectorRed(1)],[0 elementVectorRed(2)],'s-r') 
               hold on
               plot([0 1],[0 0],'k--','linewidth',2)
               plot(cos(angYX(1:(counter-2))), sin(angYX(1:(counter-2))),'b.')
               title(['Proposed Angle = ' num2str(rad2deg(angYX(counter-1)))])
               axis equal
               xlabel x; ylabel y;
               xlim([-1 1]); ylim([-1 1])
               pause(0.15)
               hold off
               
               subplot(1,3,3)
               polarhistogram(angYX(1:(counter-1)),360,'displaystyle','bar','Normalization','pdf');
               %  if mod(counter-1,100) == 0
%                pause(1)
%                end
            end
            
        else % If either start nor end is in the zone
          
       end
   end               
end

angYX(isnan(angYX)) = [];



% Create output structure which is useful for filtering operations
networkStructure.name = networkName;
networkStructure.dir = networkDir;
networkStructure.directionData = angYX;



ap = histcounts(mod(angYX,pi),linspace(0,pi,181));

networkStructure.FORatio = 0.5*(ap(1)+ap(end))/ap(90);


networkStructure.fiberProp = fiberProp;

networkStructure.lc = lc;
networkStructure.lp = lp;
networkStructure.width = width;
networkStructure.wallTkn = wallTkn;
networkStructure.curl = curl;


% Fiber orientation data
% Note that orientation is only calculated for 1 quadrant and then expanded
% via double-symmetry (A fiber oriented at 0 deg == a fiber oriented at 180
% deg).
% In-plane
networkStructure.angYX = networkStructure.directionData(:,1);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        figure();
        h = polarhistogram(networkStructure.angYX,360,'displaystyle','bar','Normalization','pdf');
        ang = h.BinEdges;
        angl = ang(1:360);
        v = h.Values;
        close;
    %   x= [v.*cos(angl) v.*cos(pi+angl)];
    %   y= [v.*sin(angl) v.*sin(pi+angl)];
    
    x= [v.*cos(angl)];
    y= [v.*sin(angl)];
    
    data = [x(:) y(:)];
    
    %==========================================================================
        %% finding the best fitted ellipse
        x=data(:,1);
        y=data(:,2);
        ellipse_t = fit_ellipse( x,y);
        X0_in= ellipse_t.X0_in;
        Y0_in= ellipse_t.Y0_in;
        a= ellipse_t.a;
        b= ellipse_t.b;
        if b>a
            a1 = a;
            b1 = b;
            a = b1;
            b = a1;
            
        end
        
        angle= ellipse_t.phi;
        

    
    
    [X,Y] = calculateEllipse(X0_in, Y0_in, a,b,angle*180/pi);
%     plot(X, Y, 'b','Linewidth',L);
    % savefig(N,['AnisoResults/EAnisoAngleStd',num2str(AngleStd(k))]);

    Aniso = 1- b/a;
    angle = angle*180/pi;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


networkStructure.Aniso = Aniso;
networkStructure.angle = angle;
networkStructure.a = a;
networkStructure.b = b;

















% Weighting module:
% Options
% 0. Not weighted
% 1. Volume weighted
% 2. Width weighted
% 3. Wall thickness weighted
% 
% weightedByWidth = weightedhistc(angYX, fiberProp(:,1), -pi/2:1/359:3*pi/2);
% weightedByNothing = weightedhistc(angYX, ones(size(fiberProp(:,3))), -pi/2:1/359:3*pi/2);
% weightedByWallThk = weightedhistc(angYX, fiberProp(:,5), 0:359:2*pi);

% if 0
%     figure; 
%     plot(-pi/2:1/359:3*pi/2,weightedByWidth./sum(fiberProp(:,1)))
%     hold on
%     plot(-pi/2:1/359:3*pi/2,weightedByNothing./sum(ones(size(fiberProp(:,3)))))
% end

% Now we actually need to convert the values into x and y, we cannot use
% the polarhistogram function.
% 
% The conversion is like this:
% 1 For each quadrant
% 2 Decompose circle into x and y component
% 3 Plot component

% cartValsNothing = [[weightedByNothing./sum(ones(size(fiberProp(:,3)))).*cos(-pi/2:1/359:3*pi/2)]' [weightedByNothing./sum(ones(size(fiberProp(:,3)))).*sin(-pi/2:1/359:3*pi/2)]'];
% cartValsWidth = [[weightedByWidth./sum(fiberProp(:,1)).*cos(-pi/2:1/359:3*pi/2)]' [weightedByWidth./sum(fiberProp(:,1)).*sin(-pi/2:1/359:3*pi/2)]'];


% networkStructure.cartValsNothing = cartValsNothing;
% networkStructure.cartValsWidth = cartValsWidth;

% 
% if 1
%     figure
%     plot(cartValsNothing(:,1),cartValsNothing(:,2),'-')
%     hold on
%     plot(cartValsWidth(:,1),cartValsWidth(:,2),'-')
% end
% 
% 
% 
% if 1
%    figure
%    subplot(1,2,1)
%    polarhistogram(angYX,360)
%    hold on
%    subplot(1,2,2)
%    polarhistogram(angYXTrunc',360)
%     
% end






end

function angYX = findAngleInPlane(elementVectorRed)
          % Find quadrant to look in
          if elementVectorRed(1) > 0 && elementVectorRed(2) > 0 % Here we can use any function
              angYX = acos(elementVectorRed(1));
          elseif elementVectorRed(1) > 0 && elementVectorRed(2) < 0 % Bottom right quadrant
              angYX = asin(elementVectorRed(2));
          elseif elementVectorRed(1) < 0 && elementVectorRed(2) > 0 % Top left quadrant
              angYX = acos(elementVectorRed(1));
          elseif elementVectorRed(1) < 0 && elementVectorRed(2) < 0 % Bottom left quadrant, here we need to 
              angYX = pi-asin(elementVectorRed(2));
          elseif elementVectorRed(1) == 0 && elementVectorRed(2) == 1
              angYX = pi/2;
          elseif elementVectorRed(1) == 1 && elementVectorRed(2) == 0
              angYX = 0;
          elseif elementVectorRed(1) == 0 && elementVectorRed(2) == -1
              angYX = -pi/2;
          elseif elementVectorRed(1) == -1 && elementVectorRed(2) == 0
              angYX = 180;
          else
              disp(stop)
          end
end

