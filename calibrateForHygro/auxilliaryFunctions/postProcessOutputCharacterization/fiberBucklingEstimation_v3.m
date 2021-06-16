% Meta instructions
clear; close all; clc
format compact; format long
[~, name] = system('hostname');
ctrl.workingDir = cd;
addpath(fullfile(ctrl.workingDir,'auxilliaryFunctions'))
addpath(fullfile(ctrl.workingDir,'@gramm'))

if ispc
    ctrl.fileSep = '\';
    ctrl.execEnvir = 'Windows';
else
    ctrl.fileSep = '/';
    
    ctrl.execEnvir = 'Linux';
end

% Decide on default settings for presentation/publication
[ctrl.refPos,custom] = initializeSettings('Publication','default');
close all;
ctrl.plotFlag = 1; 

targetPointer = 'L:\tempRadiusFebruary2019\File1\C\R_0.25';

targetPointer = 'C:\Users\augus\Desktop\tempBucklingPaperB';




disp('-> Importing nodal data.')
tic
nodalData = importNodalData(horzcat(targetPointer,ctrl.fileSep,'nodalData','.csv'));
toc

disp('-> Importing type data')
tic 
realSetData = importRealData(horzcat(targetPointer,ctrl.fileSep,'realData','.csv'));
toc

disp('-> Importing element data')
%disp('--> Observe the ordering: [INDEX CORNER 1 CORNER')
tic 
elementData = importElementData(horzcat(targetPointer,ctrl.fileSep,'elementData','.csv'));
toc

disp('-> Importing material data')
%disp('--> Observe the ordering: [INDEX CORNER 1 CORNER')
tic 
materialData = importMaterialData(horzcat(targetPointer,ctrl.fileSep,'materialData','.csv'));
toc

disp('Importing solution state energy and length')
tic 
fiberState = importFiberState(horzcat(targetPointer,ctrl.fileSep,'eState.csv'));
toc
fiberState(:,11) = 1.*(fiberState(:,6)>0); % Binary indicator of plasticity
fiberState(:,12) = 0;
fiberState(:,13) = 0;

% We need to truncate the data so that we are only looking in the region of
% interest. This is because the data close to the edges do not contain
% reliable information. 

xDim = 0.00035*[-1 1];
yDim = 0.003*[-1 1];
             
fiberData = zeros(realSetData(end,1),7);


newFiberMarker = diff(elementData(:,5));
newFiberMarker = [1 ; newFiberMarker];
newFiberIndex = find(newFiberMarker);


offsetCounter = 0;
for xLoop = 1:sum(newFiberMarker)-1
   %disp(['Working on fiber ',num2str(xLoop),'/',num2str(sum(newFiberMarker)-1)])
    
   startPos = newFiberIndex(xLoop);
   endPos = newFiberIndex(xLoop+1)-1;

   lReal = 0;
   
   yLoop = -1;
   while yLoop <=(endPos-startPos)-1
       yLoop = yLoop + 1; 
       ia = startPos+yLoop;
       
           
           lReal = lReal + sqrt((nodalData(elementData(ia,3),2)-nodalData(elementData(ia,2),2))^2 ...
                               +(nodalData(elementData(ia,3),3)-nodalData(elementData(ia,2),3))^2 ...
                               +(nodalData(elementData(ia,3),4)-nodalData(elementData(ia,2),4))^2);
       
   end
   
   lProj = sqrt((nodalData(elementData(endPos,3),2)-nodalData(elementData(startPos,2),2))^2 ...
               +(nodalData(elementData(endPos,3),3)-nodalData(elementData(startPos,2),3))^2 ...
               +(nodalData(elementData(endPos,3),4)-nodalData(elementData(startPos,2),4))^2);
   curvature = lReal/lProj;
                    
   fiberData(xLoop,:) = [lProj lReal realSetData(xLoop,2:5) curvature]; % Cross section type   
end
fiberData(fiberData(:,2)==0,:) = []; 

% This section handles the weighted data. The weighting is done by LENGTH
% or MASS.
fiberData(:,8) = 1430*fiberData(:,1).*fiberData(:,4).*fiberData(:,5)*1e-18;
fiberData(fiberData(:,3)==2,8) = 1430*fiberData(fiberData(:,3)==2,1).* ...
                                (fiberData(fiberData(:,3)==2,4)*2+fiberData(fiberData(:,3)==2,5)*2).* ...
                                 fiberData(fiberData(:,3)==2,6)*1e-18;
fiberData(:,9:15) = fiberData(:,1:7).*(fiberData(:,8)/mean(fiberData(:,8)));                   
grammage = 1000*sum(fiberData(:,8))/(range(xDim)*range(yDim))
                             
                             
                             

lambdaCalc = @(fData) [mean(fData) std(fData) median(fData) mode(fData)];

% Plot the average quantities
figure('name','fiberLength');
histogram(fiberData(:,1)*1e3,100);
hold on
histogram(fiberData(:,2)*1e3,100);
xlabel('Fiber length [mm]')
legend('Projected length','Real length','location','northeast')

tempFormatSpec = '%10s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n';
fprintf('%40s %45s\n','UNWEIGHTED','WEIGHTED BY MASS')
fprintf('%10s %10s %10s %10s %10s %10s %10s %10s %10s\n','','Mean','STD','Median','Mode','Mean','STD','Median','Mode')
fprintf(tempFormatSpec,'LProj',    lambdaCalc(fiberData(:,1)), lambdaCalc(fiberData(:,9)))
fprintf(tempFormatSpec,'LReal',    lambdaCalc(fiberData(:,2)), lambdaCalc(fiberData(:,10)))
fprintf(tempFormatSpec,'Width',    lambdaCalc(fiberData(:,4)), lambdaCalc(fiberData(:,12)))
fprintf(tempFormatSpec,'Height',   lambdaCalc(fiberData(:,5)), lambdaCalc(fiberData(:,13)))
fprintf(tempFormatSpec,'Thick',    lambdaCalc(fiberData(:,6)), lambdaCalc(fiberData(:,14)))
fprintf(tempFormatSpec,'Curvature',lambdaCalc(fiberData(:,7)),     lambdaCalc(fiberData(:,15)))

% figure('name','fiberWidth');
% histogram(fiberData(:,4),100);
% xlabel('Fiber width [\mum]')
% 
% figure('name','fiberHeight');
% histogram(fiberData(:,5),100);
% xlabel('Fiber height [\mum]')
% 
% figure('name','fiberThickness');
% histogram(fiberData(:,6),100);
% xlabel('Fiber wall thickness [\mum]')
% 
% figure('name','fiberCurl');
% histogram(fiberData(:,7),100);
% xlabel('Fiber curvature [-]')
%                                                  
%                                         




% Calculating the free span length
dexp = importBondData(horzcat(targetPointer,ctrl.fileSep,'DEXP.csv'));


% Make a comparison with analytic estimates.
nf = numel(unique(fiberData(:,1))); % Number of fibers in the network
rf = (mean(fiberData(:,4))+mean(fiberData(:,5)))/2*1e-6;     % Width (should actually be radius)
Lf = 1e-3;                          % Fiber length (hard coded because fibers used to calculate values above are truncated at the network edges)
t = 480e-6;                         % Network thickness
A = 700e-6*6000e-6;                 % Network in plane area (700 um x 6000 um)

NFC_35 = 2*Lf^2*rf*nf^2/(pi*A*t)





% The data recorded in dexp is:
% Column  1 : Contact element number
% Column  2 : Contact status (0 = unbonded, 1 = bonded)
% Column  3 : Contact damage (0 = undamaged, 1 = delaminated)
% Column  4 : X coordinate
% Column  5 : Y coordinate
% Column  6 : Z coordinate
% Column  7 : Fiber (real number) on the master side
% Column  8 : Fiber (real number) on the slave side
% Column  9 : Element number on the master side
% Column 10 : Element number on the slave side

% Basic idea to calculate the typical free span of a fiber:
%
% 1. Import all of the contact element information
% 2. Limit the analysis to bonds which are bonded (Column 2 == 1)
% 3. Limit the analysis to bonds which are unbroken (Column 3 < 1)
% 4. Tie together the bonded sites along the length of the fiber by
%    comparing their element numbers, noting that the element numbers are
%    always ascending along the length of the fiber. 
%     >>>>>>>DOESN'T SEEM TO BE TRUE?!
% 5. Extract the fiber cross sectional width through the real data 
% 6. Extract the orientation of the underlying beam elements to determine
%    whether the overlap is orthogonal or angled.
% 7. Sum the total fiber length.
% 8. Sum the total bonded length.
% 9. Compare the two using a histogram.

% Current status:
% 1. OK

% 2. + 3. Cleaning the data, removing unbonded and delaminated contact elements. 
dexpClean = dexp(dexp(:,2) == 1,:);
dexpClean = dexpClean(dexpClean(:,3) < 1,:);
disp(['Number of contacts changed by ',num2str(100*(size(dexpClean,1)-size(dexp,1))/size(dexp,1)),' %.' ])

for xLoop = 1:numel(unique(dexpClean(:,7)))
    dexpClean(dexpClean(:,7)==xLoop,11) = fiberData(xLoop,4); % Add the cross section width of the master fiber
    dexpClean(dexpClean(:,8)==xLoop,12) = fiberData(xLoop,4); % Add the cross section width of the slave fiber
end
% 4. Tie together the points


% figure('name','debugSingleFiber');
uniqueFibers = unique(dexpClean(:,7:8)); % OBS BOTH column 7 and 8 should somehow be incorporated here because both
                                       % target and contact elements count
                                       % as a contact.
                                       
                                       % Proposal for implementation:
                                       % reshape  the contact and target
                                       % columns (3 different ones) into
                                       % single columns (condensing 2 to 1
                                       % in all 3 cases).
% This is the step where I reshape dexpClean to a form where both target and contact contact is counted the same
%dexpCleanNew = [dexpClean(:,[1:6]) zeros(size(dexpClean(:,1))) zeros(size(dexpClean(:,1)))]

dexpCleanTemp = [dexpClean(:,8) dexpClean(:,7) dexpClean(:,10) dexpClean(:,9) dexpClean(:,12) dexpClean(:,11)];
dexpClean2 = [dexpClean ; dexpClean];
dexpClean2(size(dexpClean2,1)/2+1:end,7:12) = dexpCleanTemp;
dexpClean = dexpClean2;                            
                                       
typeOfContact = zeros(2,3);
tempDebug = [];

plotFlag = 0; % OBS ONLY PUT TO 1 IN DEBUG MODE!!!

% kk = figure('name','relErrorPerFiber');
for yLoop = 1:numel(uniqueFibers)-1
    evalFiber = uniqueFibers(yLoop);
    OneFiber(yLoop,1) = sum(dexpClean(:,7)==evalFiber);                     % Number of bonds, non length weighted
    OneFiber(yLoop,2) = OneFiber(yLoop,1)./fiberData(evalFiber,2)*1e-3;     % Length weighted per mm
    OneFiber(yLoop,3) = evalFiber;
    OneFiber(yLoop,4) = fiberData(evalFiber,2); % Fiber length
    OneFiber(yLoop,5) = 0; % "Free hanging length"
    
    % debug plotting
    if plotFlag
        figure();
        
        dataTemp = nodalData(elementData(elementData(:,5)==evalFiber,2:3),:);
        [~,ia,~] = unique(dataTemp(:,1));
        dataTemp = dataTemp(ia,:);
        
        plot3(dataTemp(:,2),dataTemp(:,3),dataTemp(:,4),'-o')
%         plot3(nodalData(elementData(elementData(:,5)==evalFiber,2),2), ...
%               nodalData(elementData(elementData(:,5)==evalFiber,2),3), ...
%               nodalData(elementData(elementData(:,5)==evalFiber,2),4),'-o')
        hold on
        axis equal
        %legend('Fiber','location','best')
    end
    
%     fprintf('Intersection analysis for fiber : %5d. Total length: %10.3f um. Total number of bonds: %3d.\n',evalFiber,fiberData(evalFiber,2),OneFiber(yLoop,1))

    % Figuring out the distance between different bond sites
    if OneFiber(yLoop,1) > 1 % If more than one bond along fiber
        fiberLength = fiberData(evalFiber,2); % Add fiber length
        
        % Calculate the free hanging ends of the fiber
    for zLoop = 1:OneFiber(yLoop,1)-1
      
        if yLoop == 1 && zLoop == 1
            OneFreeSpan(1,1) = evalFiber;
        else
            OneFreeSpan(end+1,1) = evalFiber; % Index of the fiber being observed
        end
        
        cordTemp1 = dexpClean(dexpClean(:,7)==evalFiber,:);
        
        % Sort the contact intersections so that we can walk along the
        % length of the fiber with monotonously increasing element number
        cordTemp1 = sortCoordinatesAlongElements(nodalData,elementData,cordTemp1,0);
        
        if plotFlag % Plot all the intersection points
            scatter3(cordTemp1(:,4),cordTemp1(:,5),cordTemp1(:,6), ...
                     10,1:size(cordTemp1,1),'filled'); axis equal;
        end        
        
        mT1 = cordTemp1(zLoop,4:6);       % Coordinate of the contact
        mT2 = cordTemp1(zLoop+1,4:6);     % Coordinate of the following contact
        eT1 = cordTemp1(zLoop,9);         % Element in the first contact
        eT2 = cordTemp1(zLoop+1,9);       % Element in the second contact
        
       if zLoop == 1
           % If we are at the first step, there is a section of the fiber before the first
           % bond that should be counted.
           
           % Step 1: Find how many elements are before the first element in
           % contact.
           elementsInFiber = elementData(elementData(:,5)==evalFiber,1);
           firstElementInFiber = elementsInFiber(1);
           lastElementInFiber = elementsInFiber(end);
           
           if firstElementInFiber == eT1
                sumLengthStart = 0;
           else % Sum the lengths of elements until you get the first bond
                sumLengthStart = sum(fiberState(firstElementInFiber:eT1-1,4));
           end
           
           elementCoordinatesT1 = nodalData(elementData(eT1,[2 3]),2:4);
           startFree = fiberState(eT1,4) - max(0,calculateElementSegment(elementCoordinatesT1,mT1));
           OneFiber(yLoop,5) = OneFiber(yLoop,5) + sumLengthStart + startFree;
       end   
        if zLoop == OneFiber(yLoop,1)-1
            if lastElementInFiber == eT2
                % If the first element has a part in contact, then just take the contribution from the fiber end that is outside that.
                sumLengthFinish = 0;
            else
                % Sum the lengths of elements until you get the first bond,
                % + the partial contact.
                sumLengthFinish = sum(fiberState(eT2+1:lastElementInFiber,4));
            end
            elementCoordinatesT2 = nodalData(elementData(eT2,[2 3]),2:4);
            endFree = sumLengthFinish + calculateElementSegment(elementCoordinatesT2,mT2);
            OneFiber(yLoop,5) = OneFiber(yLoop,5) + endFree;
        end
        
                
        if eT1 == eT2 
            % If both contacts occur over the same element
            % Method: Use 2-norm of the contact positions to calculate free
            % span.
            
            elementCoordinatesT1 = nodalData(elementData(eT1,[2 3]),2:4); % XYZ position of end nodes in eT1
            OneFreeSpan(end,2) = calculateElementSegment(elementCoordinatesT1,[mT1 ; mT2]);
            
            OneFreeSpan(end,16) = sum(fiberState(eT1,7:10))*OneFreeSpan(end,2)/fiberState(eT1,4);
                              
            typeOfContact(1,1) = typeOfContact(1,1)+1;
            typeOfContact(2,1) = typeOfContact(2,1)+OneFreeSpan(end,2);

            if plotFlag
                plot3(elementCoordinatesT1(:,1),elementCoordinatesT1(:,2),elementCoordinatesT1(:,3),'-ko');
                hold on; xlabel x; ylabel y; zlabel z; 
                plot3(mT1(:,1),mT1(:,2),mT1(:,3),'-sb');
                plot3(mT2(:,1),mT2(:,2),mT2(:,3),'-sr');
                title('Inside element')
            end
            
            TCTEMP = 1;
        elseif abs(eT1-eT2) == 1
            % If the contacts are on consecutive elements
            % Method: Calculate the element segment in one element, and the
            % element segment from the second element, and then sum.
            elementCoordinatesT1 = nodalData(elementData(eT1,[2 3]),2:4); % XYZ position of end nodes in eT1
            elementCoordinatesT2 = nodalData(elementData(eT2,[2 3]),2:4); % XYZ position of end nodes in eT1
            
            
            % Calculate the partial contributions
            freeSpanTemp1 = calculateElementSegment(elementCoordinatesT1,mT1);        %max(0,calculateElementSegment(elementCoordinatesT1,mT1));
            freeSpanTemp2 = calculateElementSegment(flip(elementCoordinatesT2,1),mT2);%   0     %max(0,calculateElementSegment(elementCoordinatesT2,mT2));   
            OneFreeSpan(end,2) = freeSpanTemp1 + freeSpanTemp2;
            typeOfContact(1,2) = typeOfContact(1,2)+1;
            typeOfContact(2,2) = typeOfContact(2,2)+OneFreeSpan(end,2);
            
            OneFreeSpan(end,16) = sum(fiberState(eT1,7:10))*freeSpanTemp1/fiberState(eT1,4) + ...
                                  sum(fiberState(eT2,7:10))*freeSpanTemp2/fiberState(eT2,4);            
            %figure(); 
            if plotFlag
                plot3(elementCoordinatesT1(:,1),elementCoordinatesT1(:,2),elementCoordinatesT1(:,3),'-o');
                hold on; xlabel x; ylabel y; zlabel z; 
                plot3(elementCoordinatesT2(:,1),elementCoordinatesT2(:,2),elementCoordinatesT2(:,3),'-o');
                plot3(mT1(:,1),mT1(:,2),mT1(:,3),'sb');
                plot3(mT2(:,1),mT2(:,2),mT2(:,3),'sr');
                axis equal
                title('Between 2 elements')
            end
            TCTEMP = 2;
        else
            % If they are not on the same element
            % Method: Count the elements between the two, and sum their
            % lengths. Furthermore, count the partial lengths coming from
            % the two elements where the contacts occur. 
            
            elementCoordinatesT1 = nodalData(elementData(eT1,[2 3]),2:4); % XYZ position of end nodes in eT1
            elementCoordinatesT2 = nodalData(elementData(eT2,[2 3]),2:4); % XYZ position of end nodes in eT2
            elementCoordinatesTX = nodalData(elementData(eT1+1:eT2-1,[2]),2:4);% XYZ position of nodes between eT1 and eT2
            elementCoordinatesTX(end+1,:) = nodalData(elementData(eT2-1,[3]),2:4);

            % Calculate the partial contributions
            freeSpanTemp1 = calculateElementSegment(elementCoordinatesT1,mT1);        %max(0,calculateElementSegment(elementCoordinatesT1,mT1));
            freeSpanTemp2 = calculateElementSegment(flip(elementCoordinatesT2,1),mT2);       
            
            sumWholeElements = sum(fiberState(eT1+1:eT2-1,4));
            
            OneFreeSpan(end,2) = sumWholeElements + freeSpanTemp1 + freeSpanTemp2;
            
            OneFreeSpan(end,16) = sum(sum(fiberState(eT1+1:eT2-1,7:10))) + ...
                                  sum(fiberState(eT1,7:10))*freeSpanTemp1/fiberState(eT1,4) + ...
                                  sum(fiberState(eT2,7:10))*freeSpanTemp2/fiberState(eT2,4); % Free span elastic energy.
            
            %figure(); 
            if plotFlag
                plot3(elementCoordinatesT1(:,1),elementCoordinatesT1(:,2),elementCoordinatesT1(:,3),'-o');
                hold on; xlabel x; ylabel y; zlabel z; 
                plot3(elementCoordinatesT2(:,1),elementCoordinatesT2(:,2),elementCoordinatesT2(:,3),'-o');
                plot3(elementCoordinatesTX(:,1),elementCoordinatesTX(:,2),elementCoordinatesTX(:,3),'-^');
                plot3(mT1(:,1),mT1(:,2),mT1(:,3),'sb');
                plot3(mT2(:,1),mT2(:,2),mT2(:,3),'sr'); 
                axis equal
                title('Many elements')
            end
           typeOfContact(1,3) = typeOfContact(1,3)+1;
           typeOfContact(2,3) = typeOfContact(2,3)+OneFreeSpan(end,2);
           TCTEMP = 3;
        end
        
        % Add in other physical characteristics of the free span and
        % calculate the second area moment of inertia, Iy.
        OneFreeSpan(end,3:6) = realSetData(evalFiber,2:5);     % This is the free span cross section height, width, wall thickness and cross section shape.
        OneFreeSpan(end,7) = calculateIy(OneFreeSpan);

         
         % This part figures out what material properties should be used in
         % the free span. To do this:
         % 
         % 1. Figure out which elements are in the free span.
         % 2. Figure out if any of those elements have dissipated plastic
         %    energy so far.
         %    IF YES:   Use the hardening modulus.
         %    IF NO :   Use the elastic modulus.
         if sum(fiberState(eT1:eT2,11)) > 0
             OneFreeSpan(end,8) = 5000e6; % Hardening modulus, in Pa
         else
             OneFreeSpan(end,8) = materialData(elementData(eT1,6),2)*1e6; % Elastic modulus, in Pa
         end
         
         
         
         % This part adjusts the estimate of the free span by accounting
         % for the fact that the beams have a finite width, and thus that
         % width needs to be subtracted from the total free span.
         %
         % 1. Find out the width of the crossing fiber.
         % 2. Find out the angle between the two elements.
         % 3. Calculate the "bond adjacent zone"
         % 4. Remove half of this amount from the free span (one half of
         %    the bond is presumably on the other side, in another free span).

         widthToRemove1 = realSetData(cordTemp1(zLoop,7),3)/2;
         widthToRemove2 = realSetData(cordTemp1(zLoop,8),3)/2;
         
         OneFreeSpan(end,9) =max(0,OneFreeSpan(end,2) - widthToRemove1 - widthToRemove2);                   % introduce MAX here, this is the edited free span now
         OneFreeSpan(end,10) = pi^2*OneFreeSpan(end,8)*OneFreeSpan(end,7)*1e-24/(1e-6*OneFreeSpan(end,9))^2;% Euler buckling load
         OneFreeSpan(end,11) = -min(fiberState(eT1:eT2,3))*1e-6;                                             % Actual load in the segment, in Newton
         OneFreeSpan(end,12) = min(1e6,OneFreeSpan(end,10));
         
         
         
         facTemp = OneFreeSpan(end,12)/OneFreeSpan(end,11);
         if facTemp < 0
             OneFreeSpan(end,13) = 1e10;
         else
             OneFreeSpan(end,13) = facTemp;
         end
         
         OneFreeSpan(end,14) = calculateCrossSectionArea(OneFreeSpan(end,:))*OneFreeSpan(end,2); % Free span mass (actually volume).
         OneFreeSpan(end,15) = calculateCrossSectionArea(OneFreeSpan(end,:))*OneFreeSpan(end,9); % Free span mass, adjusted for fiber width (actually volume).
         OneFreeSpan(end,17) = OneFreeSpan(end,16)*(OneFreeSpan(end,9)/OneFreeSpan(end,2));
         
         %OneFreeSpan(end,13) = min(1e10,OneFreeSpan(end,12)/OneFreeSpan(end,11)); % Factor increase needed to buckle
         
         % Associate the element with buckling marker (yes/no)
         fiberState(eT1:eT2,13) = max(fiberState(eT1:eT2,13),OneFreeSpan(end,11)/OneFreeSpan(end,12)); % Marker for buckled state (OBS SHOULD BE MADE CONSERVATIVE)
         
         
        %fprintf('Bond number %3d. Type: %3d. Free length: %10.3f um.\n',zLoop,TCTEMP,OneFreeSpan(end,2))
    end
    elseif OneFiber(yLoop,1) == 1
        OneFiber(yLoop,5) = OneFiber(yLoop,4); %fiberData(evalFiber,2)
    elseif OneFiber(yLoop,1) == 0
        OneFiber(yLoop,5) = OneFiber(yLoop,4); %fiberData(evalFiber,2)
    end
    
%     fprintf('Fiber %5d completed.                   Total length: %10.3f um. Diff against expected length: %7.2f percent.\n', ...
%              evalFiber,OneFiber(yLoop,5)+sum(OneFreeSpan(OneFreeSpan(:,1)==evalFiber,2)),100*(OneFiber(yLoop,5)+sum(OneFreeSpan(OneFreeSpan(:,1)==evalFiber,2))-fiberData(evalFiber,2))/fiberData(evalFiber,2))
%     
    
         
%     figure(kk)
%     plot(evalFiber,100*((OneFiber(yLoop,5)+sum(OneFreeSpan(OneFreeSpan(:,1)==evalFiber,2))-fiberData(evalFiber,2))/fiberData(evalFiber,2)),'o','color',custom.Color(1,:))
%     hold on
%     xlabel('Fiber evaluated')
%     ylabel('Relative error [%]')
%     pause(0.1)
    if mod(yLoop,100) == 0
        fprintf('%5d / %5d fibers examined, in of which %4.1f percent of the free spans have yielded.\n',yLoop,numel(uniqueFibers),100*length(OneFreeSpan(OneFreeSpan(:,8)==5000e6,9))/length(OneFreeSpan(:,9)))
        pause(0.1);
    end
end

fprintf('Length of free spans: %10.3f um. Length of free hanging: %10.3f um. Length of fibers: %10.3f um.\n',sum(OneFreeSpan(:,2)),sum(OneFiber(:,5)),sum(fiberData(:,2)))

OneFreeSpan(OneFreeSpan(:,2)==0,:) = [];


figure('name','freeSpansInNetwork');%,'Units','centimeters','Position',custom.FigureSize.*[1 1 2 1]+[10 10 0 0],'PositionMode','manual');
f3 = gramm('x',{OneFreeSpan(:,2),OneFreeSpan(:,9)},'color',{'Pointwise free spans','Accounting for fiber width'});    
%f3.stat_bin('normalization','pdf','geom','stairs','nbins',700);
f3.stat_bin('normalization','cdf','geom','stairs','edges',0:1:100);%stat_bin('geom','overlaid_bar');
f3.set_text_options('font','Times New Roman', ...
                        'interpreter','tex' , ...
                        'base_size',18);
f3.set_layout_options('legend_position',[0.4 0.2 0.3 0.4], ...
                          'redraw',true);
% f3.set_color_options('legend','merge');
f3.axe_property('xlim',[0 100],'ylim',[0 1]);
f3.set_names('x','Estimated free span between bond sites [\mum]','y','EDF','color','');                      
f3.draw();
set(gca,'XColor','k');
set(gca,'YColor','k');  


f4Legend = {horzcat('Elastic fiber (',num2str(round(100*length(OneFreeSpan(OneFreeSpan(:,8)~=5000e6,9))/length(OneFreeSpan(:,9)),1)),'% of fibers)') , ...
            horzcat('Yielding fiber (',num2str(round(100*length(OneFreeSpan(OneFreeSpan(:,8)==5000e6,9))/length(OneFreeSpan(:,9)),1)),'% of fibers)')};
figure('name','factorNecessaryToBuckle');%,'Units','centimeters','Position',custom.FigureSize.*[1 1 2 1]+[10 10 0 0],'PositionMode','manual');
f4 = gramm('x',{OneFreeSpan(OneFreeSpan(:,8)~=5000e6,13),OneFreeSpan(OneFreeSpan(:,8)==5000e6,13)},'color',f4Legend);    
f4.stat_bin('normalization','cdf','geom','stairs','edges',[0:0.1:100 110:10:1000 1100:100:1e6]);%,'edges',0:1:100);%stat_bin('geom','overlaid_bar');
f4.set_text_options('font','Times New Roman', ...
                        'interpreter','tex' , ...
                        'base_size',18);
f4.set_layout_options('legend_position',[0.2 0.8 0.3 0.2], ...
                          'redraw',true);
%f4.set_color_options('legend','merge');
f4.axe_property('xlim',[0 1e6],'ylim',[0 0.25],'XScale','log');
f4.set_names('x','\eta','y','Cumulative fraction of free spans','color','');                      
% f4.axe_property() 
f4.draw();
set(gca,'XColor','k');
set(gca,'YColor','k');  
% set(gca, 'XScale', 'log')



% Time to make a weighted histogram.
h1 = weightedhistc(OneFreeSpan(OneFreeSpan(:,8)~=5000e6,13), OneFreeSpan(OneFreeSpan(:,8)~=5000e6,15),[0:0.1:100 110:10:1000 1100:100:1e6]);
h2 = weightedhistc(OneFreeSpan(OneFreeSpan(:,8)==5000e6,13), OneFreeSpan(OneFreeSpan(:,8)==5000e6,15),[0:0.1:100 110:10:1000 1100:100:1e6]);
h3 = weightedhistc(OneFreeSpan(OneFreeSpan(:,8)~=5000e6,13), OneFreeSpan(OneFreeSpan(:,8)~=5000e6,17),[0:0.1:100 110:10:1000 1100:100:1e6]);
h4 = weightedhistc(OneFreeSpan(OneFreeSpan(:,8)==5000e6,13), OneFreeSpan(OneFreeSpan(:,8)==5000e6,17),[0:0.1:100 110:10:1000 1100:100:1e6]);

% figure(); stairs([0:0.1:100 110:10:1000 1100:100:1e6],cumsum(h1)./sum(OneFreeSpan(OneFreeSpan(:,8)~=5000e6,2)));
% hold on; stairs([0:0.1:100 110:10:1000 1100:100:1e6],cumsum(h2)./sum(OneFreeSpan(OneFreeSpan(:,8)==5000e6,2)));
% set(gca,'xscale','log')
% 

f6 = gramm('x',{[0:0.1:100 110:10:1000 1100:100:1e6],[0:0.1:100 110:10:1000 1100:100:1e6]}, ...
           'y',{cumsum(h1)./sum(OneFreeSpan(OneFreeSpan(:,8)~=5000e6,14)),cumsum(h2)./sum(OneFreeSpan(OneFreeSpan(:,8)==5000e6,14))},'color',f4Legend);%,h2./sum(OneFreeSpan(OneFreeSpan(:,8)==5000e6,2))},'color',f4Legend);
figure('name','factorNecessaryToBuckle2');%,'Units','centimeters','Position',custom.FigureSize.*[1 1 2 1]+[10 10 0 0],'PositionMode','manual');
%f6.stat_bin('geom','stairs','edges',[0:0.1:100 110:10:1000 1100:100:1e6]);%,'edges',0:1:100);%stat_bin('geom','overlaid_bar');
f6.geom_line();
f6.set_text_options('font','Times New Roman', ...
                        'interpreter','tex' , ...
                        'base_size',18);
f6.set_layout_options('legend_position',[0.2 0.8 0.3 0.2], ...
                          'redraw',true);
%f4.set_color_options('legend','merge');
f6.axe_property('xlim',[0 1e6],'ylim',[0 0.25],'XScale','log');
f6.set_names('x','\eta','y','Cumulative fraction of mass','color','');                      
% f4.axe_property() 
f6.draw();
set(gca,'XColor','k');
set(gca,'YColor','k');  



% Energy weighted histogram
f7 = gramm('x',{[0:0.1:100 110:10:1000 1100:100:1e6],[0:0.1:100 110:10:1000 1100:100:1e6]}, ...
           'y',{cumsum(h3)./sum(OneFreeSpan(OneFreeSpan(:,8)~=5000e6,16)),cumsum(h4)./sum(OneFreeSpan(OneFreeSpan(:,8)==5000e6,16))},'color',f4Legend);%,h2./sum(OneFreeSpan(OneFreeSpan(:,8)==5000e6,2))},'color',f4Legend);
figure('name','factorNecessaryToBuckleEnergy');%,'Units','centimeters','Position',custom.FigureSize.*[1 1 2 1]+[10 10 0 0],'PositionMode','manual');
%f6.stat_bin('geom','stairs','edges',[0:0.1:100 110:10:1000 1100:100:1e6]);%,'edges',0:1:100);%stat_bin('geom','overlaid_bar');
f7.geom_line();
f7.set_text_options('font','Times New Roman', ...
                        'interpreter','tex' , ...
                        'base_size',18);
f7.set_layout_options('legend_position',[0.2 0.8 0.3 0.2], ...
                          'redraw',true);
%f4.set_color_options('legend','merge');
f7.axe_property('xlim',[0 1e6],'ylim',[0 0.25],'XScale','log');
f7.set_names('x','\eta','y','Cumulative fraction of elastic energy','color','');                      
% f4.axe_property() 
f7.draw();
set(gca,'XColor','k');
set(gca,'YColor','k');  



% Back-calculate the number of defects per fiber length.
factorsToTest = [1:100 110:10:1000 1100:100:10000 11000:1000:1e5];
defectsSave = [];
for mLoop = 1:length(factorsToTest)
    facTemp = factorsToTest(mLoop);
    numberOfBuckles = sum(OneFreeSpan(:,13)<facTemp);
    totalFiberLength = sum(OneFreeSpan(:,9));%+sum(OneFiber(:,5)); % OBS decision needs to be made: All length or only free span length
    defectsPerMm = numberOfBuckles/totalFiberLength;
    defectsSave = [defectsSave ; defectsPerMm];
end
f12 = gramm('x',{factorsToTest,[0 1e5] , [0 1e5] , [0 1e5]},'y',{defectsSave',[0.048 0.048] , [0.078 0.078] , [0.085 0.085]}, ...
            'color',{'Our data','Mill kraft [67]','Lab. kraft [67]','Holocellulose [67]'});
figure('name','defectsPerMmFreespan');
f12.geom_line();
% f12.geom_hline('yintercept',0.048);
% f12.geom_hline('yintercept',0.085);
% f12.geom_hline('yintercept',0.078);
f12.set_text_options('font','Times New Roman', ...
                        'interpreter','tex' , ...
                        'base_size',18);
f12.set_layout_options('legend_position',[0.5 0.2 0.3 0.4], ...
                          'redraw',true);
f12.axe_property('xlim',[0 1e5],'ylim',[0 0.1]);
f12.set_names('x','\eta','y','Buckled segments per mm fiber','color','');                      
% f4.axe_property() 
f12.draw();
set(gca,'XColor','k');
set(gca,'YColor','k');  


% plot(factorsToTest,defectsSave,'xb-')
% hold on   
% plot([0 1e4],0.078*[1 1],'k-')
% plot([0 1e4],0.048*[1 1],'k-')
% plot([0 1e4],0.085*[1 1],'k-')
% xlabel('factor \eta')
% ylabel('Buckled segments per mm fiber')
% %set(gca,'xscale','log')
% legend('Our data','Iribarne','location','best')









% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % NEXT SECTION: INVESTIGATING WITH THE AIM OF FINDING WHERE IN THE NETWORK
% %               THE HIGHLY LOADED SEGMENTS ARE.
% %
% 
% % Axial load in each element given by the element data in fiberState
% % Coordinates given in the nodalData(elementData(:,1)==evalFiber,2:4)
% % Plot on a sheet or 
% 
% elementsToCheck = unique(fiberState(:,1));
% elementsToCheck(elementsToCheck==0) = [];
% 
% minColorVal = min(fiberState(:,3)); % Minimum value in the set
% maxColorVal = max(fiberState(:,3)); % Maximum value in the set
% meanColorVal = mean(fiberState(:,3)); % Maximum value in the set
% % discColor1 = linspace(-2.5e4,0,126);
% % discColor2 = maxColorVal;
% % discColor3 = minColorVal;
% %discColor2 = linspace(meanColorVal,maxColorVal,96);
% %discColor = [discColor3 discColor1 discColor2];
% discColor = [ linspace(0,1,127)' ; max(fiberState(:,8))];
% %discColorTop = max(fiberState(:,8));
% %discColor = [discColor1];
% indexColor = 1:128;
% colorsToUse = [flip((pink(128))) linspace(0,1,128)']; % 128 levels of colors.
% 
% figure('name','fibersAtRiskOfBuckling');
% set(gcf,'Renderer', 'openGL');
% for pLoop = 1:length(elementsToCheck)
%     evEle = elementsToCheck(pLoop);
%     
%     % Nodal coordinates
%     nX = nodalData(elementData(elementData(:,1)==evEle,2:3),2);
%     nY = nodalData(elementData(elementData(:,1)==evEle,2:3),3);
%     nZ = nodalData(elementData(elementData(:,1)==evEle,2:3),4);
%     
%     cE = round(interp1(discColor,indexColor,fiberState(evEle,8))); %fiberState(evEle,3),linspace(minColorVal,maxColorVal,128));
%     % Fiber being evaluated
%     if pLoop == 1
%         plot3(nX,nY,nZ,'-','color',colorsToUse(cE,:))
%         hold on
%         xlabel x; ylabel y; zlabel z;
%         axis equal
%     elseif fiberState(evEle,8) > 1e-3
%         line(nX,nY,nZ,'color',colorsToUse(cE,:))
%     else
%     end
%     if mod(pLoop,100) == 0
%         pause(0.05);
%     end
% end
% 
% colormap(flip(pink(128)));
% colorbar
% 
% 
% figure('name','fibersInNetwork');
% set(gcf,'Renderer', 'openGL');
% for pLoop = 1:length(elementsToCheck)
%     evEle = elementsToCheck(pLoop);
%     
%     % Nodal coordinates
%     nX = nodalData(elementData(elementData(:,1)==evEle,2:3),2);
%     nY = nodalData(elementData(elementData(:,1)==evEle,2:3),3);
%     nZ = nodalData(elementData(elementData(:,1)==evEle,2:3),4);
%     
%     %cE = round(interp1(discColor,indexColor,fiberState(evEle,8))); %fiberState(evEle,3),linspace(minColorVal,maxColorVal,128));
%     % Fiber being evaluated
%     if pLoop == 1
%         plot3(nX,nY,nZ,'-','color',colorsToUse(64,:))
%         hold on
%         xlabel x; ylabel y; zlabel z;
%         axis equal
%     elseif fiberState(evEle,8) > 1e-3
%         line(nX,nY,nZ,'color',colorsToUse(64,:))
%     else
%     end
%     if mod(pLoop,100) == 0
%         pause(0.05);
%     end
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %     [N1,edges] = histcounts(OneFreeSpan(OneFreeSpan(:,8)~=5000e6,9),0:0.01:100);
% %     [N2,edges] = histcounts(OneFreeSpan(OneFreeSpan(:,8)==5000e6,9),0:0.01:100);
% %     [N3,edges2] = histcounts(OneFreeSpan(OneFreeSpan(:,8)~=5000e6,10));
% %     [N4,edges2] = histcounts(OneFreeSpan(OneFreeSpan(:,8)==5000e6,10),edges2);
% %     subplot(1,3,1)
% %     stairs(diff(edges2)+edges2(1:end-1),1-cumsum(N3)/sum(N3))
% %     hold on
% %     stairs(diff(edges2)+edges2(1:end-1),1-cumsum(N4)/sum(N4))
% %     legend('Elastic fiber with force above this','Plastic fiber with force above this','location','best')
% %     xlabel('Force [N]')
% %     ylabel('1-EDF [-]')
% %     %set(gca, 'XScale', 'log')
% %     hold off
% %     
% %     subplot(1,3,2)
% %     %disp(['Fibers counted: ',num2str(yLoop),'/',num2str(numel(uniqueFibers))])
% %     stairs(diff(edges)+edges(1:end-1),cumsum(N1)/sum(N1))
% %     hold on
% %     stairs(diff(edges)+edges(1:end-1),cumsum(N2)/sum(N2))
% %     xlabel('Force [N]')
% %     ylabel('EDF [-]')
% %     set(gca, 'XScale', 'log')
% %     legend('Buckling force, elastic','Buckling force, plastic','location','best')
% %     hold off
% %     
% %     subplot(1,3,3)
% %     %disp([num2str(length(OneFreeSpan(OneFreeSpan(:,8)~=5000e6,9))/length(OneFreeSpan(:,9))),' of fibers have yielded'])
% %     [N5,edges5] = histcounts(OneFreeSpan(OneFreeSpan(:,8)~=5000e6,9)./OneFreeSpan(OneFreeSpan(:,8)~=5000e6,10),0:0.1:1000);
% %     [N6,edges6] = histcounts(OneFreeSpan(OneFreeSpan(:,8)==5000e6,9)./OneFreeSpan(OneFreeSpan(:,8)==5000e6,10),edges5);
% %     stairs(diff(edges6)+edges6(1:end-1),cumsum(N5)/sum(N5))
% %     hold on
% %     stairs(diff(edges6)+edges6(1:end-1),cumsum(N6)/sum(N6))
% %     legend('Factor the elastic modulus would need to be decreased to cause buckling')
% %     xlabel('Factor [-]')
% %     set(gca, 'XScale', 'log')
% %     %xlim([0 1000])
% %     hold off
% 
% % 
% % sum(OneFreeSpan(OneFreeSpan(:,13)<0.1,9))/sum(OneFreeSpan(:,9))
% % ans =
% %    0.109318467628419
% % sum(OneFreeSpan(OneFreeSpan(:,13)<1,9))/sum(OneFreeSpan(:,9))
% % ans =
% %    0.143055630766062
% % sum(OneFreeSpan(:,13)<1)/length(OneFreeSpan(:,9))
% % ans =
% %    0.116016387793391
% % sum(OneFreeSpan(:,13)<1)/length(OneFreeSpan(:,9))
% %sum(OneFreeSpan(:,9)>0)/sum(fiberState(:,4))*1e3
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(); 
% plot(1,sum(OneFreeSpan(OneFreeSpan(:,13)<1,9))./sum(OneFreeSpan(:,2)),'o')
% hold on
% plot(3,sum(OneFreeSpan(OneFreeSpan(:,13)<3,9))./sum(OneFreeSpan(:,2)),'o')
% plot(5,sum(OneFreeSpan(OneFreeSpan(:,13)<5,9))./sum(OneFreeSpan(:,2)),'o')
% plot(10,sum(OneFreeSpan(OneFreeSpan(:,13)<10,9))./sum(OneFreeSpan(:,2)),'o')
% plot(20,sum(OneFreeSpan(OneFreeSpan(:,13)<20,9))./sum(OneFreeSpan(:,2)),'o')
% plot(100,sum(OneFreeSpan(OneFreeSpan(:,13)<100,9))./sum(OneFreeSpan(:,2)),'o')
% plot(1000,sum(OneFreeSpan(OneFreeSpan(:,13)<1000,9))./sum(OneFreeSpan(:,2)),'o')
% 
% 
% plot(1,sum(OneFreeSpan(:,13)<1)./length(OneFreeSpan(:,2)),'x')
% hold on
% plot(3,sum(OneFreeSpan(:,13)<3)./length(OneFreeSpan(:,2)),'x')
% plot(5,sum(OneFreeSpan(:,13)<5)./length(OneFreeSpan(:,2)),'x')
% plot(10,sum(OneFreeSpan(:,13)<10)./length(OneFreeSpan(:,2)),'x')
% plot(20,sum(OneFreeSpan(:,13)<20)./length(OneFreeSpan(:,2)),'x')
% plot(100,sum(OneFreeSpan(:,13)<1e2)./length(OneFreeSpan(:,2)),'x')
% plot(1e3,sum(OneFreeSpan(:,13)<1e3)./length(OneFreeSpan(:,2)),'x')

savePlots(ctrl)