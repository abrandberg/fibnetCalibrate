function [refPos,custom] = initializeSettings(gfxMode,colorFlag)
%initializeSettings(gfxMode,colorFlag) sets the default figure and plot
%properties for the current MATLAB session.
%
% INPUT gfxMode:   String to specify whether to format for presentation or
%                  publication.
%       colorFlag: String to specify whether black & white or a color
%                  scheme is desired.
%
% OUTPUT RefPos: Reference position for the window. Can be used to reset
%                the window position without restarting the session (the 
%                default position is not saved by MATLAB).
%        custom: Structure containing colors and markers.
%
% REMARKS:
% - The applied changes last until the session is ended.
% - The file was originally created for LS-Dyna Paper Fibre Bonds article.
%
% TO DO:
%  
% 
% created by: August Brandberg
% DATE: 02-02-2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the reference position for plots (middle of the screen)
refPos = get(0,'defaultFigurePosition');


% Set the correct sizes for everything depending on whether graphs for
% PowerPoint or Word are required
if strcmp(gfxMode,'Presentation')
    set(groot,'defaultLineLineWidth' ,1.50);
    set(groot,'defaultAxesFontSize'  ,18  );
    set(groot,'defaultLegendFontSize',18  );
    set(groot,'defaultLineMarkerSize',8   );
    set(groot,'defaultAxesLineStyleOrder',{'-'});
    
elseif strcmp(gfxMode,'Publication')
    % Axes
    set(groot,'defaultAxesBox','off');
    set(groot,'defaultAxesFontName','Times New Roman');
    set(groot,'defaultAxesFontSize'  ,14  );
    set(groot,'defaultAxesXColor','k');
    set(groot,'defaultAxesYColor','k');
    set(groot,'defaultAxesZColor','k');
    %set(groot,'defaultAxesFontWeight','bold');
    set(groot,'defaultAxesLineStyleOrder',{'-'}); 
    %set(groot,'defaultAxesTickLabelFormat','%.2e'); Implemented in 2017a
    % Figure
    set(groot,'defaultFigureColor', 'white');
    custom.FigureSize = [0 0 12 9.5];
    set(groot,'defaultFigurePaperPosition', custom.FigureSize);
    set(groot,'defaultFigurePaperPositionMode', 'manual');
    set(groot,'defaultFigurePaperUnits', 'centimeters');
    set(groot,'defaultFigureRenderer', 'painter');
    % Legends
    set(groot,'defaultLegendEdgeColor',[1 1 1]); 
    set(groot,'defaultLegendFontSize',10  );
    set(groot,'defaultLegendColor','white');
    % Lines
    set(groot,'defaultLineLineWidth' ,2);
    set(groot,'defaultLineMarkerSize',6); %6 i vanliga fall
    %set(groot,'defaultLineMarkerFaceColor','on');
    
end

% Set color and marker properties and return the structure CUSTOM which
% contains the information
if strcmp(colorFlag,'BlackWhite')
    Color = [0 0 1; 1 0 0];
elseif strcmp(colorFlag,'default')
    Color = lines(50);
elseif strcmp(colorFlag,'jet2')
    ColorTemp = jet(20);
    Color = ColorTemp([1 20 6 16],:);
    
else
    Color = colormap(parula(20)); 
%     Color = [        0    0.4470    0.7410
%                 0.4940    0.1840    0.5560
%                 0.6350    0.0780    0.1840];
    %Color = lines(5);
    close all
end
Marker = {'o','^','s','x','*','d','+','v','>','<','p','h'};
%Marker = {'o','^','s'};%,'.','.'};
set(groot,'defaultAxesColorOrder',Color);
custom.Color = Color;
custom.Marker = Marker;

