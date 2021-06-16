function plotPulp(ctrl,lc,width,wallTkn,curl,fibrillation,lc2,width2,wallTkn2,curl2,fibrillation2)

if nargin == 6
    lc2 = nan;
    width2 = nan;
    wallTkn2 = nan;
    curl2 = nan;
    fibrillation2 = nan;
end




saveExtensions = {'-dpng','-depsc'};
ctrl.saveName = '';
saveName = [ctrl.saveDir filesep 'overview_syn' ctrl.saveName];
if not(2== exist([saveName '.' saveExtensions{1}(3:end)],'file'))
    corrMatrix = corr([lc,width,wallTkn,curl,fibrillation],'type','kendall'); % 

    histogramInstructions = {'normalization','pdf','edgecolor',[1 1 1]};
    lineInstructions = {'markersize',0.2};
    textInstructions = {'horizontalAlignment','center','fontsize',22,'interpreter',ctrl.interpreter};
    numBins = 40;



    horizontalTitles = {'$L_c$','$W_f$','Curl','$t_f$','$\varphi$'};
    if strcmp(ctrl.interpreter,'tex')
        horizontalTitles = strrep(horizontalTitles,'$','');
    end

    figure('units','centimeters','OuterPosition',[10 10 24 24]);

    % Histograms
    subplot(5,5,1)
    histogram(lc,linspace(0,7.6,numBins),histogramInstructions{:})
    hold on
    histogram(lc2,linspace(0,7.6,numBins),histogramInstructions{:})
    legend(horizontalTitles{1},'location','northeast','interpreter',ctrl.interpreter)
    xlim([0 2]) % Lc

    subplot(5,5,7)
    histogram(width,linspace(0,100,numBins),histogramInstructions{:})
    hold on
    histogram(width2,linspace(0,100,numBins),histogramInstructions{:})
    legend(horizontalTitles{2},'location','northeast','interpreter',ctrl.interpreter)
    xlim([0 50]) % Width
    
    subplot(5,5,13)
    histogram(wallTkn,linspace(0,60,numBins),histogramInstructions{:})
    hold on
    histogram(wallTkn2,linspace(0,60,numBins),histogramInstructions{:})
    legend(horizontalTitles{4},'location','northeast','interpreter',ctrl.interpreter)
    xlim([0 20]) % tf

    subplot(5,5,19)
    histogram(curl,linspace(0,1,numBins),histogramInstructions{:})
    hold on
    histogram(curl2,linspace(0,1,numBins),histogramInstructions{:})
    legend(horizontalTitles{3},'location','northeast','interpreter',ctrl.interpreter)
    xlim([0 2]) % Curl

    subplot(5,5,25)
    histogram(fibrillation,linspace(0,60,numBins),histogramInstructions{:})
    hold on
    histogram(fibrillation2,linspace(0,60,numBins),histogramInstructions{:})
    legend(horizontalTitles{5},'location','northeast','interpreter',ctrl.interpreter)
    xlim([0 50]) % Fibrillation

    subplot(5,5,2)
    plot(width,lc,'.',lineInstructions{:})
    hold on
    plot(width2,lc2,'.',lineInstructions{:})
    xlim([0 50]) % Width
    ylim([0 2]) % Lc

    subplot(5,5,3)
    plot(wallTkn,lc,'.',lineInstructions{:})
    hold on
    plot(wallTkn2,lc2,'.',lineInstructions{:})
    xlim([0 20]) % tf
    ylim([0 2]) % Lc
    
    subplot(5,5,4)
    plot(curl,lc,'.',lineInstructions{:})
    hold on
    plot(curl2,lc2,'.',lineInstructions{:})
    xlim([0 2]) % Curl
    ylim([0 2]) % Lc
    
    subplot(5,5,5)
    plot(fibrillation,lc,'.',lineInstructions{:})
    hold on
    plot(fibrillation2,lc2,'.',lineInstructions{:})
    xlim([0 50]) % Fibrillation
    ylim([0 2]) % Lc
    
    subplot(5,5,8)
    plot(wallTkn,width,'.',lineInstructions{:})
    hold on
    plot(wallTkn2,width2,'.',lineInstructions{:})
    xlim([0 20]) % tf
    ylim([0 50]) % Width
    
    subplot(5,5,9)
    plot(curl,width,'.',lineInstructions{:})
    hold on
    plot(curl2,width2,'.',lineInstructions{:})
    xlim([0 2]) % Curl
    ylim([0 50]) % Width
    
    subplot(5,5,10)
    plot(fibrillation,width,'.',lineInstructions{:})
    hold on
    plot(fibrillation2,width2,'.',lineInstructions{:})
    xlim([0 50]) % Fibrillation
    ylim([0 50]) % Width
    
    subplot(5,5,14)
    plot(curl,wallTkn,'.',lineInstructions{:})
    hold on
    plot(curl2,wallTkn2,'.',lineInstructions{:})
    xlim([0 2]) % Curl
    ylim([0 20]) % tf
    
    
    subplot(5,5,15)
    plot(fibrillation,wallTkn,'.',lineInstructions{:})
    hold on
    plot(fibrillation2,wallTkn2,'.',lineInstructions{:})
    xlim([0 50]) % Fibrillation
    ylim([0 20]) % tf   ylim([0 1]) % Curl
    
    subplot(5,5,20)
    plot(fibrillation,curl,'.',lineInstructions{:});
    hold on
    plot(fibrillation2,curl2,'.',lineInstructions{:});
    xlim([0 50]) % Fibrillation
    ylim([0 2]) % Curl
    
    corrIdx = [2 1 ;
               3 1 ;
               3 2 ;
               4 1 ;
               4 2 ;
               4 3 ;
               5 1 ;
               5 2 ;
               5 3 ;
               5 4 ];
    splotIdx = [6 11 12 16 17 18 21 22 23 24];
    for aLoop = 1:10
        subplot(5,5,splotIdx(aLoop))
    %     text(0.5,0.5,sprintf('%3.2f \n %3.2f',corrMatrix(corrIdx(aLoop,1),corrIdx(aLoop,2)),corrMatrix2(corrIdx(aLoop,1),corrIdx(aLoop,2))),textInstructions{:})
        text(0.5,0.5,sprintf('%3.2f',corrMatrix(corrIdx(aLoop,1),corrIdx(aLoop,2))),textInstructions{:})
        axis off
    end
    set(findobj(gcf,'type','axes'),'TickLabelInterpreter',ctrl.interpreter);



    for bLoop = 1:numel(saveExtensions)
            print(saveName,saveExtensions{bLoop})
    end
    close;
end

% print([ctrl.saveDir filesep 'overview_syn' ctrl.saveName(1:end-4)],'-dpng')
% print([ctrl.saveDir filesep 'overview_syn' ctrl.saveName(1:end-4)],'-dpng')
% print([ctrl.saveDir filesep 'overview_syn' ctrl.saveName(1:end-4)],'-dpng')
% print(['C:\Users\augus\Documents\softwareProjects\copulamodel\delmePlots\' 'overview_syn' strrep(challengeName(1:end-4),'_','')],'-dpng')
% print(['C:\Users\augus\Documents\softwareProjects\copulamodel\delmePlots\' 'overview_syn' strrep(challengeName(1:end-4),'_','')],'-dpdf')
% print(['C:\Users\augus\Documents\softwareProjects\copulamodel\delmePlots\' 'overview_syn' strrep(challengeName(1:end-4),'_','')],'-depsc')
