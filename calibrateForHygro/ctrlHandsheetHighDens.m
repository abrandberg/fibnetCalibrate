function ctrlHandsheetHighDens(ctrl,netgen,fibnet,simulation,materialProps)

coeffs = determineRHtoMCRatio(ctrl,'handsheetHighDens');
experimentalData = [4706.746    0.14775    
                    4057.886    0.16732
                    4559.684    0.14772
                    2703.32     0.26687
                    2764.12     0.23362
                    2609.578    0.25764    ]; % DMA data, taken from beaten
experimentalData(:,1) = experimentalData(:,1)*1e6;
experimentalData(:,2) = experimentalData(:,2)./(diff(polyval(coeffs,[33 66])));
experimentalSTD = [ 380.872     0.00807
                    731.104     0.00677
                    372.567     0.00710
                     39.097     0.01919
                     76.216     0.00616
                     52.551     0.00242]; % DMA data, taken from
experimentalSTD(:,1) = experimentalSTD(:,1)*1e6;
experimentalSTD(:,2) = experimentalSTD(:,2)./(diff(polyval(coeffs,[33 66])));

% Generate networks
PaperDir = generateNetworks(netgen,ctrl);


% fibnet runs
executeFibnetRuns(PaperDir,netgen,fibnet,simulation,ctrl,materialProps)


% Import results
[resultsCtrl,~] = importFibnetOutputs(PaperDir,netgen,fibnet,simulation,ctrl);

% Plot the results
selE = strcmp({resultsCtrl.type},'E');
selH = strcmp({resultsCtrl.type},'Hyg');

linFitExp = polyfit(experimentalData(:,1).*1e-9,experimentalData(:,2),1)
linFitModel = polyfit([resultsCtrl(selE).data].*1e-9,[resultsCtrl(selH).data],1)

fitXvec = [min([min(experimentalData(:,1).*1e-9),min([resultsCtrl(selE).data].*1e-9)]) max([max(experimentalData(:,1).*1e-9),max([resultsCtrl(selE).data].*1e-9)]) ];


% figure('color','w','units','centimeters','OuterPosition',[10 10 16 16]);

plot(fitXvec, polyval(linFitExp,fitXvec),'color',ctrl.colors(1,:))
hold on   
plot(fitXvec, polyval(linFitModel,fitXvec),'color',ctrl.colors(2,:))

% plot(experimentalData(:,1).*1e-9,experimentalData(:,2),'o','color',ctrl.colors(1,:),ctrl.lineInstructions{:})
errorbar(experimentalData(:,1).*1e-9,experimentalData(:,2), ...
         experimentalSTD(:,2), ...
         experimentalSTD(:,2), ...
         experimentalSTD(:,1).*1e-9, ...
         experimentalSTD(:,1).*1e-9, ...
        'o','color',ctrl.colors(1,:),'markersize',8)
hold on   
plot([resultsCtrl(selE).data].*1e-9,[resultsCtrl(selH).data],'s','color',ctrl.colors(2,:),ctrl.lineInstructions{:});
xlabel('Sheet modulus [GPa]','interpreter',ctrl.interpreter)
ylabel('Hygroexpansion $\beta$ [\%/\%mc]','interpreter',ctrl.interpreter)
xlim([0.8 5.0])
legend('Data (TH)','Model','location','northeast','interpreter',ctrl.interpreter)
set(gca,'TickLabelInterpreter',ctrl.interpreter,'fontsize',14)
print('cdl07_3B','-dpng','-r0')
pause(0.5)