function coeffs = determineRHtoMCRatio(ctrl,dataToImport)


if strcmp(dataToImport,'monopolX')
    monopolXData = [20.02 3.50
                    32.82 7.39
                    49.86 11.10
                    79.56 16.79
                    49.73 14.05
                    32.85 10.76
                    19.86 3.84
                    49.85 11.79
                    79.74 15.84];
     
     dataToAnalyze = monopolXData;
elseif strcmp(dataToImport,'handsheetLowDens')
    handsheetData =    [50	2.605	10.56876061
                    33	2.57	9.083191851
                    20	2.524	7.130730051
                    33	2.557	8.531409168
                    50	2.597	10.22920204
                    66	2.643	12.18166384
                    80	2.711	15.06791171
                    66	2.674	13.49745331
                    50	2.626	11.46010187];
    handsheetData(:,2) = [];   
    
    dataToAnalyze = handsheetData;
    
elseif strcmp(dataToImport,'handsheetMidDens')
    handSheetDataMid = [50	2.62	 9.761206535
                        33	2.605	 9.132802681
                        20	2.572	 7.750314202
                        33	2.601	 8.96522832
                        50	2.637	10.47339757
                        66	2.689	12.65186426
                        80	2.756	15.45873481
                        66	2.718	13.86677838
                        50	2.648	10.93422706];
    handSheetDataMid(:,2) = [];
    dataToAnalyze = handSheetDataMid;
elseif strcmp(dataToImport,'handsheetHighDens')
    handSheetDataHigh = [50     2.625	9.329446064
                         33     2.598	8.204914619
                         20     2.565	6.830487297
                         33     2.596	8.121615993
                         50     2.629	9.496043315
                         66     2.672	11.28696377
                         80     2.736	13.95251978
                         66     2.69	12.0366514
                         50     2.666	11.03706789
                         33     2.611	8.746355685
                         20     2.571	7.080383174
                         33     2.608	8.621407747
                         50     2.644	10.12078301
                         80     2.737	13.9941691
                         20     2.573	7.163681799
                         80     2.739	14.07746772
                         50     2.665	10.99541858];
   handSheetDataHigh(:,2) = [];  
   dataToAnalyze = handSheetDataHigh;              
end

% Moisture iso-therm


   




% dataToAnalyze = handsheetData;
coeffs = polyfit(dataToAnalyze(:,1),dataToAnalyze(:,2),1);



if 0%1
    figure('color','w','units','centimeters','OuterPosition',[10 10 16 16]);
    plot(dataToAnalyze(:,1),dataToAnalyze(:,2),'-o','color',[0.5 0.5 0.5],ctrl.lineInstructions{:})
    hold on
    plot(linspace(20,80,60),polyval(coeffs,linspace(20,80,60)),'-k',ctrl.lineInstructions{:},'linewidth',2)
    xlabel('Relative humidity [\%]','interpreter',ctrl.interpreter)
    ylabel('Moisture content [\%]','interpreter',ctrl.interpreter)
    xlim([0 100])
    legend('Data (TH)',['Fit ' sprintf('mc(RH) = %3.2f$\\cdot$RH + %3.2f ',coeffs(1),coeffs(2) )],'location','northwest','interpreter',ctrl.interpreter)
    set(gca,'TickLabelInterpreter',ctrl.interpreter,'fontsize',14)
    print('cdl07_2','-dpng','-r0')
end