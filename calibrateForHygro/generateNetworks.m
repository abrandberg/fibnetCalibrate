function PaperDir = generateNetworks(netgen,ctrl)

numGrammage = numel(netgen.grammage);
numThickness = numel(netgen.thickness);
numAngle = numel(netgen.angleStd);

    
maxTests = max([numGrammage numThickness numAngle]);
minTests = min([numGrammage numThickness numAngle]);

if numGrammage < maxTests && numGrammage == 1
    netgen.grammage = netgen.grammage*ones(maxTests,1);
end

if numThickness < maxTests && numThickness == 1
    netgen.thickness = netgen.thickness*ones(maxTests,1);
end



% Generate networks
PAUS = ' ';
PAUS(2:ctrl.parallelNetGetRuns)='&';
% 
fid  = fopen('ModelingData.txt','r');
f_in=fread(fid,'*char')';
fclose(fid);
f_in = strrep(f_in,'10.0d-3					NetworkLength', [num2str(netgen.length) 'd-3					NetworkLength']);
f_in = strrep(f_in,'10.0d-3					NetworkWidth',  [num2str(netgen.width) 'd-3					NetworkWidth']);
f_in = strrep(f_in,'outID.txt			        MeasurementFile',      [ctrl.pulpDataFile  '			        MeasurementFile']);
f_in = strrep(f_in,'9.0				    	InterfaceAngle', [num2str(netgen.interfaceAngle) '.0				    	InterfaceAngle']);



for i = 1:length(netgen.angleStd)

    if (netgen.angleStd(i)==0)
        PaperDir{i}  = 'P0.0';
    else
        PaperDir{i}  = sprintf(['P%' num2str(fix(log10(netgen.angleStd(i))+2)) '.1f'],netgen.angleStd(i));
    end
%     disp(PaperDir{i})
    f = strrep(f_in,'100.0d-3 				Grammage',      [num2str(netgen.grammage(i))  'd-3 				Grammage']);
    f = strrep(f,'80.0d-6					PaperThk',      [num2str(netgen.thickness(i))  'd-6					PaperThk']);
    f = strrep(f,'0.0d0					AngleStd',[num2str(netgen.angleStd(i)) '					AngleStd']);

    FName = sprintf(['file1_L%' num2str(fix(log10(netgen.length)+3)) '.1f_W%' num2str(fix(log10(netgen.width)+3)) '.1f_g%' num2str(fix(log10(netgen.grammage(i))+3)) '.1f'],netgen.length,netgen.width,netgen.grammage(i));
    
    if exist([ctrl.runDir PaperDir{i} '\' FName '.dat'],'file')~=0 ,continue,end %&& ctrl.forceRewrite == 0
    
%     if ctrl.forceRewrite == 1
%         delete([ctrl.runDir PaperDir{i} filesep '*'])
%     end
    
    system(['md ' ctrl.runDir PaperDir{i}]);
    system(['copy MyPacking.exe ' ctrl.runDir PaperDir{i} '\MyPacking.exe']);
    system(['copy ' ctrl.pulpDataFile ' ' ctrl.runDir PaperDir{i} '\' ctrl.pulpDataFile]);
    fid  = fopen([ctrl.runDir PaperDir{i} '\ModelingData.txt'],'w');
    fprintf(fid,'%s',f);
    fclose(fid);
    system(['cd ' ctrl.runDir PaperDir{i} ' && MyPacking.exe' PAUS(1+mod(i,ctrl.parallelNetGetRuns))]);

end

% Copy auxilliary files
copyfile('FDCurveBE.mac',[ctrl.runDir filesep 'FDCurveBE.mac'])

copyfile(['auxilliaryFunctions\femOutputCharacterization\exportNetworkData.mac'],[ctrl.runDir filesep 'exportNetworkData.mac'])
copyfile(['auxilliaryFunctions\femOutputCharacterization\ABdstress_energy.mac'],[ctrl.runDir filesep 'ABdstress_energy.mac'])
copyfile(['auxilliaryFunctions\femOutputCharacterization\freeSpan4.mac'],[ctrl.runDir filesep 'freeSpan4.mac'])




