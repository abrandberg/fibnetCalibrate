function [PaperDir,FName] = generateNetworks(netgen,ctrl)





% Generate networks
% PAUS = ' ';
% PAUS(2:ctrl.parallelNetGetRuns)='&';
% 
fid  = fopen(ctrl.modelingDataPointer,'r');
f_in=fread(fid,'*char')';
fclose(fid);
f_in = strrep(f_in,'10.0d-3					NetworkLength', [num2str(netgen.length) 'd-3					NetworkLength']);
f_in = strrep(f_in,'10.0d-3					NetworkWidth',  [num2str(netgen.width) 'd-3					NetworkWidth']);
f_in = strrep(f_in,'outID.txt			        MeasurementFile',      [flip(strtok(flip(ctrl.pulpDataFile),filesep))  '			        MeasurementFile']);
f_in = strrep(f_in,'9.0				    	InterfaceAngle', [num2str(netgen.interfaceAngle) '.0				    	InterfaceAngle']);



for i = 1

    if (netgen.angleStd(i)==0)
    PaperDir{i}  = 'P0.0';
    else
    PaperDir{i}  = sprintf(['P%' num2str(fix(log10(netgen.angleStd(i))+2)) '.1f'],netgen.angleStd(i));
    end

    f = strrep(f_in,'100.0d-3 				Grammage',      [num2str(netgen.grammage(i))  'd-3 				Grammage']);
    f = strrep(f,'80.0d-6					PaperThk',      [num2str(netgen.thickness(i))  'd-6					PaperThk']);
    f = strrep(f,'0.0d0					AngleStd',[num2str(netgen.angleStd(i)) '					AngleStd']);

    FName = sprintf(['file1_L%' num2str(fix(log10(netgen.length)+3)) '.1f_W%' num2str(fix(log10(netgen.width)+3)) '.1f_g%' num2str(fix(log10(netgen.grammage(i))+3)) '.1f'],netgen.length,netgen.width,netgen.grammage(i));

    if exist([ctrl.runDir PaperDir{i} filesep FName '.dat'],'file')~=0 ,continue,end 

    mkdir([ctrl.runDir PaperDir{i}])
    copyfile(ctrl.myPackingPointer,[ctrl.runDir PaperDir{i} filesep 'MyPacking.exe'])
    copyfile(ctrl.pulpDataFile,[ctrl.runDir PaperDir{i} filesep flip(strtok(flip(ctrl.pulpDataFile),filesep))])

    fid  = fopen([ctrl.runDir PaperDir{i} '\ModelingData.txt'],'w');
    fprintf(fid,'%s',f);
    fclose(fid);

    cd([ctrl.runDir PaperDir{i}])
    system('MyPacking.exe')
    cd(ctrl.workDir)


end




