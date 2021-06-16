function executeFibnetRuns(PaperDir,netgen,fibnet,simulation,ctrl,materialProps)
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


PAUS = ' ';
PAUS(2:ctrl.parallelFibNetRuns)='&';
counter = 0;


for aLoop = 1:numel(simulation)

    fid  = fopen(['HygVsEMod-' simulation(aLoop).type '.mac'],'r');
    try
    f_in=fread(fid,'*char')';
    catch
        disp('ok')
    end
    fclose(fid);

    

    f_in = strrep(f_in,'len = 5000',[ 'len = ' num2str(fibnet.length)]);
    f_in = strrep(f_in,'wid = 5000',[ 'wid = ' num2str(fibnet.width)]);

    f_in = strrep(f_in,'gxy,i,cEx/10','gxy,i,12e3');
%     f_in = strrep(f_in,'ez,i,cEx/10','ez,i,8.5e3');
    f_in = strrep(f_in,'ez,i,cEx/10',['ez,i,' num2str(materialProps.EZ)]);
    f_in = strrep(f_in,'alpx,i,0.06/100','alpx,i,0.02/100');        
    f_in = strrep(f_in,'alpz,i,0.5/100','alpz,i,1.0*0.40/100');         %0.35 for 2x2
%     f_in = strrep(f_in,'ex,i,cEx','ex,i,19e3');
    f_in = strrep(f_in,'ex,i,cEx',['ex,i,' num2str(materialProps.EX)]);
    
    f_in = strrep(f_in,'E_ct = 20*10032*kof',['E_ct = ' num2str(materialProps.Kct)]);
    f_in = strrep(f_in,'usr6,meth,1,fract,partype,0.90',['usr6,meth,1,fract,partype,' num2str(materialProps.searchRadiusMultiplier)]);
    
    
    
    for bLoop = 1:length(netgen.angleStd)
        
        FName = sprintf(['file1_L%' num2str(fix(log10(netgen.length)+3)) '.1f_W%' num2str(fix(log10(netgen.width)+3)) '.1f_g%' num2str(fix(log10(netgen.grammage(bLoop))+3)) '.1f'],netgen.length,netgen.width,netgen.grammage(bLoop));
        f = strrep(f_in,'file1_L10.0_W10.0_g100.0',FName);
        
        Project = [simulation(aLoop).type '_AngStd' num2str(netgen.angleStd(bLoop)) simulation(aLoop).dir];
        
               
        if exist([ctrl.runDir Project '.rea'],'file')~=0 && ctrl.forceRewrite == 0,continue,end
        
        
        disp(Project)

        f = strrep(f,'*get,np,runst,,rspeed,nproc',['np = ' num2str(ctrl.fibnetNP)]);
        f = strrep(f,'Paper/',[ PaperDir{bLoop} '/']);
        
        if strcmp(simulation(aLoop).dir,'_x')
            writeDir = 0;
        elseif strcmp(simulation(aLoop).dir,'_y')
            writeDir = 1;
        else
            writeDir = 0;
        end
        
        f = strrep(f,'Orient = 0 ',['Orient = ' num2str(writeDir)]);

        fid  = fopen([ctrl.runDir Project '.dat'],'w');
        fprintf(fid,'%s',f);
        fclose(fid);
        
        counter = counter+1;
        while (exist([ctrl.runDir PaperDir{bLoop} filesep FName '.dat'],'file')~=2)
        end
        
%         if ( aLoop == numel(simulation) ) & ( bLoop == numel(netgen.angleStd) )
%             system(['SET KMP_STACKSIZE=4096k & "C:\Program Files\ANSYS Inc\v150\ANSYS\bin\winx64\ansys150.exe"  -p aa_t_a -np 2 -dir "' ctrl.runDir '" -j "' Project '" -s read -l en-us -b -i "' ctrl.runDir Project '.dat" -o "' ctrl.runDir Project '.out" -custom "' ctrl.executable '" ' PAUS(1)]);
%         else
            if exist([ctrl.runDir filesep Project '.rea'],'file') == 2
                delete([ctrl.runDir filesep Project '.rea'])
            end
            system(['SET KMP_STACKSIZE=4096k & "C:\Program Files\ANSYS Inc\v150\ANSYS\bin\winx64\ansys150.exe"  -p aa_t_a -np 2 -dir "' ctrl.runDir '" -j "' Project '" -s read -l en-us -b -i "' ctrl.runDir Project '.dat" -o "' ctrl.runDir Project '.out" -custom "' ctrl.executable '" ' PAUS(rem(counter,ctrl.parallelFibNetRuns)+1) ' ']);
%         end
        
        %system(['SET KMP_STACKSIZE=4096k & "C:\Program Files\ANSYS Inc\v150\ANSYS\bin\winx64\ansys150.exe"  -p aa_t_a -np 1 -dir "' ctrl.runDir '" -j "' Project '" -s read -l en-us -b -i "' ctrl.runDir Project '.dat" -o "' ctrl.runDir Project '.out" -custom "' ctrl.executable '" ' PAUS(rem(counter,ctrl.parallelFibNetRuns)+1)]);

    end

end
% save(['HygVsE.mat']);
