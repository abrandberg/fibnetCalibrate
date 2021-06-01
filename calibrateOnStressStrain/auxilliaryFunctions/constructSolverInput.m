function constructSolverInput(x,fibnet,ctrl,paramName)


fid  = fopen(ctrl.apdlHeader,'r');
fHeader = fread(fid,'*char')';
fclose(fid);


fid  = fopen(ctrl.apdlInput,'r');
fInput = fread(fid,'*char')';
fclose(fid);



fieldsInFibnetStruct = fields(fibnet);
valuesPerFibnetField = struct2cell(fibnet);

writeStringFix{1} = '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!';
writeStringFix{2} = '! START OF FIXED BUT PARAMETRIZED VALUES';
for aLoop = 1:numel(fieldsInFibnetStruct)
    if isnumeric(valuesPerFibnetField{aLoop})
        writeStringFix{aLoop+2} = [fieldsInFibnetStruct{aLoop} ' = ' num2str(valuesPerFibnetField{aLoop})];
    else
        writeStringFix{aLoop+2} = [fieldsInFibnetStruct{aLoop} ' = ' char(strjoin(valuesPerFibnetField{aLoop},''))];
    end
 
end
writeStringFix{end+1} = '! END OF FIXED BUT PARAMETRIZED VALUES';
writeStringFix{end+1} = '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!';



assert(numel(paramName) == numel(x))
writeString{1} = '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!';
writeString{2} = '! START OF FITTING PARAMETERS';
for aLoop = 1:numel(paramName)
    
    writeString{aLoop+2} = [paramName{aLoop} ' = ' num2str(x(aLoop))];
 
end
writeString{end+1} = '! END OF FITTING PARAMETERS';
writeString{end+1} = '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!';


fileID = fopen([ctrl.runDir filesep ctrl.simName '.dat'],'w');
fprintf(fileID,'%s\n',fHeader,writeStringFix{:},writeString{:},fInput);
fclose(fileID);
