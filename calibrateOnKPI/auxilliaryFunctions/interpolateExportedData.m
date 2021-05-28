function [results] = interpolateExportedData(results)

lengthOfNewData = 200;
newEpsilon = [linspace(0,-5,lengthOfNewData)]';


for bLoop = 1:numel(results)
    % Force and debonding data
    if strcmp(results(bLoop).loadCase,'tension')
        newEpsilon =  abs(newEpsilon);
    else
        newEpsilon = -abs(newEpsilon);
    end
    
    if length(results(bLoop).reaFile(:,1))>2
        results(bLoop).RSreaForce = interp1(results(bLoop).reaFile(1:end-1,1),results(bLoop).reaFile(1:end-1,2),newEpsilon);  
        results(bLoop).RSreaDebo = interp1(results(bLoop).reaFile(1:end-1,1),results(bLoop).reaFile(1:end-1,4),newEpsilon);  
        results(bLoop).RSstrain = newEpsilon*0.01;
    else
        results(bLoop).RSreaForce = interp1(results(bLoop).reaFile(1:end,1),results(bLoop).reaFile(1:end,2),newEpsilon);  
        results(bLoop).RSreaDebo = interp1(results(bLoop).reaFile(1:end,1),results(bLoop).reaFile(1:end,4),newEpsilon);  
        results(bLoop).RSstrain = newEpsilon*0.01;        
    end
    
    
   % Stress and slenderness data    
   

   
   results(bLoop).density = results(bLoop).grammage/1e3/(results(bLoop).thickness*1e-6);
   
   disp(num2str(bLoop))
   [results(bLoop).peakForce,results(bLoop).peakStressIdx] = max(abs(results(bLoop).RSreaForce));
   results(bLoop).peakForce = results(bLoop).peakForce*sign(results(bLoop).RSreaForce(results(bLoop).peakStressIdx));
   results(bLoop).stress = results(bLoop).RSreaForce/results(bLoop).width/results(bLoop).thickness;
   results(bLoop).peakStress = results(bLoop).peakForce/results(bLoop).width/results(bLoop).thickness;
   results(bLoop).slenderness = 2*sqrt(3)*results(bLoop).length/results(bLoop).thickness; 
   results(bLoop).Ett =  (results(bLoop).reaFile(2:end,2)-results(bLoop).reaFile(1:end-1,2))/results(bLoop).width/results(bLoop).thickness ...
                        ./ (0.01*(results(bLoop).reaFile(2:end,1)-results(bLoop).reaFile(1:end-1,1)));
   indexOffset = 3;
   if length(results(bLoop).reaFile(:,1))>2
        results(bLoop).RSreaForceTrunc = results(bLoop).RSreaForce(1:min(results(bLoop).peakStressIdx+indexOffset,length(results(bLoop).RSreaForce)));  
        results(bLoop).RSreaDeboTrunc = results(bLoop).RSreaDebo(1:min(results(bLoop).peakStressIdx+indexOffset,length(results(bLoop).RSreaForce)))./results(bLoop).width; 
        results(bLoop).RSstrainTrunc = results(bLoop).RSstrain(1:min(results(bLoop).peakStressIdx+indexOffset,length(results(bLoop).RSreaForce)));
        results(bLoop).stressTrunc = results(bLoop).stress(1:min(results(bLoop).peakStressIdx+indexOffset,length(results(bLoop).RSreaForce)));

      results(bLoop).EMax = mean(results(bLoop).Ett(2));
   results(bLoop).drop = results(bLoop).reaFile(end,2)/max(results(bLoop).reaFile(:,2));        
        
   else
        results(bLoop).RSreaForceTrunc = results(bLoop).RSreaForce(1:min(results(bLoop).peakStressIdx+indexOffset,length(results(bLoop).RSreaForce))); 
        results(bLoop).RSreaDeboTrunc = results(bLoop).RSreaDebo(1:min(results(bLoop).peakStressIdx+indexOffset,length(results(bLoop).RSreaForce)))./results(bLoop).width; 
        results(bLoop).RSstrainTrunc = results(bLoop).RSstrain(1:min(results(bLoop).peakStressIdx+indexOffset,length(results(bLoop).RSreaForce)));%0.01*newEpsilon((1:results(bLoop).peakStressIdx+indexOffset));  
        results(bLoop).stressTrunc = results(bLoop).stress(1:min(results(bLoop).peakStressIdx+indexOffset,length(results(bLoop).RSreaForce)));
        
           results(bLoop).EMax = nan;
   results(bLoop).drop = nan;        
   
   end                
   
            
end