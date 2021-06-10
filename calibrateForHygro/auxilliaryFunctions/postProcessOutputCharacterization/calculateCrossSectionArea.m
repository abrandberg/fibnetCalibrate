function Area = calculateCrossSectionArea(OneFreeSpan)
% This part uses the actual fiber properties to calculate the
% area. IF clause to account for the two
% types of cross section (hollow and solid rectangle) present in
% the network.
if OneFreeSpan(end,3) == 1 
     Area = OneFreeSpan(end,4)*OneFreeSpan(end,5);
else
     Area = OneFreeSpan(end,4)*OneFreeSpan(end,5) - ...
            OneFreeSpan(end,4)-OneFreeSpan(end,6)*(OneFreeSpan(end,5)-OneFreeSpan(end,6));
end