function Iy = calculateIy(OneFreeSpan)
% This part uses the actual fiber properties to calculate the
% second area moment of inertia. IF clause to account for the two
% types of cross section (hollow and solid rectangle) present in
% the network.
if OneFreeSpan(end,3) == 1 
 Iy = (OneFreeSpan(end,4)*OneFreeSpan(end,5)^3)/12;
else
 Iy = (OneFreeSpan(end,4)*OneFreeSpan(end,5)^3)/12 - ...
      (OneFreeSpan(end,4)-OneFreeSpan(end,6)*(OneFreeSpan(end,5)-OneFreeSpan(end,6))^3)/12;
end