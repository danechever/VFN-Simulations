function couplingEffMap = generateCouplingMap( fiberMode, E_PSF, totalPower, cropsize)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[rowsMode,colsMode] = size(fiberMode); 

croprows = rowsMode/2-round(cropsize/2)+1:rowsMode/2+round(cropsize/2);
cropcols = colsMode/2-round(cropsize/2)+1:colsMode/2+round(cropsize/2);

fibermode_cropped = fiberMode(croprows,cropcols);

couplingEffMap = abs(conv2(E_PSF,fibermode_cropped,'same')).^2/totalPower;

end

