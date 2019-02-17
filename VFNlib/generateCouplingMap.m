function couplingEffMap = generateCouplingMap( fiberMode, E_PSF, totalPower, cropsize)
%couplingEffMap = generateCouplingMap( fiberMode, E_PSF, totalPower, cropsize)
%   Generates a 2D map of the coupling efficiency versus source position.
%   
%   Inputs:
%       fiberMode: 2D array with the fiber mode
%       E_PSF: 2D array withoint spread function electric field 
%       totalPower: Total power/energy in the PSF; i.e. sum(abs(E_PSF(:)).^2))
%       cropsize: Size to crop the fiberMode to for convolution. Smaller is
%           faster but less accurate. Several times lambdaOverD is recommended.
%   Outputs:
%       couplingEffMap: 2D array with coupling efficiency versus source
%           position. Same size as E_PSF. 

    [rowsMode,colsMode] = size(fiberMode); 

    croprows = rowsMode/2-round(cropsize/2)+1:rowsMode/2+round(cropsize/2);
    cropcols = colsMode/2-round(cropsize/2)+1:colsMode/2+round(cropsize/2);

    fibermode_cropped = fiberMode(croprows,cropcols);

    couplingEffMap = abs(conv2(E_PSF,fibermode_cropped,'same')).^2/totalPower;

end

