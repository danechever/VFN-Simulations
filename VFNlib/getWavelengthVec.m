function lams = getWavelengthVec( lambda0, fracBW, numWavelengths )
%lams = getWavelengthVec( lambda0, fracBW, numWavelengths )
%   Generates wavelength vector 
%
%   Inputs:
%       lambda0: central wavelength (meters)
%       fracBW: fractional bandwidth (\Delta\lambda/\lambda)
%       numWavelengths: number of wavelengths
%   Outputs:
%       lams: wavelengths (meters)
    
    if numWavelengths == 1
        lams = lambda0;
    else
        lams = lambda0*linspace(1-fracBW/2,1+fracBW/2,numWavelengths);
    end

end