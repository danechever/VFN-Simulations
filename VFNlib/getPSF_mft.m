function [iPSF, PSF] = getPSF_mft( Epup, lambdas, foc, coordsPP, coordsFP )
%iPSF_BB = getPSF( Epup, lambdas, foc, coordsPP, coordsFP )
%   Generates the point spread function (E-field itself) 
%   - uses falco's propcustom_mft_PtoF() function 
%   - Returns the wavelength-averaged PSF (Broadband intensity)
%   - Also returns a cube with the E-field at each wavelength
%
%   Inputs:
%       Epup: the entrance pupil field (NxNxnumWavelengths)
%       lambdas: wavelengths (meters)
%       foc: focal length of focusing optic in meters
%       coordsPP: Pupil-plane coordinates structure
%               - Must contain: N (num pix), dx (pix resolution in m)
%       coordsFP: Focal-plane coordinates structure
%               - Must contain: N (num pix), dx (pix resolution in m)
%   Outputs:
%       iPSF: 2D array with point spread function (normalized irradiance) 
%       PSF:  cube with complex-valued electric field in focal plane

    numWavelengths = length(lambdas);
    PSF = nan(coordsFP.N,coordsFP.N,numWavelengths);
    
    for ch = 1:numWavelengths
        % Get PSF at each wavelength (channel)
        PSF(:,:,ch) = propcustom_mft_PtoF(Epup(:,:,ch),foc, lambdas(ch), coordsPP.dx, coordsFP.dx, coordsFP.N, coordsFP.dx, coordsFP.N);
    end

    iPSF = mean(abs(PSF).^2,3);
    iPSF = iPSF/max(iPSF(:));
end