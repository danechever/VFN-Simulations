function iPSF = getallPSF( Epup, lambda0, lambdas, normI, coords )
%iPSF_BB = getPSF( Epup, lambda0, lambdas, normI, coords )
%   Generates the point spread function (square magnitude of the field) 
%
%   Inputs:
%       Epup: the entrance pupil field (NxNxnumWavelengths)
%       lambda0: central wavelength (meters)
%       lambdas: wavelengths (meters)
%       normI: irradiance normalization factor
%       coords: coordinates structure 
%   Outputs:
%       iPSF: 2D array with point spread function (normalized irradiance) 

    myfft2 = @(x) fftshift(fft2(fftshift(x)));
    
    iPSF = zeros(coords.N,coords.N,length(lambdas)); 
    
    ch = 1;% channel index
    for lam = lambdas 
        PSF = myfft2(Epup(:,:,ch)); % PSF (complex field)
        iPSF_lam = abs(PSF).^2/normI;
        lam_frac = lam/lambda0;
        iPSF(:,:,ch) = (1/lam_frac)^2*interp2(coords.X,coords.Y,iPSF_lam,coords.X/lam_frac,coords.Y/lam_frac,'linear',0);
        ch = ch + 1;
    end

end