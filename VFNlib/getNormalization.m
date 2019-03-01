function [normI, normP] = getNormalization( PUPIL_SUPPORT )
%[normI, normP] = getNormalization( PUPIL_SUPPORT )
%   Generates the normalization factors for irradiance and total power 
%
%   Inputs:
%       PUPIL_SUPPORT: support of the entrance pupil
%   Outputs:
%       normI: Normalization factor for irradiance
%       normI: Normalization factor for total power

    myfft2 = @(x) fftshift(fft2(fftshift(x)));
    PSF = myfft2(PUPIL_SUPPORT); % PSF with a flat wavefront 
    iPSF = abs(PSF).^2;
    normI = max(iPSF(:));% Normalization for PSF plots
    normP = sum(iPSF(:));% Normalization for coupling fractions

end