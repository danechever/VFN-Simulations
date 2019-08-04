function coupling_eff_map_cube = generateCouplingMap_polychromatic( Epup, fiber_props, lambda0, Fnum, lambdas, totalPower, lambdaOverD, cropsize, coords, fibmodes)
%coupling_eff_map_cube = generateCouplingMap_polychromatic( Epup, fiber_props, lambda0, Fnum, lambdas, totalPower, lambdaOverD, cropsize, coords, fibmodes)
%   Generates a 2D maps of the coupling efficiency versus source position.
%   Returns cube of maps with dimensions (NxNxnumWavelengths)
%   
%   Inputs:
%       Epup: the entrance pupil field (NxNxnumWavelengths)
%       fiber_props: structure of fiber properties 
%       lambda0: central wavelength (meters)
%       Fnum: focal ratio at the fiber (f/Dbeam)
%       lambdas: wavelengths (meters)
%       totalPower: Total power/energy in the PSF at lambda0; i.e. sum(abs(E_PSF(:)).^2))
%       lambdaOverD: samples per lambda/D in the image plane
%       cropsize: Size to crop the fiberMode to for convolution. Smaller is
%           faster but less accurate. Several times lambdaOverD is recommended.
%       coords: coordinates structure 
%       fibmodes: [OPTIONAL] Fiber modes - reduces runtime when looping on this func
%                  Note: fibmodes should contain the mode at each wavelength
%                   
%   Outputs:
%       coupling_eff_map_cube: Cube where each slice is 2D array 
%           with coupling efficiency versus source position. The slices
%           corresponds to lambdas.

    myfft2 = @(x) fftshift(fft2(fftshift(x)));

    coupling_eff_map_cube = zeros(coords.N,coords.N,numel(lambdas));

    ch = 1;% wavelength channel number
    for lam = lambdas
        
        lam_frac = lam/lambda0;% fractional wavelength (\lambda/\lambda_0)
        
        % Compute PSF field
        PSFv = myfft2(Epup(:,:,ch)); 
        
        if nargin >= 10
            %-- Fiber modes provided via fibmodes; pick out mode at lam
            fibermode = fibmodes(:,:,ch);
        else
            %-- Fiber modes not provided; generate them
            % Generate the fiber mode for this wavelength with proper scaling
            fibermode = generateSMFmode( fiber_props, lam, Fnum, lambdaOverD, coords);
        end
        
        % Compute the monochromatic coupling map (spatial units of lambda/D)
        map = generateCouplingMap( fibermode, PSFv, totalPower, cropsize, coords);
    
        % Scale to units of lambda0/D
        map_rescaled = interp2(coords.X,coords.Y,map,coords.X/lam_frac,coords.Y/lam_frac,'linear',0);
        
        % Add it the the cube
        coupling_eff_map_cube(:,:,ch) = map_rescaled;
        
        ch = ch + 1;% increment wavelength channel number
    end 
end

