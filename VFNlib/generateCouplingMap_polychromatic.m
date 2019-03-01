function coupling_eff_map_cube = generateCouplingMap_polychromatic( Epup, fiber_props, lambda0, lambdas, totalPower, lambdaOverD, cropsize, coords)
%coupling_eff_map_cube = generateCouplingMap_polychromatic( Epup, fiber_props, lambda0, lambdas, totalPower, lambdaOverD, cropsize, coords)
%   Generates a 2D maps of the coupling efficiency versus source position.
%   Returns cube of maps with dimensions (NxNxnumWavelengths)
%   
%   Inputs:
%       Epup: the entrance pupil field (NxNxnumWavelengths)
%       fiber_props: structure of fiber properties 
%       lambda0: central wavelength (meters)
%       lambdas: wavelengths (meters)
%       totalPower: Total power/energy in the PSF at lambda0; i.e. sum(abs(E_PSF(:)).^2))
%       lambdaOverD: samples per lambda/D in the image plane
%       cropsize: Size to crop the fiberMode to for convolution. Smaller is
%           faster but less accurate. Several times lambdaOverD is recommended.
%       coords: coordinates structure 
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
        
        % Get fiber mode
        if(strcmpi('bessel',fiber_props.type))
            fibermode = generateSMFmode( fiber_props.n_core, fiber_props.n_clad, fiber_props.core_rad, lam, lam*fiber_props.Fnum/lambdaOverD, coords );
        elseif(strcmpi('gaussian',fiber_props.type))
            fibermode = generateSMFmode_gaussian(fiber_props.fiberDiam*lambdaOverD/lam_frac,coords);
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

