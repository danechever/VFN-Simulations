function mode = generateSMFmode( fiber_props, lambda, Fnum, lambdaFnum_samps, coords)
%mode = generateSMFmode( fiber_props, lambda, Fnum, lambdaFnum_samps, coords)
%   Generates the fundamental mode of a single mode fiber. 
%
%   Inputs:
%       fiber_props: Structure with fiber properties 
%           fiber_props.n_core - core index 
%           fiber_props.n_clad - cladding index 
%           fiber_props.NA - numerical aperture (can be defined instead of
%                            n_core and n_clad).
%           fiber_props.core_rad - core radius  
%
%       lambda: wavelength (meters)
%       Fnum: Focal ratio (f/D)
%       lambdaFnum_samps: lambda times the focal ratio in samples 
%       coords: coordinates structure 
%
%   Outputs:
%       mode: 2D array with normalized mode field 

    if(~isfield(fiber_props,'NA'))
        fiber_props.NA = sqrt(fiber_props.n_core^2-fiber_props.n_clad^2); 
    end

    % Get fiber mode
    if(strcmpi('bessel',fiber_props.type))
        mode = generateSMFmode_bessel( fiber_props.NA, fiber_props.core_rad, lambda, lambda*Fnum/lambdaFnum_samps, coords );
    elseif(strcmpi('gaussian',fiber_props.type))
        MFD = getMFD(fiber_props,lambda);% meters
        MFD_lamoverd = MFD/(lambda*Fnum);% units of lambda/D
        mode = generateSMFmode_gaussian(MFD_lamoverd*lambdaFnum_samps,coords);
    end

end