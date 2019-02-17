function mode = generateSMFmode( n_core, n_clad, core_rad, lambda, dx, coords )
%mode = generateSMFmode( n_core, n_clad, core_rad, lambda, dx, coords )
%   Generates the fundamental mode of a single mode fiber. 
%
%   Inputs:
%       n_core: the refractive index of the fiber core. 
%       n_clad: the refractive index of the fiber cladding. 
%       core_rad: radius of the fiber core (meters) 
%       lambda: wavelength (meters)
%       dx: sample spacing in the image plane (meters) 
%       coords: coordinates structure 
%   Outputs:
%       mode: 2D array with normalized mode field 

    % Normalized frequency 
    V = 2*pi*core_rad/lambda*sqrt(n_core^2-n_clad^2); 
    
    % mode parameters (see any fiber optics textbook)
    u = fzero(@(x) besselj(0,x)./(x.*besselj(1,x)) - besselk(0,sqrt(V^2 - x.^2))./(sqrt(V^2 - x.^2).*besselk(1,sqrt(V^2 - x.^2))),1);
    w = sqrt(V^2 - u^2); 
    
    % mode field distributions inside and outside core
    mode_in_core = besselj(0,u*coords.RHO*dx/core_rad)/besselj(0,u);
    mode_outside_core = besselk(0,w*coords.RHO*dx/core_rad)/besselk(0,w);

    % Combine inner and outer fields into a single array 
    mode = mode_in_core;
    mode(coords.RHO*dx>core_rad) = mode_outside_core(coords.RHO*dx>core_rad);
    
    % Return peak normalized mode field
    mode = mode/sqrt(sum(mode(:).^2));

end