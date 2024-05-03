function mask = generatePhaseStripeMask( phzStep, clockAng_deg, period_samps, oddOrEven, coords, offsets_samps )
%mask = generatePhaseStripeMask( phzStep, clockAng, coords, offsets )
%   Generates array with representation of the stripe mask 
%
%   Inputs: 
%       phzStep - the phase at each wavelength sample [nWvls x 1]
%       clockAng_deg - The clocking angle in degrees 
%       period - The period of the stripes in units of pupil diameters 
%       oddOrEven - 'odd' makes an odd phase function, 'even' makes an even
%       coords - The coordinate system structure 
%       offsets - x and y offsets to the center boundary 
%       
%   Outputs:
%       mask - 2D array with phase mask 
%
%

    mask = zeros(coords.N,coords.N,numel(phzStep));
    
    Xoffset = coords.X-offsets_samps(1);
    Yoffset = coords.Y-offsets_samps(2);
    
    [THETA,RHO] = cart2pol(Xoffset,Yoffset);
    
    clockAng_rad = deg2rad(clockAng_deg);
    
    Xp = RHO.*cos(THETA-clockAng_rad);
    % Yp = RHO.*sin(THETA-clockAng_rad);

    for index = 1:numel(phzStep)
        switch lower(oddOrEven)
            case 'odd'
                phz = phzStep(index)*sign(sin(2*pi*Xp/period_samps))/2;
            case 'even'
                phz = phzStep(index)*sign(cos(2*pi*Xp/period_samps))/2;
        end
        % phz(phz==0) = phzStep(index)/2;
        phz(:,coords.N/2+1) = 0;
        mask(:,:,index) = exp(1i*phz);
    end

end

