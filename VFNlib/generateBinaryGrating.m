function grating = generateBinaryGrating( gratingAmp, period_samp, clockAng_deg, coords, offsets )
%grating = generateBinaryGrating( gratingAmp, period_samp, clockAng_deg, coords, offsets )
%   Generates array with representation of a binary phase grating 
%

% original example:
% gratingPhz = pi/2*(sign(sin(2*pi*(coordsFP.RHO.*cos(coordsFP.THETA - deg2rad(clockAng) )/lambda0Fnum_samp-shiftLam0OverD)*cyclesPerLam0OverD))+1);

    grating = zeros(coords.N,coords.N,numel(gratingAmp));
    
    Xoffset = coords.X-offsets(1);
    Yoffset = coords.Y-offsets(2);
    
    [THETA,RHO] = cart2pol(Xoffset,Yoffset);
    
    clockAng_rad = deg2rad(clockAng_deg);
    
    for index = 1:numel(gratingAmp)
        grating0 = 0.5*(sign(sin(2*pi*RHO.*cos(THETA -clockAng_rad )/period_samp))+1);
        grating(:,:,index) = gratingAmp(index)*grating0;
    end

end

