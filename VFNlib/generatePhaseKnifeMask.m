function mask = generatePhaseKnifeMask( phzStep, clockAng_deg, coords, offsets )
%mask = generatePhaseKnifeMask( phzStep, clockAng, coords, offsets )
%   Generates array with representation of a phase knife mask
%
    
    mask = zeros(coords.N,coords.N,numel(phzStep));
    
    Xoffset = coords.X-offsets(1);
    Yoffset = coords.Y-offsets(2);
    
    [THETA,RHO] = cart2pol(Xoffset,Yoffset);
    
    clockAng_rad = deg2rad(clockAng_deg);
    
    for index = 1:numel(phzStep)
        phz = phzStep(index)*(RHO.*cos(THETA-clockAng_rad)>0);
        mask(:,:,index) = exp(1i*phz);
    end

end

