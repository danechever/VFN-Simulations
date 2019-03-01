function mask = generateVortexMask( charge, coords, offsets )
%mask = generateVortexMask( charge, coords, offsets )
%   Generates array with representation of a vortex phase mask
%
    
    mask = zeros(coords.N,coords.N,numel(charge));
    
    for index = 1:numel(charge)
        mask(:,:,index) = exp(1i*charge(index)*atan2(coords.Y-offsets(2),coords.X-offsets(1)));
    end

end

