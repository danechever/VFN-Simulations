function mask = generateVortexMask( charge, coords, offsets )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

mask = exp(1i*charge*atan2(coords.Y-offsets(2),coords.X-offsets(1)));

end

