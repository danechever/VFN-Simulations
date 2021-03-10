% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
function [ arrayOut ] = hexSegMirror_addHexSegment_Modified( cenrow, cencol, apDia,piston, tiltx, tilty, arrayIn, X, Y, hexAmpConst, hexPhzConst)
%hexSegMirror_addHexSegment Adds hexagonal mirror segment to arrayIn, 
% centered at (cenrow, cencol). The full mirror have numRings rings of 
% hexagonal segments, flat-to-flat diameter (in samples) of apDia, 
% wGap (in samples) between segments, and piston, tiltx, and tilty
% phase offsets. Piston is in units of waves. tiltx and tilty are waves
% across the full flat-to-flat pupil diameter apDia. 
%
%   Inputs:
%   cenrow - row of hexagon center (samples)
%   cencol - column of hexagon center (samples)
%   numRings - number of rings in the segmented mirror (samples)
%   apDia - flat to flat aperture diameter (samples)
%   wGap - width of the gap between segments (samples)
%   piston - Segment piston in waves
%   tiltx - Tilt on segment in horizontal direction (waves/apDia)
%   tilty - Tilt on segment in vertical direction (waves/apDia)
%   arrayIn - Input array
%   X,Y - cartesian coordinates of the plane that each segment lies in
%   hexAmpConst - the constant coefficient corresponding to amplitude
%   initialized at top of stack
%   hexPhzConst - the constant coefficient corresponding to phase
%   initialized at top of stack
%   
%   Coordinate system origin: (rows/2+1, cols/2+1)

    HEXamp = hexAmpConst;
    
    HEXphz = hexPhzConst...
             .*exp(1i*2*pi.*(piston+ ...
                             tiltx/apDia*(X-cencol)+...
                             tilty/apDia*(Y-cenrow)));
                         
%    HEXphzw = HEXphz.*exp(1i*2*pi.*(tiltw/apDia*(Y-cenrow)));
                         
%             .*exp(1i*2*pi*piston)...
%             .*exp(1i*2*pi*tiltx/apDia*(X-cencol))...
%             .*exp(1i*2*pi*tilty/apDia*(Y-cenrow));
        
    arrayOut = arrayIn + HEXamp.*HEXphz;
end

