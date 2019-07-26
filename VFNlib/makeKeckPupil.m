function [PUPIL,circum_diam] = makeKeckPupil( N_flat2flat, Ngrid )
%PUPIL = makeKeckPupil( N_flat2flat, Ngrid )
%   Generates the Keck pupil.
%
%   Inputs: 
%       N_flat2flat - Number of samples across the keck aperture. Measured 
%                     flat to flat. 
%       Ngrid - Number of samples in the computational grid (NxN).
%
%   Returns:
%       PUPIL - 2D array containing the Keck pupil
%       circum_diam - Circumscribed diameter in pixels
    
    addpath([fileparts(which(mfilename)),'/segMirrorFunctions']);
    apDia0 = N_flat2flat; % flat to flat diameter in samples 
    gapWidth = 25.4/10916*apDia0/2;% gap width in samples 
    numRings = 3;% number of rings of hexagonal segments (for Keck, numRings=3)

    [X,Y] = meshgrid(-Ngrid/2:Ngrid/2-1);
    [THETA,RHO] = cart2pol(X,Y);
    
    hexMirror.apDia = apDia0; % flat to flat aperture diameter (samples)
    hexMirror.wGap = gapWidth; % samples
    hexMirror.numRings = numRings;% Number of rings in hexagonally segmented mirror 
    hexMirror.Npad = Ngrid;

    PUPIL = hexSegMirror_getSupport( hexMirror ); 

    [rowslist,colslist] = find(round(PUPIL));
    apRad = max(sqrt((rowslist-Ngrid/2-1).^2 + (colslist-Ngrid/2-1).^2));
    circum_diam = 2*apRad;

    secRad = apRad*2600/10949;
    spwidth = 2*apRad*25.4/10949;
    PUPIL = PUPIL.*(1-exp(-(RHO/secRad).^1000));
    PUPIL = PUPIL.*(1-exp(-(Y/(spwidth/2)).^1000));
    PUPIL = PUPIL.*(1-exp(-(RHO.*cos(THETA-30*pi/180)/(spwidth/2)).^1000));
    PUPIL = PUPIL.*(1-exp(-(RHO.*cos(THETA+30*pi/180)/(spwidth/2)).^1000));

end

