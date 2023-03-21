%HOME PATH FOR THIS FILE - \VFN-Simulations\ScalarOnPrimary_Simulations

%-- Coordinates in the focal plane
coordsPP = generateCoordinates(inpar.N);% Creates NxN arrays with coordinates 
% Helpful for plotting: get axes coords in pupil radii
inpar.xvalsPP = coordsPP.xvals/inpar.apRad;
inpar.yvalsPP = coordsPP.yvals/inpar.apRad;
% Used in generateVortexMaskKeckPrimary
inpar.xvals = inpar.xvalsPP;
inpar.yvals = inpar.yvalsPP;

%-- Coordinates in the pupil plane
Nxi = im_size*lambda0Fnum_samp;
coordsFP = generateCoordinates( Nxi );
% Helpful for plotting: get axes coords in lam/D 
inpar.xvalsFP = coordsFP.xvals/lambda0Fnum_samp;
inpar.yvalsFP = coordsFP.yvals/lambda0Fnum_samp;

%-- Key values for getting scaling right using falco MFT propagator
% Lambda over D of central wavelength in meters at final focal plane 
lambda0Fnum_meters = inpar.lambda0*foc/DPup;     
% "pixel" size in pupil and image planes
coordsPP.dx  = DPup/(2*inpar.apRad);   % meters/sample
coordsFP.dx = lambda0Fnum_meters/lambda0Fnum_samp;     % meters/sample