clear; %close all; 
addpath('VFNlib');

%% Input parameters 

% Define smapling info
N = 2^12; % Size of computational grid (NxN samples) 
apRad = 256; % Aperture radius in samples 

% Define wavelength info
lambda0 = 635e-9; %central wavelength
fracBW = 0.4; %\Delta\lambda/\lambda
numWavelengths = 5;% number of discrete wavelengths 
lambdas = getWavelengthVec(lambda0,fracBW,numWavelengths);% array of wavelengths (meters)

% Define charge of the vortex mask at each wavelength
%charge = ones(1,numWavelengths); % achromatic
charge = lambda0./lambdas; % simple scalar model

% Define wavefront error at the central wavelength
nolls = 4:8;
coeffs = [0.01,-0.0099,-0.0095,-0.0008,0.0033];

% Give offsets for the vortex mask
offsetX = 0;%0.0952*apRad;
offsetY = 0;%0.0524*apRad; 

%% Generate the coordinate system

coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
xvals = coords.xvals;% Helpful for plotting
yvals = coords.yvals;

%% Create array with pupil function

PUPIL = makeCircularPupil( apRad, N );
[normI, totalPower0] = getNormalization(PUPIL);% Normalization factors
lambdaOverD = N/apRad/2; % lam/D in units of samples in the image plane

% figure(1)
% imagesc(xvals/apRad,yvals/apRad,PUPIL);
% axis image; 
% axis([-1 1 -1 1]);
% title('Pupil');
% colorbar; 
% colormap(parula(256));
% drawnow;

%% Define pupil field

phz = generateZernike_fromList( nolls, coeffs, PUPIL, apRad, coords); 

figure(2);
for ch = 1:numWavelengths
    Epup(:,:,ch) = exp(1i*phz*lambda0/lambdas(ch)).*PUPIL;
    
    subplot(1,numWavelengths,ch);
    imagesc(xvals/apRad,yvals/apRad,angle(Epup(:,:,ch)));
    axis image; 
    axis([-1 1 -1 1]);
    title(['Phase at ',num2str(lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(parula(256));
   
end
drawnow;

%% Get PSF without vortex mask

% Get broadband PSF
iPSF_BB = getPSF(Epup,lambda0,lambdas,normI,coords);

figure(3)
imagesc(xvals/lambdaOverD,yvals/lambdaOverD,iPSF_BB);
axis image; 
axis([-3 3 -3 3]);
title('broadband PSF w/o vortex');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Make vortex mask 

EPM = generateVortexMask( charge, coords, [offsetX offsetY] );

central_band_index = ceil(numWavelengths/2);

figure(4)
imagesc(xvals/apRad,yvals/apRad,angle(EPM(:,:,central_band_index).*Epup(:,:,central_band_index)));
axis image; 
axis([-1 1 -1 1]);
title('Pupil phase at \lambda_0');
colormap(hsv(256));
colorbar; 
drawnow;

%% Get PSF with vortex mask

iPSFv_BB = getPSF(Epup.*EPM,lambda0,lambdas,normI,coords);

figure(5)
imagesc(xvals/lambdaOverD,yvals/lambdaOverD,iPSFv_BB);
axis image; 
axis([-3 3 -3 3]);
title('broadband PSF w/ vortex');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Generate coupling maps for each wavelength

% Parameters for Thorlabs SM600
    % link: https://www.thorlabs.com/NewGroupPage9_PF.cfm?ObjectGroup_ID=949
fiber_props.core_rad = 4.3e-6/2;% Core radius [um]
fiber_props.n_core = 1.4606;% core index (interpolated from linear fit to 3 points)
fiber_props.n_clad = 1.4571;% cladding index (interpolated from linear fit to 3 points)
fiber_props.Fnum = 5; % focal ratio of the beam at the fiber
fiber_props.type = 'bessel';

eta_maps = generateCouplingMap_polychromatic( Epup.*EPM, fiber_props, lambda0, lambdas, totalPower0, lambdaOverD, 3*lambdaOverD, coords);

figure(6);
for ch = 1:numWavelengths
    subplot(1,numWavelengths,ch);
    imagesc(xvals/lambdaOverD,yvals/lambdaOverD,eta_maps(:,:,ch));
    axis image; 
    axis([-3 3 -3 3]);
    title(['\eta at ',num2str(lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(gray(256));
end