% Demo of a simple wedge implementation

clear; close all;

% Set Paths
addpath(genpath(fullfile('..','VFNlib')));
addpath(genpath(fullfile('..','..','falco-matlab')));
addpath(genpath(fullfile('..','..','VFN-Lab')));

%% Input parameters 

%-- Provide regular parameters
% Define sampling info
N = 2^10; % Size of computational grid (NxN samples) 
apRad = N/2-4; % [samples] Aperture radius

%-- Define wavelength info
lambda0 = 2.2e-6; % [m] central wavelength
fracBW = 0.1818; % [\Delta\lambda/\lambda] fractional bandwidth
numWavelengths = 5; % number of discrete wavelengths 
lambdas = getWavelengthVec(lambda0, fracBW, numWavelengths);% [m] array of wavelengths

%-- Define parameters for falco MFT propagator
    % these don't matter in themselves as long as they are consistent w/ each 
    % other and with lambda
% OPTION 1: define focal length and pup diameter manually
%foc = 11e-3;      %[m] final focal length
keckD = 10.949;   %[m] pupil size 
% OPTION 2: solve for focal length based on ideal Fnum and pup diameter
Fnum = 4.6; % focal ratio of the beam at the final focal plane
foc = Fnum*keckD;
fprintf('Focus in use: %f [mm]\n',foc*1e3)

mag = 1;
beamD = mag*keckD;

% Since using falco MFT propogator, can choose our final image sampling
lambda0Fnum_samp = 50; %[samples/ (lam0 F#)] in units of samples in the image plane
im_size = 10;    %[lam0/D] field of view in final image plane

%% Wedge & ADC initial parameters

wType = 'wedge';    % options: 'wedge', 'adc'

wedge_mat = 'caf2';     % wedge material ('NBK7', 'CaF2', 'BaF2', 'ZnSe')
%n_wedge = getRefractiveIndex(wedge_mat, lambda0);   % refractive index at central wavelength
%beam_dev = 1.1*lambda0/beamD;    %[rads] Beam deviation to simulate at central wavelength
wedge_angle = deg2rad(-75.9266);%atan2(beam_dev,(n_wedge-1));   % [rads] wedge angle to use
        % equation from max: atan(890.16*(inpar.p*inpar.lambda0)/(inpar.keckD*(getRefractiveIndex(inpar.wedge_mat ,1e6*inpar.lambda0) - 1)));

%% Generate the coordinate system

%-- Coordinates in the focal plane
coordsPP = generateCoordinates(N);% Creates NxN arrays with coordinates 
% Helpful for plotting: get axes coords in pupil radii
xvalsPP = coordsPP.xvals/apRad;
yvalsPP = coordsPP.yvals/apRad;

%-- Coordinates in the pupil plane
Nxi = im_size*lambda0Fnum_samp;
%inpar.Nwedge = Nxi;
coordsFP = generateCoordinates( Nxi );
% Helpful for plotting: get axes coords in lam/D 
xvalsFP = coordsFP.xvals/lambda0Fnum_samp;
yvalsFP = coordsFP.yvals/lambda0Fnum_samp;

%-- Key values for getting scaling right using falco MFT propagator
% Lambda over D of central wavelength in meters at final focal plane
lambda0Fnum_meters = lambda0*foc/keckD; % [m]

% "pixel" size in pupil and image planes
coordsPP.dx = keckD/(2*apRad);   % [meters/sample]
coordsFP.dx = lambda0Fnum_meters/lambda0Fnum_samp;     % [meters/sample]

%% Create array with pupil function

PUPIL = makeCircularPupil(apRad, N);

figure();
imagesc(xvalsPP,yvalsPP,PUPIL); 
axis image;
title('Pupil');
colorbar; 
colormap(parula(256));
drawnow;

%% Get PSF without wedge
% Replicate pupil at each wavelength for this since all have the same field
Epup = repmat(PUPIL, [1,size(lambdas)]);
[~, PSF] = getPSF_mft(Epup, lambdas, foc, coordsPP, coordsFP);


figure();
pause(0.1);
for ch = 1:numWavelengths    
    subplot(1,numWavelengths,ch);
    iPSF = abs(PSF(:,:,ch)).^2; % Compute intensity (ie. image)
    imagesc(xvalsFP,yvalsFP,iPSF);
    axis image; 
    axis xy;
    axis([-3 3 -3 3]);
    title(['PSF at ',num2str(lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(parula(256));
    drawnow;
end
%% Generate the wedge

wphz = nan(size(Epup));
figure()
for ch = 1:numWavelengths
    % Get refractive index of wedge at this wavelength
    n_wedge_ch = getRefractiveIndex(wedge_mat, lambdas(ch));   % refractive index at central wavelength
    
    % Get wedge at this wavelength
    wphz(:,:,ch) = simpleWedge_WIP(wedge_angle, n_wedge_ch, lambdas(ch), coordsPP);
    
    % Plot wedge phase 
    subplot(1,numWavelengths,ch);
    imagesc(xvalsPP,yvalsPP,wphz(:,:,ch)/2/pi.*PUPIL);   % Use .*PUPIL to get within beam
    axis image; 
    axis xy;
    title(['Wedge Phase at ',num2str(lambdas(ch)*1e9),'nm']);
    cb = colorbar; 
    cb.Label.String = 'Fractional wavelength';
    colormap(hsv(256));
    drawnow;
    disp(lambdas(ch));
    disp(max(max(wphz(:,:,ch)/2/pi))-min(min(wphz(:,:,ch)/2/pi)));
end


wedge = exp(1i*wphz);

%% Get PSF with the wedge

[~, PSFw] = getPSF_mft(Epup.*wedge, lambdas, foc, coordsPP, coordsFP);

figure();
pause(0.1);
for ch = 1:numWavelengths    
    subplot(1,numWavelengths,ch);
    iPSFw = abs(PSFw(:,:,ch)).^2; % Compute intensity (ie. image)
    imagesc(xvalsFP,yvalsFP,iPSFw);
    axis image; 
    axis xy;
    axis([-3 3 -3 3]);
    title(['PSF at ',num2str(lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(parula(256));
    drawnow;
end