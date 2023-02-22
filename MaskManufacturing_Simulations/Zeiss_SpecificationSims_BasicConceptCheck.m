%{
    Script for simulating a step-wise vortex pattern as was manufactured by Zeiss for HCST SVC projects

    The goal was to confirm that a step-wise mask still works conceptually for VFN. It does a normal set 
    of VFN simulations with the usual plots (Null vs wavelength, peak coupling vs. wavelength, rad profile, etc.)

    Note: This script is also a helpful basis for other simple scripts wishing to test basic pupil-VFN performance.
    With few mods, you will see all the key metrics for whatever system youd like to test.
%}

clear; close all; 
addpath(genpath(fullfile('..','VFNlib')));
addpath(genpath(fullfile('..','..','falco-matlab')))

%% Input parameters 

%-- Provide regular parameters
% Define smapling info
N = 2^10; % Size of computational grid (NxN samples) 
apRad = N/2-4; % Aperture radius in samples 

%-- Define wavelength info
lambda0 = 2.2e-6; %central wavelength
fracBW = 0.2; % \Delta\lambda/\lambda
numWavelengths = 3; % number of discrete wavelengths 
lambdas = getWavelengthVec(lambda0,fracBW,numWavelengths);% array of wavelengths (meters)

%-- Define charge of the vortex mask 
charge = 1;     % Since using falco, provide a single wavelength
phaseScaleFac = lambda0./lambdas;    % Since using falco, provide the scaling factor separately
Nsteps = 10;    % Number of steps to use in the staircase model
vortClocking = 0;   % [deg] Clocking angle for the vortex 
isStaircase = false;  % Optional flag to use staircase model instead of normal scalar model used by VFN-Sims library

%-- Define wavefront error at the central wavelength
  % 1) Pist, Tip, Tilt, defoc, ob.astig, ver.astig, ver.coma, hor.coma,
  % 9) ver.tref, ob.tref, spher
nolls = 4:8;
coeffs = [0.0, 0.0, 0.0, 0.0, 0.0];  % [Waves RMS]

%-- Give offsets for the vortex mask
offsetX = 0;    % [samples in the pupil-plane]
offsetY = 0; 

%-- Parameters for SMF (Thorlabs SM2000 in this case)
%     % link: https://www.thorlabs.com/NewGroupPage9_PF.cfm?ObjectGroup_ID=949
fiber_props.core_rad = 11e-6/2;% Core radius [um]
fiber_props.n_core = 1.4436;% core index 
fiber_props.n_clad = 1.4381;% cladding index 
fiber_props.type = 'gaussian';  

%-- Define parameters for falco MFT propagator
    % these don't matter in themselves as long as they are consistent w/ each 
    % other and with lambda
% OPTION 1: define focal length and pup diameter manually
%foc = 11e-3;      %[m] final focal length
DPup = 13e-3;    %[m] pupil size 
% OPTION 2: solve for focal length based on ideal Fnum and pup diameter
Fnum = getMFD(fiber_props,lambda0)/(lambda0*1.42); % focal ratio of the beam at the fiber
foc = Fnum*DPup;
fprintf('Focus in use: %f [mm]\n',foc*1e3)
fprintf('F/# in use:      %f \n',Fnum)
fprintf('Pup Diam in use: %f [mm]\n',DPup*1e3)

% Since using falco MFT propogator, can choose our final image sampling
lambda0Fnum_samp = 50; %[samples/ (lam0 F#)] in units of samples in the image plane
im_size = 10;    %[lam0/D] field of view in final image plane

%% Generate the coordinate system

%-- Coordinates in the focal plane
coordsPP = generateCoordinates(N);% Creates NxN arrays with coordinates 
% Helpful for plotting: get axes coords in pupil radii
xvalsPP = coordsPP.xvals/apRad;
yvalsPP = coordsPP.yvals/apRad;

%-- Coordinates in the pupil plane
Nxi = im_size*lambda0Fnum_samp;
coordsFP = generateCoordinates( Nxi );
% Helpful for plotting: get axes coords in lam/D 
xvalsFP = coordsFP.xvals/lambda0Fnum_samp;
yvalsFP = coordsFP.yvals/lambda0Fnum_samp;

%-- Key values for getting scaling right using falco MFT propagator
% Lambda over D of central wavelength in meters at final focal plane 
lambda0Fnum_meters = lambda0*foc/DPup;     
% "pixel" size in pupil and image planes
coordsPP.dx  = DPup/(2*apRad);   % meters/sample
coordsFP.dx = lambda0Fnum_meters/lambda0Fnum_samp;     % meters/sample

%% Create array with pupil function

PUPIL = makeKeckPupil( 2*apRad, N );

%-- Get norm for coupling fractions (simple sum since using MFT propagator)
totalPower0 = sum(abs(PUPIL(:)));

figure(1); 
imagesc(xvalsPP,yvalsPP,PUPIL); 
axis image;
axis([-1 1 -1 1]);
title('Pupil');
colorbar; 
colormap(parula(256));
drawnow;


%% Define pupil field
phz = generateZernike_fromList( nolls, coeffs, PUPIL, apRad, coordsPP); 

figure(2);
Epup = nan(N,N,numWavelengths);
for ch = 1:numWavelengths
    Epup(:,:,ch) = exp(1i*phz*lambda0/lambdas(ch)).*PUPIL;
    
    subplot(1,numWavelengths,ch);
    imagesc(xvalsPP,yvalsPP,angle(Epup(:,:,ch)));
    axis image; 
    axis([-1 1 -1 1]);
    title(['Phase at ',num2str(lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(parula(256));
   
end
drawnow;

%% Get PSF without vortex mask
[iPSF_BB, PSF] = getPSF_mft(Epup, lambdas, foc, coordsPP, coordsFP);

figure(3)
imagesc(xvalsFP,yvalsFP,iPSF_BB);
axis image; 
%axis([-3 3 -3 3]);
title('broadband PSF w/o vortex');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Make vortex mask (using staircase from falco)
if isStaircase
    fprintf('-- Using Staircase Model\n\n')
    % Required Inputs
    inputs.type = 'staircase';
    inputs.N = N; % number of pixels across the array
    inputs.charge = charge; % charge of the mask (makes most sense for vortex)
    inputs.Nsteps = Nsteps; % number of steps per 2*pi radians. For 'staircase' only

    % Optional Inputs
    inputs.centering = 'pixel';
    inputs.xOffset = offsetX; % [pixels]
    inputs.yOffset = offsetY; % [pixels]
    inputs.clocking = vortClocking; % [degrees]

    % Generate the vortex (iteratively at each wavelength)
    EPM = zeros(N,N,numWavelengths);
    for ch = 1:numWavelengths
        % Set the scaling factor for the given wavelength
        inputs.phaseScaleFac = phaseScaleFac(ch); % Factor to apply uniformly to the phase. Used to add chromaticity.
        EPM(:,:,ch) = falco_gen_azimuthal_phase_mask(inputs);
    end
else
    fprintf('-- Using Simple Model\n\n')
    EPM = generateVortexMask( charge.*phaseScaleFac, coordsPP, [offsetX offsetY] );
end
    
central_band_index = ceil(numWavelengths/2);

figure(4)
imagesc(xvalsPP,yvalsPP,angle(EPM(:,:,central_band_index).*Epup(:,:,central_band_index)));
axis image; 
axis([-1 1 -1 1]);
title('Pupil phase at \lambda_0');
colormap(hsv(256));
colorbar; 
drawnow;

%% Get PSF with vortex mask
[iPSFv_BB, PSFv] = getPSF_mft(Epup.*EPM, lambdas, foc, coordsPP, coordsFP);


figure(5)
imagesc(xvalsFP,yvalsFP,iPSFv_BB);
axis image; 
axis([-3 3 -3 3]);
title('broadband PSF w/ vortex');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;


%% Generate fibermode at each lambda     
%-- Iterate through wavelengths generating modes
fibmode = nan(Nxi, Nxi, numWavelengths);
for ch = 1:numWavelengths
    % Generate the fiber mode for this wavelength with proper scaling
	fibmode(:,:,ch) = generateSMFmode_mft( fiber_props, lambdas(ch), coordsFP.dx, coordsFP);
    
end

%% Calculate Throughputs
%-- Get null depth (using overlap integral)
eta_onAx = nan(1,numWavelengths);
for ch = 1:numWavelengths
    eta_onAx(ch) = (abs(sum(sum(PSFv(:,:,ch).*fibmode(:,:,ch)))).^2)/totalPower0;
end

%-- Compute 2D coupling map
eta_maps = zeros(Nxi,Nxi,numWavelengths);
for ch = 1:numWavelengths
    % Compute the monochromatic coupling map (spatial units of lambda/D)
    eta_maps(:,:,ch) = generateCouplingMap( fibmode(:,:,ch), PSFv(:,:,ch), totalPower0, 5*lambda0Fnum_samp, coordsFP);

end

%-- Compute azimuthally-averaged performance
for ch = 1:numWavelengths
    % Compute average
    [tmpAvg, qvec] = radialAverage(eta_maps(:,:,ch), [Nxi/2+1 Nxi/2+1]);
    
    % Yes, I know it's not efficient to recreate matrix on every itr., I
    % don't care in this particular application...
    eta_pAvgs(:,ch) = tmpAvg;
end
% Since all eta_maps have the same spatial scale, they all have the same qvec
qvec = qvec/lambda0Fnum_samp; % [lam/D] (converted from pix to lam/D directly)

disp('Key Coupling Points:')
for ch = 1:numWavelengths
    fprintf('lambda = %f nm,    on-axis coup = %e,    peak radAvg = %f %%\n',lambdas(ch)*1e9, eta_onAx(ch), max(eta_pAvgs(:,ch))*100);
end

%% Display coupling maps
%-- Linear scale
figure(6);
for ch = 1:numWavelengths
    subplot(1,numWavelengths,ch);
    imagesc(xvalsFP,yvalsFP,eta_maps(:,:,ch));
    axis image; 
    axis([-3 3 -3 3]);
    title(['\eta at ',num2str(lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(gray(256));
end
drawnow

%-- Log Scale
figure(7);
for ch = 1:numWavelengths
    subplot(1,numWavelengths,ch);
    imagesc(xvalsFP,yvalsFP,log10(eta_maps(:,:,ch)));
    axis image; 
    axis([-3 3 -3 3]);
    title(['log10(\eta) at ',num2str(lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(gray(256));
end
drawnow

%-- Log Scale (focusing on the null region)
figure(8);
for ch = 1:numWavelengths
    subplot(1,numWavelengths,ch);
    imagesc(xvalsFP,yvalsFP,log10(eta_maps(:,:,ch)));
    axis image; 
    axis([-0.25 0.25 -0.25 0.25]);
    title(['log10(\eta) at ',num2str(lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(gray(256));
end
drawnow


%% Display broadband coupling map if applicable
% This is what the photodiode should see when broadband light is applied

if numWavelengths > 1
    % Compute BB coupling map as average of all coupling maps
      % (This assumes the input spectrum is is flat --> ie. equal power at
      % all wavelengths)
    eta_map_BB = mean(eta_maps,3);
    
    % Compute broadband average radial coupling
    eta_pAvg_BB = mean(eta_pAvgs, 2);
    
    % Plot Linear Scale
    figure(9);
    imagesc(xvalsFP, yvalsFP, eta_map_BB);
    axis image;
    axis([-3 3 -3 3]);
    title('Broadband \eta');
    colorbar;
    colormap(gray(256));
    
    % Plot Log Scale
    figure(10);
    imagesc(xvalsFP, yvalsFP, log10(eta_map_BB));
    axis image;
    axis([-3 3 -3 3]);
    title('Broadband log10(\eta)');
    colorbar;
    colormap(gray(256));
    
    % Print key values
    disp('---')
    fprintf('Broadband Performance:     on-axis coup = %e,    peak radAvg = %f %%\n', mean(eta_onAx), max(eta_pAvg_BB)*100);
end