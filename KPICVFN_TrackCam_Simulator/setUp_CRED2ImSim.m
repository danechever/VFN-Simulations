function cred2sim = setUp_CRED2ImSim(Npup, Nim, charge, lambda0, starHmag)
%cred2sim = setUp_CRED2ImSim(N, charge, lambda0, starHmag)
%   Setup the KPIC Phase II CRED2 Image simulator
%
%   Inputs:
%       Npup: size of computational grid in pupil plane (NxN samples)
%       Nxi: size of computation grid in focal plane (Nxi x Nxi samples)
%       charge: vortex charge at central wavelength
%       lambda0: central wavelength (meters)
%       starHmag: H-band magnitude of the star
%   Outputs:
%       cred2sim: struct with elements needed for the simulator

%-- Provide regular parameters
% Define smapling info
%N = 2^7; % Size of computational grid (NxN samples) 
cred2sim.apRad = Npup/2-4; % Aperture radius in samples 

%-- Define wavelength info
%lambda0 = 1640e-9; %central wavelength
cred2sim.lambda0 = lambda0;
fracBW = 1e-9;%0.2; % \Delta\lambda/\lambda
cred2sim.numWavelengths = 1; % number of discrete wavelengths 
cred2sim.lambdas = getWavelengthVec(cred2sim.lambda0,fracBW,cred2sim.numWavelengths);% array of wavelengths (meters)

%-- Define charge of the vortex mask at each wavelength
cred2sim.charge = charge*ones(1,cred2sim.numWavelengths); % achromatic

%-- Give offsets for the vortex mask
offsetX = 0;    % [samples in the pupil-plane]
offsetY = 0; 

%-- Define parameters for falco MFT propagator
    % these don't matter in themselves as long as they are consistent w/ each 
    % other and with lambda
fnum = 38.9; % from Mitsuko Zemax measurement (link=https://caltech.sharepoint.com/sites/coo/Shared%20Documents/OIR%20-%20Preprojects/KPIC/KPIC%20-%20Systems%20Engineering%20%5BL2%5D/KPIC%20-%20L2%20Optical%20Design/Optical_simulation_reports/KPIC%20phase%20II%20TTM%20and%20focal%20plane%20motion.pptx?d=w16f8e6de681241769a41208ed9388f89&csf=1&web=1&e=fAYtnq)
             % (Alternate value of unkown origin: 23.9)
cred2sim.DPup = 1e-3;    %[m] pupil size 
cred2sim.foc = fnum*cred2sim.DPup;      %[m] final focal length

%-- Key values for getting scaling right using falco MFT propagator
% Lambda over D of central wavelength in meters at final focal plane 
lambda0Fnum_meters = cred2sim.lambda0*cred2sim.foc/cred2sim.DPup;     
% "pixel" size in pupil planes
PUPIL_dx  = cred2sim.DPup/(2*cred2sim.apRad);   % meters/sample
%-- Set the final image size since Falco lets us set any value manually
% OPTION 1: Provide value in lambda/D and solve for number of pixels
%im_size = 10;    %[lam0/D] field of view in final image plane
  % NOTE: ceil() is important to guarantee square matrix
%Nxi = im_size*ceil(lambda0Fnum_samp);
% OPTION 2: Provide number of pixels directly
%Nxi = 64; % 64x64 is a common cropping window I use on the CRED2 for TT tests

% OPTION 1: set lambda0Fnum sampling and compute pixel size
% lambda0Fnum_samp = 50; %[samples/ (lam0 F#)] in units of samples in the image plane
% CRED2_dx = lambda0Fnum_meters/lambda0Fnum_samp;     % meters/sample
% OPTION 2: set pixel size and compute lambda0Fnum sampling
cred2sim.CRED2_dx = 15e-6;      % [m] physical pixel size
cred2sim.lambda0Fnum_samp = lambda0Fnum_meters/cred2sim.CRED2_dx;

% Compute lambda0OverD_rad for Keck
cred2sim.Keck_D = 10.949; % [m] circumscribed diameter of Keck pupil
lambda0OverKeckD_rad = lambda0/cred2sim.Keck_D;
% Now convert to mas
cred2sim.lambda0OverKeckD_mas = lambda0OverKeckD_rad*206265*1e3;    


%-- Define remaining CRED2 properties
cred2sim.CRED2_QE = 0.85;  % [e-/photon]
cred2sim.CRED2_DarkCurrent = 449;    % [e-/s/pix] (value @-40C)
cred2sim.CRED2_ReadNoiseRMS = 39.6;  % [e- rms (/pix?)]
cred2sim.CRED2_ADU = 2.26;   % [e-/count]

%% Generate the coordinate system

%-- Coordinates in the focal plane
coordsPP = generateCoordinates(Npup);% Creates NxN arrays with coordinates 
% Helpful for plotting: get axes coords in pupil radii
cred2sim.xvalsPP = coordsPP.xvals/cred2sim.apRad;
cred2sim.yvalsPP = coordsPP.yvals/cred2sim.apRad;

%-- Coordinates in the pupil plane
coordsFP = generateCoordinates( Nim );
% Helpful for plotting: get axes coords in lam/D 
cred2sim.xvalsFP = coordsFP.xvals/cred2sim.lambda0Fnum_samp;
cred2sim.yvalsFP = coordsFP.yvals/cred2sim.lambda0Fnum_samp;

coordsPP.dx = PUPIL_dx;
coordsFP.dx = cred2sim.CRED2_dx;

% Store into simulator struct
cred2sim.coordsPP = coordsPP;
cred2sim.coordsFP = coordsFP;

%% Create array with pupil function

cred2sim.PUPIL = makeCircularPupil( cred2sim.apRad, Npup );

% figure(1); 
% imagesc(cred2sim.xvalsPP,cred2sim.yvalsPP,cred2sim.PUPIL); 
% axis image;
% axis([-1 1 -1 1]);
% title('Pupil');
% colorbar; 
% colormap(parula(256));
% drawnow;

%% Deal with flux stuff

%-- Get magnitude and corresponding Flux
% Define star magnitude
%starHmag = 4.26;
cred2sim.starHmag = starHmag;

% Define spectral band
Hlam_cen = 1.64e-6;     % [m]   central wavelength for Astro band
Hlam_del = 0.33e-6;    % [m]    bandwidth for Astro band
Hlam_low = Hlam_cen - Hlam_del/2;   % [m] low-end cutoff for Astro band
Hlam_hi = Hlam_cen + Hlam_del/2;    % [m] hi-end cutoff for Astro band
CRED2_Hlam_hi = 1.7e-6; % [m] Upper cutoff of CRED2 detector
%Compute the bands
Hlam_astro_band = Hlam_hi-Hlam_low;     % [m] width of Astro band
CRED2_Hband_fracBand = (CRED2_Hlam_hi-Hlam_low)/Hlam_astro_band; % fractional width of effective CRED2 band

% Define H-band flux for 0'th mag star
Hflux_mag0 = 3.22e9;    % [ph/s/m2]     Flux in Astro band
Hflux_mag0 = CRED2_Hband_fracBand * Hflux_mag0;     % [ph/s/m2] Flux in effective CRED2 band

% Scale the flux from 0'th mag to our magnitude
flux_at_mag = Hflux_mag0 * 10^(-cred2sim.starHmag/2.5); 


%-- Distribute flux over the pupil
% Determine the flux collected by Keck aperture
Keck_A = 72.341;    % [m2] Keck Primary collecting area
flux_in_beam = Keck_A * flux_at_mag;    % [ph/s] flux 
% Account for instrument throughput losses
cred2sim.KPIC_thpt = 0.01;%0.0268;     % Keck+KPIC throughput to the CRED2
flux_in_beam = cred2sim.KPIC_thpt * flux_in_beam;   % [ph/s]

% Count number of "pixels" within the pupil in sim
    % NOTE: due to apodization at edges of pupil, this is a non-integer but
    % it's good enough and a valid representation of flux in the pupil
litPupilPix = sum(abs(cred2sim.PUPIL(:)));   % [pix]
% Distribute flux over all lit "pixels" 
flux_per_Pupil_pix = flux_in_beam / litPupilPix;    % [ph/s/pix]

% Take sqrt() of flux per pixel to get E-field amplitude per pixel 
    % REMEMBER THAT THE E-field is the sqrt() of the flux!!
EPup_amp = sqrt(flux_per_Pupil_pix); % [sqrt(ph/s/pix)]


% Apply to PUPIL variable so that E-field is scaled correctly
    % NOTE: this doesn't affect the calls to logical() in some functions
    % (eg. genZern_...) since logical() will convert any non-0 value to a 1
cred2sim.PUPIL = cred2sim.PUPIL * EPup_amp;

fprintf('Flux per pupil pix: %f [ph/s/pix]\n',flux_per_Pupil_pix);
fprintf('E-field Amp at Pup: %f [sqrt(ph/s/pix)]\n',EPup_amp);

%% Make vortex mask
cred2sim.EPM = generateVortexMask( cred2sim.charge, cred2sim.coordsPP, [offsetX offsetY] );

% central_band_index = ceil(cred2sim.numWavelengths/2);
% 
% figure(4)
% imagesc(cred2sim.xvalsPP,cred2sim.yvalsPP,angle(cred2sim.EPM(:,:,central_band_index).*cred2sim.PUPIL));
% axis image; 
% axis([-1 1 -1 1]);
% title('Pupil phase at \lambda_0');
% colormap(hsv(256));
% colorbar; 
% drawnow;


end