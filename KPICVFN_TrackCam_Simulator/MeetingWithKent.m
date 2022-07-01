%{
%}

%clear; close all; 
addpath(genpath(fullfile('..','VFNlib')));
addpath(genpath(fullfile('..','..','falco-matlab')))

%% Input parameters 

%-- Provide regular parameters
% Define smapling info
N = 2^7; % Size of computational grid (NxN samples) 
apRad = N/2-4; % Aperture radius in samples 

%-- Define wavelength info
lambda0 = 1640e-9; %central wavelength
fracBW = 1e-9;%0.2; % \Delta\lambda/\lambda
numWavelengths = 1; % number of discrete wavelengths 
lambdas = getWavelengthVec(lambda0,fracBW,numWavelengths);% array of wavelengths (meters)

%-- Define charge of the vortex mask at each wavelength
charge = 0*ones(1,numWavelengths); % achromatic

%-- Define wavefront error at the central wavelength
  % 1) Pist, Tip, Tilt, defoc, ob.astig, ver.astig, ver.coma, hor.coma,
  % 9) ver.tref, ob.tref, spher
nolls = [9, 50:150];
coeffs = [0.09, 0.04*rand(1,length(nolls)-1)];  % [Waves RMS]

%-- Define tilts 
tilts = [0, -30.8955]; % (x, y) [mas] 

%-- Give offsets for the vortex mask
offsetX = 0;    % [samples in the pupil-plane]
offsetY = 0; 

%-- Define parameters for falco MFT propagator
    % these don't matter in themselves as long as they are consistent w/ each 
    % other and with lambda
fnum = 38.9; % from Mitsuko Zemax measurement (link=https://caltech.sharepoint.com/sites/coo/Shared%20Documents/OIR%20-%20Preprojects/KPIC/KPIC%20-%20Systems%20Engineering%20%5BL2%5D/KPIC%20-%20L2%20Optical%20Design/Optical_simulation_reports/KPIC%20phase%20II%20TTM%20and%20focal%20plane%20motion.pptx?d=w16f8e6de681241769a41208ed9388f89&csf=1&web=1&e=fAYtnq)
             % (Alternate value of unkown origin: 23.9)
             % (Alternate value from Kent: 32.87)
DPup = 1e-3;    %[m] pupil size 
foc = fnum*DPup;      %[m] final focal length

%-- Key values for getting scaling right using falco MFT propagator
% Lambda over D of central wavelength in meters at final focal plane 
lambda0Fnum_meters = lambda0*foc/DPup;     
% "pixel" size in pupil planes
PUPIL_dx  = DPup/(2*apRad);   % meters/sample
%-- Set the final image size since Falco lets us set any value manually
% OPTION 1: Provide value in lambda/D and solve for number of pixels
%im_size = 10;    %[lam0/D] field of view in final image plane
  % NOTE: ceil() is important to guarantee square matrix
%Nxi = im_size*ceil(lambda0Fnum_samp);
% OPTION 2: Provide number of pixels directly
Nxi = 90; % 64x64 is a common cropping window I use on the CRED2 for TT tests

% OPTION 1: set lambda0Fnum sampling and compute pixel size
% lambda0Fnum_samp = 50; %[samples/ (lam0 F#)] in units of samples in the image plane
% CRED2_dx = lambda0Fnum_meters/lambda0Fnum_samp;     % meters/sample
% OPTION 2: set pixel size and compute lambda0Fnum sampling
CRED2_dx = 15e-6;      % [m] physical pixel size
lambda0Fnum_samp = lambda0Fnum_meters/CRED2_dx;

%% Generate the coordinate system

%-- Coordinates in the focal plane
coordsPP = generateCoordinates(N);% Creates NxN arrays with coordinates 
% Helpful for plotting: get axes coords in pupil radii
xvalsPP = coordsPP.xvals/apRad;
yvalsPP = coordsPP.yvals/apRad;

%-- Coordinates in the pupil plane
coordsFP = generateCoordinates( Nxi );
% Helpful for plotting: get axes coords in lam/D 
xvalsFP = coordsFP.xvals/lambda0Fnum_samp;
yvalsFP = coordsFP.yvals/lambda0Fnum_samp;

coordsPP.dx = PUPIL_dx;
coordsFP.dx = CRED2_dx;

%% Create array with pupil function
% TODO::: Switch to keck pupil: makeKeckPupil
PUPIL = makeCircularPupil( apRad, N );

figure(1); 
imagesc(xvalsPP,yvalsPP,PUPIL); 
axis image;
axis([-1 1 -1 1]);
title('Pupil');
colorbar; 
colormap(parula(256));
drawnow;

%% Deal with flux stuff

%-- Get magnitude and corresponding Flux
% Define star magnitude
star_Hmag = 4.26;

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
flux_at_mag = Hflux_mag0 * 10^(-star_Hmag/2.5); 


%-- Distribute flux over the pupil
% Determine the flux collected by Keck aperture
Keck_A = 72.341;    % [m2] Keck Primary collecting area
flux_in_beam = Keck_A * flux_at_mag;    % [ph/s] flux 
% Account for instrument throughput losses
KPIC_thpt = 0.01; %0.0268;     % Keck+KPIC throughput to the CRED2
flux_in_beam = KPIC_thpt * flux_in_beam;   % [ph/s]

% Count number of "pixels" within the pupil in sim
    % NOTE: due to apodization at edges of pupil, this is a non-integer but
    % it's good enough and a valid representation of flux in the pupil
litPupilPix = sum(abs(PUPIL(:)));   % [pix]
% Distribute flux over all lit "pixels" 
flux_per_Pupil_pix = flux_in_beam / litPupilPix;    % [ph/s/pix]

% Take sqrt() of flux per pixel to get E-field amplitude per pixel 
    % REMEMBER THAT THE E-field is the sqrt() of the flux!!
EPup_amp = sqrt(flux_per_Pupil_pix); % [sqrt(ph/s/pix)]


% Apply to PUPIL variable so that E-field is scaled correctly
    % NOTE: this doesn't affect the calls to logical() in some functions
    % (eg. genZern_...) since logical() will convert any non-0 value to a 1
PUPIL = PUPIL * EPup_amp;

fprintf('Flux per pupil pix: %f [ph/s/pix]\n',flux_per_Pupil_pix);
fprintf('E-field Amp at Pup: %f [sqrt(ph/s/pix)]\n',EPup_amp);

%%

% For proper airy function:
%  - centroiding accuracy
%  
%  [FWHM (PSF width = lam/D)] / c * SNR        :: c ~ 2 ::
%  
%  ---> To centroid to within 1/10 lam/D, you need an SNR of 10

% Can use this double check that things are working well by working with
% the Airy pattern and making sure that my code matches that equation
% roughly before moving to donut PSF.



% TODO list:
% - finish working out flux stuff
% - Check things are working with Airy PSF:
%   - Create a time-series of PSF's with T/T offsets
%   - Centroid on that PSF and back out where the centroiding algorithm
%       thinks the PSF is
%   - Plot the two (this will let me check how well the centroiding
%       algorithm is detecting the T/T offsets
%   - Increase the noise in the image and repeat. 
%   - See if the results agree with the equation Kent provided above.
%
% ----->>> FIGURE OUT HOW TO DISTRIBUTE FLUX BY WAVELENGTH (fracBW)
%         I think to do this, I can just divide PUPIL in the Epup(:,:,ch) = ...
%         definition line below by the number of wavelengths being sampled.
%         This would effectively distribute the power among all the sampled
%         wavelengths equally. To simulate a non-constant amount of power 
%         (ie. a non-flat spectrum), can divide by a value proportional to
%         the amount of power present in the given spectral channel. One
%         caveat would be to confirm that all requested wavelengths lie
%         within the CRED2 sensitivity window. 
%      ** In fact, this brings up a clever way to deal with the wavelengths. 
%         Instead of doing the flux allocation by computing the relative
%         CRED2 window compared to the astronomical band, I could provide
%         a vector mapping wavelength to power. In the simple limiting
%         case currently implemented, the power in the astronomical band
%         would be divided by the number of wavelengths and then the vector
%         would be all 1's for the region where the CRED2 is sensitive and
%         0's elsewhere. Then the QE could also be a similar vector but
%         actually wavelength-dependent instead of just a flat 85% across
%         the whole window. 
%      ** The only thing to figure out now is how to correctly combine the
%         individual wavelengths into a single image. Would it be a direct 
%         sum of the complex-valued PSF or a mean or...? I think it would
%         be a direct sum do to superposition and the fact that the
%         individual wavelengths are incoherent.

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

%% Add tilts
% Compute lambda0OverD_rad for Keck
Keck_D = 10.949; % [m] circumscribed diameter of Keck pupil
lambda0OverKeckD_rad = lambda0/Keck_D;
% Now convert to mas
lambda0OverKeckD_mas = lambda0OverKeckD_rad*206265*1e3;    


%%%%%% TODO:::::::::::::: Check if this should be /N or /(2*apRad)
% Get tilt phase at central wavelength
phz = 2*pi*tilts(1)/lambda0OverKeckD_mas*coordsPP.X/(2*apRad);%N;
phz = phz + 2*pi*tilts(2)/lambda0OverKeckD_mas*coordsPP.Y/(2*apRad);%N;

figure(10)
for ch = 1:numWavelengths
    Epup(:,:,ch) = Epup(:,:,ch).*exp(1i*phz*lambda0/lambdas(ch));
    
    subplot(1,numWavelengths,ch);
    imagesc(xvalsPP,yvalsPP,angle(Epup(:,:,ch)));
    axis image; 
    axis([-1 1 -1 1]);
    title(['Phase With Tilt at ',num2str(lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(parula(256));
end


%% Get PSF without vortex mask
[~, PSF] = getPSF_mft(Epup, lambdas, foc, coordsPP, coordsFP);

% Compute broadband image intensity (without peak-normalization)
iPSF_BB = sum(abs(PSF).^2,3);
peakFlux = max(iPSF_BB(:));
fprintf('Peak Intensity in PSF: %f [ph/s]\n',peakFlux);

figure(3)
imagesc(xvalsFP,yvalsFP,iPSF_BB);
axis image; 
%axis([-3 3 -3 3]);
title('Broadband PSF w/o vortex');
cb = colorbar;
cb.Label.String = 'Intensity [ph/s]';
%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Make vortex mask
EPM = generateVortexMask( charge, coordsPP, [offsetX offsetY] );

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
[~, PSFv] = getPSF_mft(Epup.*EPM, lambdas, foc, coordsPP, coordsFP);

% Compute broadband image intensity (without peak-normalization)
iPSFv_BB = mean(abs(PSFv).^2,3);

figure(5)
imagesc(xvalsFP,yvalsFP,iPSFv_BB);
axis image; 
%axis([-3 3 -3 3]);
title('broadband PSF w/ vortex');
cb = colorbar;
cb.Label.String = 'Intensity [ph/s]';
%caxis([-3 0])
colormap(parula(256));
drawnow;

%% simulate the image with noise
% 1. get photons per pixel in the pupil plane
    % photons/second/m2/time_sample
% 2. multiply by system throughput/transmission
% 3. take sqrt() and treat that as the magnitude of the electric field

%-- 0) Square the PSF (E-field) to get intensity 
iPSFv = abs(PSFv).^2;   % [photons/s/pixel]

%-- 4) Account for integration time 
Tint = 0.028998518;
iPSFv = Tint * iPSFv; % [ph/pix/t_sample]

%-- 1) Given intensity in focal plane, add photon noise (in photons/s/pix)
% draw from poisson distribution to get shot noise (on per-pixel basis)
noisy_im = poissrnd(iPSFv);     % [photons/t_sample/pix]
    
    % TODO::: CONFIRM THAT THIS IS INDEED THE NOISY IMAGE AND NOT TGHE
    % NOISE ITSELF

%-- 2) Convert from photons/s/pix to electrons/s/pix using QE
CRED2_QE = 0.85;
noisy_im = CRED2_QE * noisy_im;     % [e-/t_sample/pixel]

%-- 3) Add Dark current (in electrons/t/pix)
CRED2_DarkCurrent = 449 * Tint;    % [e-/t_sample/pix] (value @-40C)
% Sample poisson distribution to get dark current per pixel
dark_noise = poissrnd(CRED2_DarkCurrent,size(noisy_im));
% Add noise
noisy_im = dark_noise + noisy_im;   % [e-/s/pix]

%-- 5) Add Read noise (in e-/pix)
CRED2_ReadNoiseRMS = 39.6;  % [e- rms (/pix?)]
% Sample gaussian distribution to get read noise per pixel
   % First 0 is to set noise to be mean-0
read_noise = normrnd(0, CRED2_ReadNoiseRMS, size(noisy_im));
noisy_im = read_noise + noisy_im;

%-- 6) Convert to DN
CRED2_ADU = 2.26;   % [e-/count]
noisy_im = noisy_im / CRED2_ADU;  % [count/pix/t_sample]



% Display the noisy image!
figure(100);
imagesc(noisy_im);
axis image;
cb = colorbar;
% Add marker showing where tilt ~should~ put the PSF centroid
xt = tilts(1)/(lambda0OverKeckD_mas/lambda0Fnum_samp) + Nxi/2 + 1;
yt = tilts(2)/(lambda0OverKeckD_mas/lambda0Fnum_samp) + Nxi/2 + 1;
hold on; scatter(xt, yt, 40, 'r+'); hold off;
xlabel('[pix]'); ylabel('[pix]');
title('Simulated Image with Noise')
ylabel(cb, 'ADU');


% OPTIONAl: Display some sample images for comparison
im_93_bkgd = fitsread('/media/Data_Drive/VFN/VFN_OnSkyData/220413/TCReferences/bias_100459_0032.00_00.02900_01_-40.00000_000_000_000_000.fits');
im_93_clean = mean(fitsread('/media/Data_Drive/VFN/VFN_OnSkyData/220413/CRED2/cred2_0093.fits'),3);
im_93 = im_93_clean - im_93_bkgd;
im_93 = im_93(2:end,:);  % remove tag row
figure(1000); imagesc(im_93)
xlabel('[pix]'); ylabel('[pix]');
title('Sample Real Image #93')
max(im_93(:))
axis image
axis([330, 420, 195, 285])
cb = colorbar();
ylabel(cb, 'ADU')

im_63_bkgd = fitsread('/media/Data_Drive/VFN/VFN_OnSkyData/220413/TCReferences/bias_100459_0032.00_00.02900_01_-40.00000_000_000_000_000.fits');
im_63_clean = mean(fitsread('/media/Data_Drive/VFN/VFN_OnSkyData/220413/CRED2/cred2_0063.fits'),3);
im_63 = im_63_clean - im_63_bkgd;
im_63 = im_63(2:end,:);  % remove tag row
figure(10001); imagesc(im_63)
xlabel('[pix]'); ylabel('[pix]');
title('Sample Real Image #63')
max(im_63(:))
axis image
axis([328, 418, 196, 286])
cb = colorbar();
ylabel(cb, 'ADU')


%%

% 
% %% Working with Poisson 
% %    -- Photon noise: with individual per pixel based on flux per pix
% %               (photons/time_sample)
% %    -- Dark current: same value for all pix but detector specific
% %               (electrons/time_sample)
% 
% % photons --> electrons via QE:  photons*QE = electron
% % Conversion to DN: computed via a value which is electron/DN
% 
% % Define Poisson distribution parameter
% poiss_pam = 1;
% 
% % Create image-size matrix of poisson samples
% poiss_noise = poissrnd(poiss_pam,Nxi,Nxi);
% figure();
% imagesc(poiss_noise);
% axis image;
% colormap(gray(256));
% colorbar;
% title('Image of Poisson Noise')
% 
% % Display the poisson sampling
% % OPTION 1: generate new vector of poisson sampling
% % poiss_samps = 1000;
% % poiss_noise_vec = poissrnd(poiss_pam,poiss_samps);
% % OPTION 2: display the values from the previous image
% poiss_noise_vec = poiss_noise(:);
% % Display
% figure();
% histogram(poiss_noise_vec);
% title('Poisson Sampling')
% xlabel('Value')
% ylabel('Occurrence')
% 
% 
% %-- The poiss_pam value for EACH pixel will be set by the number of counts
% %   on that pixel. ie. resample for each pixel.
% 
% %% Working with Gaussian
% %    -- Read noise: same value for all pix but detector specific
% %           (RMS value of electrons) - if mean-zero then RMS=sigma
% %           added as a value at the end regardless of integration time
% 
% % Define the Gaussian distribution parameters
% gauss_pam_mean  = 1;
% gauss_pam_sigma = 0.2;
% 
% % Create image-size matrix of poisson samples
% gauss_noise = normrnd(gauss_pam_mean,gauss_pam_sigma,Nxi,Nxi);
% figure();
% imagesc(gauss_noise);
% axis image;
% colormap(gray(256));
% colorbar;
% title('Image of GaussianNoise')
% 
% % Display the gaussian sampling;
% % OPTION 1: generate new vector of gaussian sampling
% % gauss_samps = poiss_samps;
% % gauss_noise_vec = normrnd(gauss_pam_mean, gauss_pam_sigma, gauss_samps);
% % OPTION 2: display the values from the previous image
% gauss_noise_vec = gauss_noise(:);
% % Display
% figure(); 
% histogram(gauss_noise_vec);
% title('Gaussian Sampling')
% xlabel('Value')
% ylabel('Occurrence')
% 
% %% Photon counting procedure:
% % 1. get photons per pixel in the pupil plane
%     % photons/second/m2/time_sample
% % 2. multiply by system throughput/transmission
% % 3. take sqrt() and treat that as the magnitude of the electric field
% 
% %% Make a noisy image
% 
% % Rescale the PSF to be some reasonable number of peak counts
% psf_count_scaling = 10000;
% iPSFv_BB = psf_count_scaling*iPSFv_BB;
% 
% % Add the Poisson noise
% noisy_PSFv_poiss = iPSFv_BB + poiss_noise;
% 
% % Display image with poisson noise
% figure(); 
% imagesc(noisy_PSFv_poiss);
% axis image;
% colormap(gray(256));
% colorbar;
% title('PSF Image with Poisson Noise')