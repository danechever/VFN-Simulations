%{
    Script for simulating a simple phase knife nuller as manufactured by Zeiss for HISPEC
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
lambda0 = 1.55e-6; % [meters] central wavelength
fracBW = 0.2; % \Delta\lambda/\lambda
numWavelengths = 5; % number of discrete wavelengths 
lambdas = getWavelengthVec(lambda0,fracBW,numWavelengths);% array of wavelengths (meters)

%-- Define phase knife properties 
    % To model ideal phase knife:
phzStep = pi * ones(size(lambdas));    % phase shift for knife. Use vector for each wavelength element
    % For true, non-wedged phase knife:
% phzStep = pi * lambda0./lambdas;     % phase shift for knife. Use vector for each wavelength element
clockAng_deg = 0;   % [deg] Clocking angle for the phase knife

%-- Give offsets for the phase knife
offsetX = 0;    % [samples in the pupil-plane]
offsetY = 0; 

%-- Define Pupil type
isKeck = false; % aperture to model. True = Keck, False = Circular

%-- Define wavefront error at the central wavelength
  % 1) Pist, Tip, Tilt, defoc, ob.astig, ver.astig, ver.coma, hor.coma,
  % 9) ver.tref, ob.tref, spher
nolls = 4:8;
coeffs = [0.0, 0.0, 0.0, 0.0, 0.0];  % [Waves RMS]

%-- Parameters for SMF 
    % Check generateSMFmode_mft for what parameters can be provided
    % Params below are for SMF28 from:
    % - https://media.thorlabs.com/globalassets/items/s/sm/smf/smf-28-j9/ttn015884-s01.pdf?v=0116030005
    % - https://focenter.com/media/wysiwyg/docs/Corning-Corning-SMF-28E-Single-Mode-Bare-Fiber-Fiber-Optic-Center.pdf
fiber_props.core_rad = 8.2e-6/2;% [meters] Core radius 
fiber_props.NA = 0.14;
fiber_props.type = 'bessel';  

%-- Define parameters for falco MFT propagator
    % these don't matter in themselves as long as they are consistent w/ each 
    % other and with lambda
% OPTION 1: define focal length and pup diameter manually
% foc = 36.07e-3;     %[m] final focal length
% DPup = 10.76e-3;    %[m] pupil size 
% OPTION 2: solve for pupil diameter based on ideal Fnum and focal length
    % Note: q = 1.33 seems to be optimal for phase knife 
q = 1.33;   % "q" term from Gary's 2019 SPIE paper ( q = MFD / (lam*F#) )
Fnum = getMFD(fiber_props,lambda0)/(lambda0*q); % focal ratio of the beam at the fiber
foc = 36.07e-3;         %[m] final focal length
DPup = foc / Fnum;      %[m] pupil size 
fprintf('Focus in use: %f [mm]\n',foc*1e3)
fprintf('F/# in use:      %f \n', foc/DPup)
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

if isKeck
    PUPIL = makeKeckPupil( 2*apRad, N );
else
    PUPIL = makeCircularPupil( apRad, N );
end

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

%% Get PSF without phase mask
[iPSF_BB, PSF] = getPSF_mft(Epup, lambdas, foc, coordsPP, coordsFP);

figure(3)
imagesc(xvalsFP,yvalsFP,iPSF_BB);
axis image; 
%axis([-3 3 -3 3]);
title('broadband PSF w/o knife');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Make phase mask (using staircase from falco)
EPM = generatePhaseKnifeMask( phzStep, clockAng_deg, coordsPP, [offsetX offsetY] );
    
central_band_index = ceil(numWavelengths/2);

figure(4)
imagesc(xvalsPP,yvalsPP,angle(EPM(:,:,central_band_index).*Epup(:,:,central_band_index)));
axis image; 
axis([-1 1 -1 1]);
title('Pupil phase at \lambda_0');
colormap(parula(256));
colorbar; 
drawnow;

%% Get PSF with phase mask
[iPSFm_BB, PSFm] = getPSF_mft(Epup.*EPM, lambdas, foc, coordsPP, coordsFP);

figure(5)
imagesc(xvalsFP,yvalsFP,iPSFm_BB);
axis image; 
axis([-3 3 -3 3]);
title('broadband PSF w/ knife');
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
    eta_onAx(ch) = (abs(sum(sum(PSFm(:,:,ch).*fibmode(:,:,ch)))).^2)/totalPower0;
end

%-- Compute 2D coupling map
eta_maps = zeros(Nxi,Nxi,numWavelengths);
for ch = 1:numWavelengths
    % Compute the monochromatic coupling map (spatial units of lambda/D)
    eta_maps(:,:,ch) = generateCouplingMap( fibmode(:,:,ch), PSFm(:,:,ch), totalPower0, 6*lambda0Fnum_samp, coordsFP);

end

% --- OLD Vortex planet coupling computation - Left here just for reference
% %-- Compute azimuthally-averaged performance
% for ch = 1:numWavelengths
%     % Compute average
%     [tmpAvg, qvec] = radialAverage(eta_maps(:,:,ch), [Nxi/2+1 Nxi/2+1]);
% 
%     % Yes, I know it's not efficient to recreate matrix on every itr., I
%     % don't care in this particular application...
%     eta_pAvgs(:,ch) = tmpAvg;
% end
% % Since all eta_maps have the same spatial scale, they all have the same qvec
% qvec = qvec/lambda0Fnum_samp; % [lam/D] (converted from pix to lam/D directly)

% Phase knife planet coupling
    % This assumes we only care about maximizing the peak
for ch = 1:numWavelengths
    % Yes, I know it's not efficient to recreate matrix on every itr., I
    % don't care in this particular application...
    [eta_pMaxs(ch), idx] = max(eta_maps(:,:,ch), [], 'all');   % Find global max
    [idx1, idx2] = ind2sub(size(eta_maps(:,:,ch)), idx);        % convert to 2D inds
    eta_pPos(ch) = coordsFP.RHO(idx1, idx2) / lambda0Fnum_samp; % [lam0/D] convert from pix to lam0/D
end

disp('Key Coupling Points:')
for ch = 1:numWavelengths
    fprintf('lambda = %0.2f nm,    on-axis coup = %0.2e,    peak coup = %0.4f %%,  peak pos = %0.3f lambda0/D\n',lambdas(ch)*1e9, eta_onAx(ch), eta_pMaxs(ch)*100, eta_pPos(ch));
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
    
    % Compute peak broadband coupling
    [eta_pMax_BB, idx] = max(eta_map_BB, [], 'all');   % Find global max
    [idx1, idx2] = ind2sub(size(eta_map_BB), idx);     % convert to 2D inds
    eta_pMax_posBB = coordsFP.RHO(idx1, idx2) / lambda0Fnum_samp; % [lam0/D] convert from pix to lam0/D
    
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
    fprintf('Broadband Performance:     on-axis coup = %0.2e,    peak coup = %0.4f %%,   peak pos = %0.3f lambda0/D \n', mean(eta_onAx), eta_pMax_BB*100, eta_pMax_posBB);
end