%{
    Script to simulate defects of specific quantity and size in a pupil-plane VFN system
    * This script tests a single parmaeter pair. For 2D grid, use the _Scattering_GridSearch script

    This introduces a specified number of defects with specified size into the pupil plane and then
    computes the null and peak planet coupling. It can be used to test the expected performance 
    from a specific mask given those parameters. This is helpful for testing the predicted performance
    of a BeamCo mask if they give us one with defects.
%}

clear; close all; 
addpath(genpath(fullfile('..','VFNlib')));
addpath(genpath(fullfile('..','..','falco-matlab')))

%% Input parameters 

%-- Scattering element parameters
defects_number = 100;   % Total number of defects
defects_maxSize = 7;   % Size along an edge (in pupil samples) of each square defect

%-- Provide regular parameters
% Define smapling info
N = 2^10; % Size of computational grid (NxN samples) 
apRad = N/2-4; % Aperture radius in samples 

%-- Define wavelength info
lambda0 = 2.2e-6; %central wavelength
fracBW = 1e-3; %0.2; % \Delta\lambda/\lambda
numWavelengths = 1; % number of discrete wavelengths 
lambdas = getWavelengthVec(lambda0,fracBW,numWavelengths);% array of wavelengths (meters)

%-- Define charge of the vortex mask at each wavelength
charge = 1*ones(1,numWavelengths); % achromatic

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

%% Get PSF with normal vortex mask
[iPSFv_BB, PSFvG] = getPSF_mft(Epup.*EPM, lambdas, foc, coordsPP, coordsFP);


figure(5)
imagesc(xvalsFP,yvalsFP,iPSFv_BB);
axis image; 
axis([-3 3 -3 3]);
title('broadband PSF w/ vortex');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;


%% Generate defects

%-- Generate scattering elements
% Get indices of all pixels inside the pupil area
indsInside = find(PUPIL >0.5 );

% Randomly select defects_number indices at which to place defects
picks = rand(defects_number,1);
picks = picks * length(indsInside); % rescale from 0-1 to the size of the indices array
picks = floor(picks);   % Round to integers so that we can use them for indexing
picks = max(picks,1);   % Makse sure no values are 0 so we can use for indexing
defects_locations = indsInside(picks);
%Convert to 2D indexing rather than linear
[defects_Rloc, defects_Cloc] = ind2sub(size(PUPIL), defects_locations);

%-- OPTION 1: 
% Randomly set sizes of scattering elements up to a defects_maxSize
%sizes = rand(defects_number,1);
%sizes = sizes * defects_maxSize;    % Scale to maximum size
%sizes = floor(sizes/2); % /2 to set size as radius and then floor to get integers for indexing
%-- OPTION 2:
% Set constant size for all defects
sizes = floor(defects_maxSize*ones(defects_number,1)/2);   % /2 to set size as radius then floot to get integers for indexing

%-- Get start and end indices in both dimensions for each defect
defects_coords = [defects_Rloc-sizes, ...
                  defects_Rloc+sizes, ...
                  defects_Cloc-sizes, ... 
                  defects_Cloc+sizes];
% Make sure no indices are out of range
defects_coords = max(defects_coords,1);
defects_coords = min(defects_coords,N);

%-- Randomly assign certain defects to be opaque (1) vs. fully transmissive (2)
defects_opacity = randi([1,2], defects_number,1);

%-- Get sorted indices so that we do the opaque defects last and hence when
%there are overlaps, the opaque defects will block the phase defects in
%overlap region.
[~, sort_inds] = sort(defects_opacity, 'descend');

%-- Generate mask identifying defects and type of defect
defects_mask = zeros(N);
for dind = 1:defects_number
    sind = sort_inds(dind);
    defects_mask(defects_coords(sind,1):defects_coords(sind,2), ...
                 defects_coords(sind,3):defects_coords(sind,4)) ...
                 = defects_opacity(sind);
end

%-- Generate the amplitude mask
defects_mask_Amp = PUPIL;
defects_mask_Amp(defects_mask == 1) = 0;     % opaque defects are blocked in the pupil

%-- Generate the amplitude mask 
defects_mask_Phz = EPM;
% Replicate defects_mask at each wavelength 
defects_mask = repmat(defects_mask,1,1,numWavelengths);
% transmissive defects, set the phase to 0
defects_mask_Phz(defects_mask == 2) = 1;

%-- Show the defects in phase
figure(4)
for ch = 1:numWavelengths
    subplot(1,numWavelengths, ch);
    imagesc(xvalsPP,yvalsPP,angle(defects_mask_Phz(:,:,ch).*PUPIL));    % use PUPIL to show region inside pupil without including amplitude defects
    axis image; 
    axis([-1 1 -1 1]);
    title([num2str(lambdas(ch)*1e9) ' nm']);
    colormap(hsv(256));
    colorbar; 
end
sgtitle('Phase in Pupil with Defects')
drawnow;
%-- Show the defects in amplitude
figure(5)
for ch = 1:numWavelengths
    subplot(1,numWavelengths, ch);
    imagesc(xvalsPP,yvalsPP,abs(defects_mask_Phz(:,:,ch).*defects_mask_Amp));
    axis image; 
    axis([-1 1 -1 1]);
    title([num2str(lambdas(ch)*1e9) ' nm']);
    colormap(parula(256));
    colorbar; 
end
sgtitle('Amplitude in Pupil with Defects')
drawnow;

%-- Report some statistics on the defects
defects_physicalSize = defects_maxSize * DPup/(2*apRad);
fprintf('Defect Physical Size: %0.2f um\n', defects_physicalSize*1e6) % *1e6 to get to um
fprintf('Defect Physical Area: %0.2f um2\n', (defects_physicalSize*1e6)^2);
fprintf('Number of Defects in CA: %d\n', defects_number);
fprintf('    -- %d amplitude defects\n', sum(defects_opacity == 1));
fprintf('    -- %d phase defects\n', sum(defects_opacity == 2));
defects_density = defects_number / (pi * (DPup/2)^2);
fprintf('Defect Density: %0.2f defects/mm2\n', defects_density/1e6)   % /1e6 to get to mm2

%% Get PSF with defective vortex mask
[iPSFv_BB, PSFvD] = getPSF_mft(defects_mask_Amp.*defects_mask_Phz.*Epup, lambdas, foc, coordsPP, coordsFP);


figure(6)
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
    eta_onAx(ch) = (abs(sum(sum(PSFvD(:,:,ch).*fibmode(:,:,ch)))).^2)/totalPower0;
end

%-- Compute 2D coupling map
eta_maps = zeros(Nxi,Nxi,numWavelengths);
for ch = 1:numWavelengths
    % Compute the monochromatic coupling map (spatial units of lambda/D)
    eta_maps(:,:,ch) = generateCouplingMap( fibmode(:,:,ch), PSFvD(:,:,ch), totalPower0, 5*lambda0Fnum_samp, coordsFP);

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
    figure(8);
    imagesc(xvalsFP, yvalsFP, eta_map_BB);
    axis image;
    axis([-3 3 -3 3]);
    title('Broadband \eta');
    colorbar;
    colormap(gray(256));
    
    % Plot Log Scale
    figure(9);
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