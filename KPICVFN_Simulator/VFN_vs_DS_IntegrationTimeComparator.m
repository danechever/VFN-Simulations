%{
Script to generate plot showing where VFN provides integration-time reduction over DS-mode.

** This assumes integration time reduction is purely due to eta_s/eta_p^2. It uses simulated performance
    values. For a similar script that instead uses as-measured on-sky curves from KPIC VFN and also computes
    SNR's, SNR ratios, etc., look on VFN server under paper_results/ directory where there is a python script
    that does these computations from the measured throughputs rather than E2E simulated coupling performance.
%}

clear; close all; 
addpath(genpath(fullfile('..','VFNlib')));
addpath(genpath(fullfile('..','..','falco-matlab')))

%% Input parameters 
%-- Define spatial sampling in pupil plane
N = 2^10;

%-- Define wavelength simulation info
lambda0 = 2.2e-6;

%-- Vortex charge
charge = 2;     

%-- Define wavefront error at the central wavelength
  % 1: Pist, 2: Tip, 3: Tilt, defoc, ob.astig, ver.astig, ver.coma, hor.coma,
  % 9: ver.tref, ob.tref, spher
nolls = 4:8;    % Noll indices to simulate
coeffs = 0*[0.0, 0.0, 0.05, 0.05, 0.05];  % amplitude in [Waves RMS]

% Since using falco MFT propogator, can choose our final image sampling
lambda0Fnum_samp = 50; %[samples/ (lam0 F#)] in units of samples in the image plane
im_size = 10;    %[lam0/D] field of view in final image plane

%% Finish setting up constants/workspace using input values
%-- Spatial sampling
apRad = N/2-4; % Aperture radius in samples 

%-- Parameters for SMF (Thorlabs SM2000 in this case)
%     % link: https://www.thorlabs.com/NewGroupPage9_PF.cfm?ObjectGroup_ID=949
fiber_props.core_rad = 11e-6/2;% Core radius [um]
fiber_props.n_core = 1.4436;% core index 
fiber_props.n_clad = 1.4381;% cladding index 
fiber_props.type = 'gaussian';  % Alternate option 'bessel' works and is slightly higher fidelity but will take longer to run

%-- Define parameters for falco MFT propagator
    % these don't matter in themselves as long as they are consistent w/ each 
    % other and with lambda
% OPTION 1: define focal length and pup diameter manually
%foc = 11e-3;      %[m] final focal length
DPup = 12.5e-3;    %[m] pupil size 
% OPTION 2: solve for focal length based on ideal Fnum and pup diameter
Fnum = getMFD(fiber_props,lambda0)/(lambda0*1.42); % focal ratio of the beam at the fiber
foc = Fnum*DPup;
fprintf('Focus in use:    %f [mm]\n',foc*1e3)
fprintf('F/# in use:      %f \n',Fnum)
fprintf('Pup Diam in use: %f [mm]\n',DPup*1e3)

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

%% Define pupil field phase
phz = generateZernike_fromList( nolls, coeffs, PUPIL, apRad, coordsPP); 

figure(2);
Epup = exp(1i*phz).*PUPIL;
imagesc(xvalsPP,yvalsPP,angle(Epup));
axis image; 
axis([-1 1 -1 1]);
title(['Phase at ',num2str(lambda0*1e9),'nm']);
colorbar; 
colormap(parula(256));
   
drawnow;

%% Get PSF without vortex mask
[iPSF_BB, PSF] = getPSF_mft(Epup, lambda0, foc, coordsPP, coordsFP);

figure(3)
imagesc(xvalsFP,yvalsFP,iPSF_BB);
axis image; 
%axis([-3 3 -3 3]);
title('PSF w/o vortex');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Make vortex mask
EPM = generateVortexMask( charge, coordsPP, [0 0] );
fprintf('Generated a theoretically ideal vortex\n')

figure(4);
% Show phase
imagesc(xvalsPP,yvalsPP,angle(EPM.*Epup));
axis image; 
axis([-1 1 -1 1]);
title(['Phase at ',num2str(lambda0*1e9),'nm']);
colormap(gca(), hsv(256));
colorbar; 

drawnow;

%% Get PSF with vortex mask
[iPSFv_BB, PSFv] = getPSF_mft(Epup.*EPM, lambda0, foc, coordsPP, coordsFP);


figure(5)
imagesc(xvalsFP,yvalsFP,iPSFv_BB);
axis image; 
axis([-3 3 -3 3]);
title('PSF w/ vortex');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Generate fibermode at each lambda     
% Generate the fiber mode for this wavelength with proper scaling
fibmode = generateSMFmode_mft( fiber_props, lambda0, coordsFP.dx, coordsFP);

%% Calculate Throughputs (DS-mode)
%-- Get on-axis coupling (using overlap integral)
eta_onAx_D = (abs(sum(sum(PSF.*fibmode))).^2)/totalPower0;

%-- Compute 2D coupling map
eta_map_D = generateCouplingMap( fibmode, PSF, totalPower0, 5*lambda0Fnum_samp, coordsFP);

%-- Compute azimuthally-averaged performance
% Compute average
[eta_D, sep_D] = radialAverage(eta_map_D, [Nxi/2+1 Nxi/2+1]);
sep_D = sep_D/lambda0Fnum_samp; % [lam/D] (converted from pix to lam/D directly)

disp('Key Coupling Points (DS):')
fprintf('lambda = %0.1f nm,    on-axis coup = %0.2e\n',lambda0*1e9, eta_onAx_D);
    
%% Calculate Throughputs (VFN-mode)
%-- Get null depth (using overlap integral)
eta_onAx_V = (abs(sum(sum(PSFv.*fibmode))).^2)/totalPower0;

%-- Compute 2D coupling map
eta_map_V = generateCouplingMap( fibmode, PSFv, totalPower0, 5*lambda0Fnum_samp, coordsFP);

%-- Compute azimuthally-averaged performance
% Compute average
[eta_V, sep_V] = radialAverage(eta_map_V, [Nxi/2+1 Nxi/2+1]);
sep_V = sep_V/lambda0Fnum_samp; % [lam/D] (converted from pix to lam/D directly)

disp('Key Coupling Points (VFN):')
% Find peak (value and separation at which it occurs) in radial average
[peak, pind] = max(eta_V,[],'all', 'linear');
fprintf('lambda = %0.1f nm,    on-axis coup = %0.2e,    peak = %0.2f %% (at %0.2f lam0/D)\n',lambda0*1e9, eta_onAx_V, peak*100, sep_V(pind));

%% Compute integration time reduction (Use eta_s/eta_p^s formula)
% For DS-mode, assume we align w/ planet perfectly s.t. eta_p=max coup
% in such a scenario, eta_s is given by coupling vector directly
tau_D = eta_D / (eta_D(1).^2);
% ALTERNATIVE: Provide on-axis planet coupling
%tau_D = eta_D / (0.3.^2);

% For VFN-mode, assume we align w/ star perfectly s.t. eta_s=on-axis value
% in such a scenario, eta_p is given by coupling vector directly
tau_V = eta_V(1) ./ (eta_V.^2);
% ALTERNATIVE: Provide on-axis VFN null
%tau_V = 1e-3 ./ (eta_V.^2);

%% Display Plots
lw = 2; % linewidth for plots
Dcolor = "#0072BD";
Vcolor = "#D95319";

% Plot of coupling vs. separation
figure(10);
plot(sep_D, eta_D, 'Color', Dcolor, 'linewidth', lw);
hold on;
plot(sep_V, eta_V, 'Color', Vcolor, 'linewidth', lw);
hold off;
legend('DS', 'VFN')
xlabel('Separation [\lambda/D]')
ylabel('Coupling [fractional]')
title('Coupling Efficiency')

% Plot of integration time reduction vs. separation
figure(11);
semilogy(sep_D, tau_D, 'Color', Dcolor, 'linewidth', lw);
hold on;
semilogy(sep_V, tau_V, 'Color', Vcolor, 'linewidth', lw);
hold off;
legend('DS', 'VFN')
ylim([0, 1])    % limit to only show when system is improving performance
xlabel('Separation [\lambda/D]')
ylabel('\tau=\eta_s/(\eta_p^2) [fractional]')
title('Integration Time Reduction')

% Plot VFN 2D map
% Coup map linear scale
figure(12);
subplot(1, 2, 1);
imagesc(xvalsFP,yvalsFP,eta_map_V);
axis image; 
axis([-3 3 -3 3]);
xlabel('\lambda_0/D')
ylabel('\lambda_0/D')
colorbar; 
colormap(gray(256));
% Coup map log scale
subplot(1, 2, 2);
imagesc(xvalsFP,yvalsFP,eta_map_V);
ax = gca();
ax.ColorScale='Log';
axis image; 
axis([-3 3 -3 3]);
xlabel('\lambda_0/D')
ylabel('\lambda_0/D')
colorbar; 
colormap(gray(256));
% figure format
sgtitle('VFN-Mode Coupling');

drawnow

% Plot DS 2D map
% Coup map linear scale
figure(13);
subplot(1, 2, 1);
imagesc(xvalsFP,yvalsFP,eta_map_D);
axis image; 
axis([-3 3 -3 3]);
xlabel('\lambda_0/D')
ylabel('\lambda_0/D')
colorbar; 
colormap(gray(256));
% Coup map log scale
subplot(1, 2, 2);
imagesc(xvalsFP,yvalsFP,eta_map_D);
ax = gca();
ax.ColorScale='Log';
axis image; 
axis([-3 3 -3 3]);
xlabel('\lambda_0/D')
ylabel('\lambda_0/D')
colorbar; 
colormap(gray(256));
% figure format
sgtitle('DS-Mode Coupling');

drawnow