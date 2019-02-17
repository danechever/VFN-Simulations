clear; %close all; 
addpath('VFNlib');

N = 2^11; % Size of computational grid (NxN samples) 
apRad = 250; % Aperture radius in samples 
charge = 1; % Charge of the vortex mask 

coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
xvals = coords.xvals;% Helpful for plotting
yvals = coords.yvals;

%% Create array with pupil function

EP = makeCircularPupil( apRad, N );

figure(1)
imagesc(xvals/apRad,yvals/apRad,EP);
axis image; 
axis([-1 1 -1 1]);
title('Pupil');
colorbar; 
colormap(parula(256));

%% Get PSF without vortex mask and normalization factors 

% Keeps the origin at the center of the array (keep N even)
myfft2 = @(x) fftshift(fft2(fftshift(x)));
myifft2 = @(x) fftshift(ifft2(fftshift(x)));

PSF = myfft2(EP); % PSF with a flat wavefront 

normI = max(max(abs(PSF).^2));% Normalization for PSF plots
totalPower0 = sum(sum(abs(PSF).^2));% Normalization for coupling fractions
lambdaOverD = N/apRad/2; % lam/D in units of samples in the image plane

figure(2)
imagesc(xvals/lambdaOverD,yvals/lambdaOverD,abs(PSF).^2/normI);
axis image; 
axis([-3 3 -3 3]);
title('PSF w/o vortex');
colorbar; 
colormap(parula(256));

%% Make vortex mask 

EPM = generateVortexMask( charge, coords, [0 0] );

figure(3)
imagesc(xvals/apRad,yvals/apRad,angle(EPM.*EP));
axis image; 
axis([-1 1 -1 1]);
title('Pupil phase');
colormap(hsv(256));
colorbar; 


%% Get PSF with vortex mask

PSFv = myfft2(EP.*EPM); % PSF with vortex mask and aberration 

figure(4)
imagesc(xvals/lambdaOverD,yvals/lambdaOverD,abs(PSFv).^2/normI);
axis image; 
axis([-3 3 -3 3]);
title('log10(PSF) w/ vortex');
colorbar; 
colormap(parula(256));


%% Coupling map for Gaussian mode 

fiberDiam = 1.4; % units of lambda/D

fibermode1 = generateSMFmode_gaussian(fiberDiam*lambdaOverD,coords);
peak_coupling1 = abs(sum(sum(fibermode1.*PSFv))).^2/totalPower0
coupling_eff_map = generateCouplingMap( fibermode1, PSFv, totalPower0, 3*lambdaOverD);

figure(5)
imagesc(xvals/lambdaOverD,yvals/lambdaOverD,fibermode1);
axis image; 
axis([-3 3 -3 3]);
title('Gaussian mode');
colorbar; 
colormap(parula(256));

figure(6)
imagesc(xvals/lambdaOverD,yvals/lambdaOverD,coupling_eff_map);
axis image; 
axis([-3 3 -3 3]);
title('Coupling map for Gaussian');
colorbar; 
colormap(parula(256));

%% Coupling map for true SMF

% Parameters for Thorlabs SM2000
core_rad = 11e-6/2;% Core radius [um]
lambda = 2e-6;% wavelength [um]
n_core = 1.4436;% core index
n_clad = 1.4381;% cladding index
Fnum = 4.58; % optimal focal ratio

fibermode2 = generateSMFmode( n_core, n_clad, core_rad, lambda, lambda*Fnum/lambdaOverD, coords );
peak_coupling2 = abs(sum(sum(fibermode2.*PSFv))).^2/totalPower0
coupling_eff_map = generateCouplingMap( fibermode2, PSFv, totalPower0, 3*lambdaOverD);

figure(7)
imagesc(xvals/lambdaOverD,yvals/lambdaOverD,fibermode2);
axis image; 
axis([-3 3 -3 3]);
title('Fiber mode');
colorbar; 
colormap(parula(256))

figure(8)
imagesc(xvals/lambdaOverD,yvals/lambdaOverD,coupling_eff_map);
axis image; 
axis([-3 3 -3 3]);
title('Coupling map for SMF');
colorbar; 
colormap(parula(256));

