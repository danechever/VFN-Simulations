clear; %close all; 
addpath('VFNlib');

N = 2^12; % Size of computational grid (NxN samples) 
apRad = 150; % Aperture radius in samples 
charge = 1; % Charge of the vortex mask 

coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
xvals = coords.xvals;% Helpful for plotting
yvals = coords.yvals;

produceThptCurve = false; % False to skip the time consuming throughput curve calculation at the end
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

offsetX = 0;%0.095;%1.6147*apRad; % <-- PSF % 0.0952*apRad; % <-- best null position % 
offsetY = 0;%0.052;%0.0011*apRad; % <-- PSF % 0.0524*apRad; % <-- best null position % 

EPM = generateVortexMask( charge, coords, [offsetX offsetY] );

figure(3)
imagesc(xvals/apRad,yvals/apRad,angle(EPM.*EP));
axis image; 
axis([-1 1 -1 1]);
title('Pupil phase');
colormap(hsv(256));
colorbar; 

%% Add Zernike aberration 

noll_index = 7; % Noll index (Coma - vertical 90*)
coeff7 = -0.0006; %<--VortOut % -0.0008; %<--VortIn % % waves rms
[Z,n,m] = generateZernike(noll_index,apRad,coords.RHO,coords.THETA);
Z = Z/sqrt(mean(Z(logical(EP)).^2)); % Re-normalize (useful when pupil is not a circle)
ABER = exp(1i*2*pi*coeff7*Z);

noll_index = 8; % Noll index (Coma - horizontal 0*)
coeff8 = 0.0030; %<--VortOut % 0.0033; %<--VortIn % % waves rms
[Z,n,m] = generateZernike(noll_index,apRad,coords.RHO,coords.THETA);
Z = Z/sqrt(mean(Z(logical(EP)).^2)); % Re-normalize (useful when pupil is not a circle)
ABER = ABER.*exp(1i*2*pi*coeff8*Z);

noll_index = 5; % Noll index (Astig - Oblique 45*) 
coeff5 =  -0.0104; %<--VortOut % -0.0099; %<--VortIn % % waves rms
[Z,n,m] = generateZernike(noll_index,apRad,coords.RHO,coords.THETA);
Z = Z/sqrt(mean(Z(logical(EP)).^2)); % Re-normalize (useful when pupil is not a circle)
ABER = ABER.*exp(1i*2*pi*coeff5*Z);

noll_index = 6; % Noll index (Astig - Vertical 0*) 
coeff6 =  -0.0095; %<--VortOut % -0.0095; %<--VortIn % %waves rms
[Z,n,m] = generateZernike(noll_index,apRad,coords.RHO,coords.THETA);
Z = Z/sqrt(mean(Z(logical(EP)).^2)); % Re-normalize (useful when pupil is not a circle)
ABER = ABER.*exp(1i*2*pi*coeff6*Z);

noll_index = 11; % Noll index (Primary Spherical) 
coeff11 = 0;% -0.0194; % waves rms
[Z,n,m] = generateZernike(noll_index,apRad,coords.RHO,coords.THETA);
Z = Z/sqrt(mean(Z(logical(EP)).^2)); % Re-normalize (useful when pupil is not a circle)
ABER = ABER.*exp(1i*2*pi*coeff11*Z);

noll_index = 4; % Noll index (Focus) 
coeff4 = 0;%0.01; % waves rms
[Z,n,m] = generateZernike(noll_index,apRad,coords.RHO,coords.THETA);
Z = Z/sqrt(mean(Z(logical(EP)).^2)); % Re-normalize (useful when pupil is not a circle)
ABER = ABER.*exp(1i*2*pi*coeff4*Z);

figure(3)
imagesc(xvals/apRad,yvals/apRad,angle(EPM.*EP.*ABER));
axis image; 
axis([-1 1 -1 1]);
title('Pupil phase');
colormap(hsv(256));
colorbar; 

%% Get PSF with vortex mask

PSFv = myfft2(EP.*EPM.*ABER); % PSF with vortex mask and aberration 

figure(4)
imagesc(xvals/lambdaOverD,yvals/lambdaOverD,abs(PSFv).^2/normI);
axis image; 
axis([-3 3 -3 3]);
title('log10(PSF) w/ vortex');
colorbar; 
colormap(parula(256));


%% Coupling map
% denominator = 1.4*lambda*F# <-- ideal Df for that F#
fiberDiam = (4.64/4.64)*1.4; % units of lambda/D  

% Parameters for Thorlabs SM600
    % link: https://www.thorlabs.com/NewGroupPage9_PF.cfm?ObjectGroup_ID=949
fiber_props.type = 'gaussian';
fiber_props.core_rad = 4.3e-6/2;% Core radius [um]
lambda = 635e-9;% wavelength [um]
fiber_props.n_core = 1.4606;% core index (interpolated from linear fit to 3 points)
fiber_props.n_clad = 1.4571;% cladding index (interpolated from linear fit to 3 points)
Fnum = 5.22; % optimal focal ratio

%fibermode0 = generateSMFmode( n_core, n_clad, core_rad, lambda, lambda*Fnum/lambdaOverD, coords );
fibermode0 = generateSMFmode_gaussian(fiberDiam*lambdaOverD,coords);
%fibermode0 = generateSMFmode( fiber_props, lambda, Fnum, lambdaOverD, coords);

abs(sum(sum(fibermode0.*PSFv))).^2/totalPower0
%coupling_eff_map = generateCouplingMap( fibermode0, PSFv, totalPower0, 10*lambdaOverD, coords);

figure(5)
imagesc(xvals/lambdaOverD,yvals/lambdaOverD,coupling_eff_map);
axis image; 
axis([-3 3 -3 3]);
title('Coupling map');
colorbar; 
colormap(parula(256));

%% Save fits of coupling map
% % Import library for managing fits headers
% import matlab.io.*
% % Create filename with path
% fldnm = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\MonoWrapUp\NewAnalysis\';
% flnm  = sprintf('Sim_Gauss_PSF_DF%0.2f_FN%0.2f_WFEOn.fits', fiberDiam, Fnum);
% % Create fits
% fitmap= fits.createFile([fldnm flnm]);
% fits.createImg(fitmap,'double',size(coupling_eff_map));
% %header data
% fits.writeKey(fitmap, 'apRad   ', apRad);
% fits.writeKey(fitmap, 'VortX   ', offsetX/apRad);
% fits.writeKey(fitmap, 'VortY   ', offsetY/apRad);
% fits.writeKey(fitmap, 'Noll4WF ', coeff4);
% fits.writeKey(fitmap, 'Noll5WF ', coeff5);
% fits.writeKey(fitmap, 'Noll6WF ', coeff6);
% fits.writeKey(fitmap, 'Noll7WF ', coeff7);
% fits.writeKey(fitmap, 'Noll8WF ', coeff8);
% fits.writeKey(fitmap, 'Noll11WF', coeff11);
% fits.writeKey(fitmap, 'fibDiam ', fiberDiam);
% fits.writeKey(fitmap, 'Fnumber ', Fnum);
% % Save
% fits.setCompressionType(fitmap,'NOCOMPRESS');
% fits.writeImg(fitmap,coupling_eff_map);
% fits.closeFile(fitmap);

%% Get Tip/Tilt sensitivity 
% Optional section to calculate null depth given a fiber misplacement
% val = 9e-9/(lambda*Fnum);
% tilt = exp(1i*2*pi*val*lambdaOverD*coords.X/N);
% FPp_t = myfft2(EP.*EPM.*tilt);
% 
% coup = abs(sum(sum(fibermode0.*FPp_t))).^2/totalPower0

%% Get throughput curve for different F numbers 

if(produceThptCurve)
    angSeps = (0:0.1:2)*lambdaOverD;
    fiberDiams = [1 1.4 1.8];
    figure(5); 
    for fiberDiam = fiberDiams

%         fibermode0 = sqrt(2/(pi*(fiberDiam*lambdaOverD/2)^2))* ...
%             exp(-(RHO/(fiberDiam*lambdaOverD/2)).^2);
        fibermode0 = generateSMFmode_gaussian(fiberDiam*lambdaOverD,coords);

        coupling_planet = [];
        for angSep = angSeps
            tilt = exp(1i*2*pi*angSep*coords.X/N);

            FPp = myfft2(EP.*EPM.*tilt);
            coupling_planet = [coupling_planet,abs(sum(sum(fibermode0.*FPp))).^2/totalPower0];
        end

        plot(angSeps/lambdaOverD,coupling_planet);hold on;
        xlabel('Angular separation (\lambda/D)');
        ylabel('\eta_p');
        drawnow;
    end
    hold off;
    legLabels = {};
    for index = 1:numel(fiberDiams)
        legLabels{index} = ['D_f = ',num2str(fiberDiams(index)),' \lambda F#'];
    end
    legend(legLabels);
end