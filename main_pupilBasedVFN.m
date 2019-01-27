clear; %close all; 
addpath('VFNlib');

N = 2^11; % Size of computational grid (NxN samples) 
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

offsetX = 0;
offsetY = 0; 

EPM = generateVortexMask( charge, coords, [offsetX offsetY] );

figure(3)
imagesc(xvals/apRad,yvals/apRad,angle(EPM.*EP));
axis image; 
axis([-1 1 -1 1]);
title('Pupil phase');
colormap(hsv(256));
colorbar; 

%% Add Zernike aberration 

noll_index = 2; % Noll index
coeff = 0.5; % waves rms
[Z,n,m] = generateZernike(noll_index,apRad,coords.RHO,coords.THETA);
Z = Z/sqrt(mean(Z(logical(EP)).^2)); % Re-normalize (useful when pupil is not a circle)
ABER = exp(1i*2*pi*coeff*Z);

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

fiberDiam = 1.4; % units of lambda/D

fibermode0 = generateFiberMode(fiberDiam*lambdaOverD,coords);

coupling_eff_map = generateCouplingMap( fibermode0, PSFv, totalPower0, 3*lambdaOverD);

figure(5)
imagesc(xvals/lambdaOverD,yvals/lambdaOverD,coupling_eff_map);
axis image; 
axis([-3 3 -3 3]);
title('Coupling map');
colorbar; 
colormap(parula(256));

%% Get throughput curve for different F numbers 

if(produceThptCurve)
    angSeps = (0:0.1:2)*lambdaOverD;
    fiberDiams = [1 1.4 1.8];
    figure(5); 
    for fiberDiam = fiberDiams

%         fibermode0 = sqrt(2/(pi*(fiberDiam*lambdaOverD/2)^2))* ...
%             exp(-(RHO/(fiberDiam*lambdaOverD/2)).^2);
        fibermode0 = generateFiberMode(fiberDiam*lambdaOverD,coords);

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