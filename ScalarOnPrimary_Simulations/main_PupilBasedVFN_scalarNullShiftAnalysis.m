%This file recreates the results in Ruane et. al. 2019 for a circular pupil

clear; close all; 
addpath('VFNlib');

%% Input parameters 

% Define smapling info
N = 2^10; % Size of computational grid (NxN samples) 
apRad = 64; % Aperture radius in samples 

% Define wavelength info
lambda0 = 2.2e-6; %central wavelength
fracBW = 0.1818; %\Delta\lambda/\lambda
numWavelengths = 5;% number of discrete wavelengths 
lambdas = getWavelengthVec(lambda0,fracBW,numWavelengths);% array of wavelengths (meters)

% Define charge of the vortex mask at each wavelength
%charge = ones(1,numWavelengths); % achromatic
%chargeC = 1;
charge = 1 * lambda0./lambdas; % simple scalar model

% Define wavefront error at the central wavelength
nolls = 4:8;
coeffs = 0*[0.01,-0.0099,-0.0095,-0.0008,0.0033];

% Give offsets for the vortex mask
offsetX = 0;%0.0952*apRad;
offsetY = 0;%0.0524*apRad; 

%% Generate the coordinate system

coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
xvals = coords.xvals;% Helpful for plotting
yvals = coords.yvals;

%% Create array with pupil function

PUPIL = makeCircularPupil(apRad, N );
[normI, totalPower0] = getNormalization(PUPIL);% Normalization factors
lambdaOverD = N/apRad/2; % lam/D in units of samples in the image plane

 figure(1)
 imagesc(xvals/apRad,yvals/apRad,PUPIL);
 axis image; 
 axis([-1 1 -1 1]);
 title('Pupil');
 colorbar; 
 colormap(parula(256));
 grid on;
 drawnow;

%% Define pupil field

phz = generateZernike_fromList( nolls, coeffs, PUPIL, apRad, coords);
% for ch = 1:numWavelengths
%     phz(:,:,ch) = angle(makeKeckPupilPhase(2*apRad,N,charge(ch)));
% end
%phz = angle(makeKeckPupilPhase(2*apRad,N,chargeC));
% phz2 = angle(makeKeckPupilField(2*apRad,N));

figure(2);
for ch = 1:numWavelengths
    Epup(:,:,ch) = exp(1i*phz*lambda0/lambdas(ch)).*PUPIL;
    
    subplot(1,numWavelengths,ch);
    imagesc(xvals/apRad,yvals/apRad,angle(Epup(:,:,ch)));
    axis image; 
    axis([-1 1 -1 1]);
    title(['Phase at ',num2str(lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(hsv(256));
   
end
drawnow;

%% Get PSF without vortex mask

%Get broadband PSF
iPSF_BB = getPSF(Epup,lambda0,lambdas,normI,coords);

figure(3)
imagesc(xvals/lambdaOverD,yvals/lambdaOverD,iPSF_BB);
axis image; 
axis([-2 2 -2 2]);
title('broadband PSF w/o vortex');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Make vortex mask 

EPM = generateVortexMask( charge, coords, [offsetX offsetY] );

central_band_index = ceil(numWavelengths/2);

figure(4)
for ch = 1:numWavelengths
    subplot(1,numWavelengths,ch);
    imagesc(xvals/apRad,yvals/apRad,angle(EPM(:,:,ch).*Epup(:,:,ch)));
    axis image; 
    axis([-1 1 -1 1]);
    title(['Pupil phase at ', num2str(lambdas(ch)*1e9) 'nm']);
    colorbar;
    colormap(hsv(256));
end
drawnow;

%% Get PSF with vortex mask

iPSFv_BB = getPSF(Epup.*EPM,lambda0,lambdas,normI,coords);

figure(5)
imagesc(xvals/lambdaOverD,yvals/lambdaOverD,iPSFv_BB);
axis image; 
axis([-2 2 -2 2]);
title('broadband PSF w/ vortex');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Generate coupling maps for each wavelength

% Parameters for Thorlabs SM600
%Using SM2000 in this version of the code
    % link: https://www.thorlabs.com/NewGroupPage9_PF.cfm?ObjectGroup_ID=949
fiber_props.core_rad = 5.5e-6;% Core radius [um]
fiber_props.n_core = 1.4436; %1.4571; core index (interpolated from linear fit to 3 points)
fiber_props.n_clad = 1.4381; %1.4558; cladding index (interpolated from linear fit to 3 points)
fiber_props.type = 'bessel';
Fnum = getMFD(fiber_props,lambda0)/(lambda0*1.4); % focal ratio of the beam at the fiber

eta_maps = generateCouplingMap_polychromatic( Epup.*EPM, fiber_props, lambda0, Fnum, lambdas, totalPower0, lambdaOverD, 3*lambdaOverD, coords);

figure(6);
for ch = 1:numWavelengths
    subplot(1,numWavelengths,ch);
    imagesc(xvals/lambdaOverD,yvals/lambdaOverD,log10(eta_maps(:,:,ch)));
    axis image; 
    axis([-2 2 -2 2]);
    caxis([-3 -0.5])
    title(['\eta at ',num2str(lambdas(ch)*1e9),'nm']);
    colorbar;
    colormap(gray(256));
end

% find the centroid of eta_maps(:,:,ch)
% find min in centroid
% find coordinates of this min in eta_maps
% plot the difference in these coordinates in x and y position from the
% origin, with lambda corresponding to ch on the bottom and offset on the
% left
% 
% Finds the radius of the centroid to crop the image to depending on the maximum eta in 1 slice
% of eta_maps(:,:,ch)
Xshift = zeros(numWavelengths,1);
Yshift = zeros(numWavelengths,1);
etas = zeros(numWavelengths, 1);

for ch = 1:numWavelengths 
    map = eta_maps(:,:,ch); %one slice of the eta_maps cube
    map_max = max(map,[],'all'); %the maximum value in cmap
    [max_ind(1),max_ind(2)] = find(map_max==map,1); %linear coordinates of max value
    max_rho = sqrt(((N/2+1)-max_ind(1))^2 + ((N/2+1)-max_ind(2))^2);
    
    crp = 2*max_rho; %The length of one side of the cube to crop the image to

    cmap = map(end/2+1-floor(crp/2):end/2+1+floor(crp/2),end/2+1-floor(crp/2):end/2+1+floor(crp/2)); %the centroid
    cmap_min = min(cmap,[],'all'); %minimum value in the centroid
    [min_ind(1),min_ind(2)] = find(cmap_min==cmap); %indices of minimum value in the centroid
    min_ind = min_ind + (N/2-floor(crp/2)); %adjust min values to reflect position in map, not cmap
   
    Xshift(ch) = N/2-min_ind(2); %x value is wavelength, y value is offset
    Yshift(ch) = N/2-min_ind(1);%^
    
    etas(ch) = cmap_min;
end

%Null shift plots for X, Y, and the trendlines that result
figure(7);
subplot(2,2,1);
plot(lambdas/lambda0,Xshift/lambdaOverD, '-o', 'Color', 'r');
title('Xshift');
xlabel('\lambda/\lambda0')
ylabel('\eta')
subplot(2,2,2);
plot(lambdas/lambda0,Yshift/lambdaOverD, '-o', 'Color', 'b');
title('Yshift');
xlabel('\lambda/\lambda0')
ylabel('\eta')

px = polyfit(lambdas/lambda0,Xshift'/lambdaOverD,1);
pxy = polyval(px,lambdas/lambda0);
subplot(2,2,3);
plot(lambdas/lambda0,pxy,'-o','Color','m')
title('X Offset Trend')
xlabel('\lambda/\lambda0')
ylabel('\eta')
txt = ['p value: ' num2str(px)];
text(mean(lambdas/lambda0),mean(pxy),txt);

py = polyfit(lambdas/lambda0,Yshift'/lambdaOverD,1);
pyy = polyval(py,lambdas/lambda0);
subplot(2,2,4);
plot(lambdas/lambda0,pyy,'-o','Color','g');
title('Y Offset Trend')
xlabel('\lambda/\lambda0')
ylabel('\eta')
txt = ['p value: ' num2str(py)];
text(mean(lambdas/lambda0),mean(px),txt);

%Null value vs wavelength offset from central wavelength
figure(8);
subplot(1,1,1);
semilogy(lambdas/lambda0,etas,'-o','Color','r'); %lambdas/lambda0,,'-o','Color','r'
title('Null Value vs \lambda/\lambda0')
xlabel('\lambda/\lambda0')
ylabel('\eta')
grid on

%Actual positional offset of null for x and y overlayed
figure(9);
plot(lambdas/lambda0,Xshift/lambdaOverD, '-o', 'Color', 'r');
hold on
plot(lambdas/lambda0,Yshift/lambdaOverD, '-o', 'Color', 'b');
legend({'Xshift', 'Yshift'}, 'Location', 'SouthEast');
title('Null Movement')
xlabel('\lambda/\lambda_0')
ylabel('Null Shift [\lambda/D]')
grid on

%Overlay of trends in x and y null positional offset
figure(12);
plot(lambdas/lambda0,pxy,'-o','Color','r');
hold on
plot(lambdas/lambda0,pyy,'-o','Color','b');
legend({'Xshift', 'Yshift'}, 'Location', 'SouthEast');
title('Null Movement')
xlabel('\lambda/\lambda_0')
ylabel('Null Shift [\lambda/D]')
txt = ['y trend p value: ' num2str(py) newline 'x trend p value: ' num2str(px)];
text(mean(lambdas/lambda0),mean(px),txt);
grid on





