%This script simulates pupil-plane scalar VFN by generating a vortex phase
%pattern using the segments of the Keck primary. Conducts analysis of the
%resulting coupling maps in a manner similar to that shown for a circular
%pupil in Ruane et. al. 2019. This can be done by either loading a 'seed'
%array of tilt and piston from a matlab file, or using the analytical
%solution that is generated by default.
%Updated 12/2021 to use MFT propagation instead of FFT
%Updated 06/2022 to improve readability and clean up the UX

clc;
clear; 
close all; 

% Set Paths
addpath(genpath(fullfile('..','VFNlib')));
addpath(genpath(fullfile('..','..','falco-matlab')));
addpath(genpath(fullfile('..','..','VFN-Lab')));

% Import Python functions
py.importlib.import_module('sellmeir1');
py.importlib.import_module('triple_prism');
py.importlib.import_module('singleprism');

% Provide path where the DM Basis file is located
inpar.dmBasisPath = '/media/Data_Drive/VFN/ScalarOnPrimaryData/';
%%
sellmeir1coeffs;
%% Load Input Parameters

% Edit this .m file to configure the basic parameters of the simulation
sppvfnConfig;

%% Wedge and ADC Parameters

inpar.magfactor = 890.16;

cType = 'wedge';

%% Generate Coordinate System

sppvfnCoordSysConfig;

inpar.coordsPP = coordsPP;
inpar.coordsFP = coordsFP;

%% Create array with pupil function

%Edit pupilType to use a different pupil
%***CAUTION***
%May need to modify input parameters
pupilType = 'keck';

PUPIL = genPup(pupilType, inpar);

totalPower0 = sum(abs(PUPIL(:)));

inpar.lam0OverD_meters = lambda0Fnum_meters/foc;

%Displays the selected pupil
figure();
imagesc(inpar.xvalsPP,inpar.yvalsPP,PUPIL); 
axis image;
axis([-1 1 -1 1]);
title('Pupil');
colorbar; 
colormap(parula(256));
drawnow;

addpath(['..' filesep '..' filesep 'falco-matlab' filesep 'lib' filesep 'utils']);

 %% Define pupil field with vortex mask
 
addpath(['..' filesep '..' filesep 'falco-matlab' filesep 'lib' filesep 'masks']);
 
vType = 'spp';

[phz, pType] = vortexSelection(vType, inpar);
 
% Phase pattern plot
figure();
imagesc(phz); 
axis image; 
colorbar;
colormap(hsv);
set(gca,'ydir','normal')
title(['Charge 1 ' pType ' Vortex Phase (radians)'])
caxis([-pi pi]);

% Phase pattern w/ pupil
figure();
Epup = nan(inpar.N,inpar.N,inpar.numWavelengths);
for ch = 1:inpar.numWavelengths
    
    Epup(:,:,ch) = exp(1i*phz*inpar.lambda0/inpar.lambdas(ch)).*PUPIL;
    
    subplot(1,inpar.numWavelengths,ch);
    imagesc(inpar.xvalsPP,inpar.yvalsPP,angle(Epup(:,:,ch)));
    axis image; 
    axis xy;
    axis([-1 1 -1 1]);
    title(['Phase at ',num2str(inpar.lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(hsv);
end
drawnow;

figure();
imagesc(inpar.xvalsPP,inpar.yvalsPP,angle(Epup(:,:,ceil(inpar.numWavelengths/2))));
axis image; 
axis xy;
axis([-1 1 -1 1]);
title([pType ' vortex on Keck primary pupil at ',num2str(inpar.lambdas(ch)*1e9),'nm']);
colorbar; 
colormap(hsv);

%% Get PSF with vortex mask
addpath(['..' filesep '..' filesep 'falco-matlab' filesep 'lib' filesep 'propcustom']);

[iPSFv_BB, PSFv] = getPSF_mft(Epup, inpar.lambdas, foc, coordsPP, coordsFP);

figure();
imagesc(inpar.xvalsFP,inpar.yvalsFP,iPSFv_BB);
axis image; 
axis([-3 3 -3 3]);
title('Charge 1 Broadband PSF w/ Vortex');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Generate fibermode at each lambda     
%-- Iterate through wavelengths generating modes
fibmode = nan(Nxi, Nxi, inpar.numWavelengths);
for ch = 1:inpar.numWavelengths
    % Generate the fiber mode for this wavelength with proper scaling
	fibmode(:,:,ch) = generateSMFmode_mft( fiber_props, inpar.lambdas(ch), coordsFP.dx, coordsFP);
    
end

%% Calculate Throughputs
%-- Get null depth (using overlap integral)
eta_onAx = nan(1,inpar.numWavelengths);
for ch = 1:inpar.numWavelengths
    eta_onAx(ch) = (abs(sum(sum(PSFv(:,:,ch).*fibmode(:,:,ch)))).^2)/totalPower0;
end

%-- Compute 2D coupling map
eta_maps = zeros(Nxi,Nxi,inpar.numWavelengths);
for ch = 1:inpar.numWavelengths
    % Compute the monochromatic coupling map (spatial units of lambda/D)
    eta_maps(:,:,ch) = generateCouplingMap( fibmode(:,:,ch), PSFv(:,:,ch), totalPower0, 5*lambda0Fnum_samp, coordsFP);

end

disp('Key Coupling Points:')
for ch = 1:inpar.numWavelengths
    fprintf('lambda = %f nm,    on-axis coup = %e,    max = %f %%\n',inpar.lambdas(ch)*1e9, eta_onAx(ch), max(eta_maps(:,:,ch),[],'all')*100);
end

%% Display coupling maps
%-- Linear scale
figure();
for ch = 1:inpar.numWavelengths
    subplot(1,inpar.numWavelengths,ch);
    imagesc(inpar.xvalsFP,inpar.yvalsFP,eta_maps(:,:,ch));
    axis image; 
    axis([-3 3 -3 3]);
    title(['\eta at ',num2str(inpar.lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(gray(256));
end
drawnow

%-- Log Scale
figure();
for ch = 1:inpar.numWavelengths
    subplot(1,inpar.numWavelengths,ch);
    imagesc(inpar.xvalsFP,inpar.yvalsFP,log10(eta_maps(:,:,ch)));
    axis image; 
    axis([-3 3 -3 3]);
    title(['log10(\eta) at ',num2str(inpar.lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(gray(256));
end
drawnow

%% Display broadband coupling map if applicable
% This is what the photodiode should see when broadband light is applied

if inpar.numWavelengths > 1
    % Compute BB coupling map as average of all coupling maps
      % (This assumes the input spectrum is is flat --> ie. equal power at
      % all wavelengths)
    eta_map_BB = mean(eta_maps,3);
    
    % Plot Linear Scale
    figure();
    imagesc(inpar.xvalsFP, inpar.yvalsFP, eta_map_BB);
    axis image;
    axis([-3 3 -3 3]);
    title('Charge 1 Broadband \eta');
    colorbar;
    colormap(gray(256));
    
    % Plot Log Scale
    figure();
    imagesc(inpar.xvalsFP, inpar.yvalsFP, log10(eta_map_BB));
    axis image;
    axis([-3 3 -3 3]);
    title('Broadband log10(\eta)');
    colorbar;
    colormap(gray(256));
    
    % Print key values
    disp('---')
    fprintf('Broadband Performance:     on-axis coup = %e,    max = %f %%\n', mean(eta_onAx), max(eta_map_BB,[],'all')*100);
end

%% Characterize Null Location
%-- find the centroid of eta_maps(:,:,ch)
%-- find min in centroid
%-- find coordinates of this min in eta_maps
%-- plot the difference in these coordinates in x and y position from the
%   origin, with lambda corresponding to ch on the bottom and offset on the
%   left
% 
%-- Finds the radius of the centroid to crop the image to depending on the maximum eta in 1 slice
%   of eta_maps(:,:,ch)

Xshift = zeros(inpar.numWavelengths,1);
Yshift = zeros(inpar.numWavelengths,1);
etas = zeros(inpar.numWavelengths, 1);

for ch = 1:inpar.numWavelengths 
    map = eta_maps(:,:,ch); %one slice of the eta_maps cube
    map_max = max(map,[],'all'); %the maximum value in cmap
    [max_ind(1),max_ind(2)] = find(map_max==map,1); %linear coordinates of max value
    max_rho = sqrt(((Nxi/2+1)-max_ind(1))^2 + ((Nxi/2+1)-max_ind(2))^2);
    
    crp = 1.5*max_rho; %The length of one side of the cube to crop the image to

    cmap = map(end/2+1-floor(crp/2):end/2+1+floor(crp/2),end/2+1-floor(crp/2):end/2+1+floor(crp/2)); %the centroid
    cmap_min = min(cmap,[],'all'); %minimum value in the centroid
    [min_ind(1),min_ind(2)] = find(cmap_min==cmap); %indices of minimum value in the centroid
    min_ind = min_ind + (Nxi/2-floor(crp/2)); %adjust min values to reflect position in map, not cmap
   
    Xshift(ch) = Nxi/2+1-min_ind(2); %x value is wavelength, y value is offset
    Yshift(ch) = Nxi/2+1-min_ind(1);%^
    
    etas(ch) = cmap_min;
end

%change lambdaOverD here 
inpar.lambdaOverD = lambda0Fnum_samp; %inpar.N/inpar.apRad/2; % lam/D in units of samples in the image plane

%-- Null shift plots for X, Y, and the trendlines that result
figure();
subplot(2,2,1);
plot(inpar.lambdas/inpar.lambda0,Xshift/inpar.lambdaOverD, '-o', 'Color', 'r');
title('Xshift');
subplot(2,2,2);
plot(inpar.lambdas/inpar.lambda0,Yshift/inpar.lambdaOverD, '-o', 'Color', 'b');
title('Yshift');

px = polyfit(inpar.lambdas/inpar.lambda0,Xshift'/inpar.lambdaOverD,1);
pxy = polyval(px,inpar.lambdas/inpar.lambda0);
subplot(2,2,3);
plot(inpar.lambdas/inpar.lambda0,pxy,'-o','Color','m')
title('X Offset Trend')
txt = ['p value: ' num2str(px(1))];
text(mean(inpar.lambdas/inpar.lambda0),mean(pxy),txt);

py_1 = polyfit(inpar.lambdas/inpar.lambda0,Yshift'/inpar.lambdaOverD,1);
pyy = polyval(py_1,inpar.lambdas/inpar.lambda0);
subplot(2,2,4);
plot(inpar.lambdas/inpar.lambda0,pyy,'-o','Color','g');
title('Y Offset Trend')
txt = ['p value: ' num2str(py_1(1))];
text(mean(inpar.lambdas/inpar.lambda0),mean(pyy),txt);

%-- Null value vs wavelength offset from central wavelength
figure();
subplot(1,1,1);
semilogy(inpar.lambdas/inpar.lambda0,eta_onAx','-o','Color','r'); %lambdas/lambda0,,'-o','Color','r'
ylim([1e-8 1]);
title(['Charge 1 ' pType ' VFN Null Value vs \lambda/\lambda0'])
xlabel('\lambda/\lambda0')
ylabel('\eta')
grid on

%-- Actual positional offset of null for x and y overlayed
figure();
plot(inpar.lambdas/inpar.lambda0,Xshift/inpar.lambdaOverD, '-o', 'Color', 'r');
hold on
plot(inpar.lambdas/inpar.lambda0,Yshift/inpar.lambdaOverD, '-o', 'Color', 'b');
legend({'Xshift', 'Yshift'}, 'Location', 'SouthEast');
title('Null Movement')
xlabel('\lambda/\lambda_0')
ylabel('Null Shift [\lambda/D]')
grid on

%-- Overlay of trends in x and y null positional offset
figure();
plot(inpar.lambdas/inpar.lambda0,pxy,'-o','Color','r');
hold on
plot(inpar.lambdas/inpar.lambda0,pyy,'-o','Color','b');
legend({'Xshift', 'Yshift'}, 'Location', 'SouthWest');
title(['Charge 1 ' pType ' Vortex Induced Null Movement'])
xlabel('\lambda/\lambda_0')
ylabel('Null Shift [\lambda/D]')
txt = ['Y Shift Slope p = ' num2str(py_1(1)) newline 'X Shift Slope p = ' num2str(px(1))]; %newline 'x trend p value: ' num2str(px)];
text(1,0.022,txt);
grid on

% Show the Deepest Null
% figure();
% imagesc(inpar.xvalsFP, inpar.yvalsFP, eta_map_BB);
% plot(min_ind_
% axis image;
% axis([-3 3 -3 3]);
% title('Broadband \eta');
% colorbar;
% colormap(gray(256));

%% Compute the relative integration time
etamean = mean(eta_maps,3);

map_max = max(etamean,[],'all'); %the maximum value in cmap
[max_ind(1),max_ind(2)] = find(map_max==etamean,1); %linear coordinates of max value
max_rho = sqrt(((Nxi/2+1)-max_ind(1))^2 + ((Nxi/2+1)-max_ind(2))^2);
    
crp = 1.5*max_rho; %The length of one side of the cube to crop the image to

cmapmean = etamean(end/2+1-floor(crp/2):end/2+1+floor(crp/2),end/2+1-floor(crp/2):end/2+1+floor(crp/2)); %the centroid
eta_s = min(cmapmean,[],'all'); %minimum value in the centroid
[min_ind(1),min_ind(2)] = find(eta_s==cmapmean); %indices of minimum value in the centroid
min_ind = min_ind + (Nxi/2-floor(crp/2)); %adjust min values to reflect position in map, not cmap
        
eta_sX = min_ind(2);
eta_sY = min_ind(1);

[radProf2, ~] = VFN_An_radAverage(map,[eta_sX, eta_sY]);
[radMx, ~] = max(radProf2);

%Compute the relative integration time here
reltints = eta_s./(radMx.^2);
fprintf(['Relative Integration Time: ' num2str(reltints) '\n']);

%% Duplicate and Save Variables
etasCopy = etas;
pyyCopy = pyy;
YshiftCopy = Yshift;
eta_onAxCopy = eta_onAx;
py_1Copy = py_1;
save('oldetas.mat','etasCopy');
save('pyyo.mat','pyyCopy');
save('oldYshift.mat','YshiftCopy');
save('oldetasonax.mat','eta_onAxCopy');
save('py1o.mat','py_1Copy');
inpar.etasCopy = etas;
inpar.pyyCopy = pyy;
inpar.YshiftCopy = Yshift;
inpar.eta_onAxCopy = eta_onAx;
inpar.py_1Copy = py_1;

%% Correct Wavelength Dependent Null Shift

%ADC Wedge Data
wedge_ADC_1 = 'BaF2';
wedge_ADC_2 = 'CaF2';
wedge_ADC_3 = 'ZnSe';

% inpar.wedge_ADC = 0.5 * pi/180; %Thorlabs 30 arcminute CaF2 wedge
% inpar.n1_ADC = getRefractiveIndex(wedge_ADC_1, 1e6*inpar.lambda0);
% inpar.n2_ADC = getRefractiveIndex(wedge_ADC_2, 1e6*inpar.lambda0);
% inpar.n3_ADC = getRefractiveIndex(wedge_ADC_3, 1e6*inpar.lambda0);

inpar.phi1_ADC = 7.0516 * pi/180;
inpar.phi2_ADC = 3.8050 * pi/180;
inpar.phi3_ADC = 1.1465 * pi/180;

wvs = py.numpy.array(1e6*inpar.lambdas);
inpar.n0 = 1.000293;
inpar.n1 = py.sellmeir1.sellmeir1(wvs, 273.15 + 3, baf2_args(1,:),baf2_args(2,:),baf2_args(3,:),baf2_args(4,:),baf2_args(5,:),baf2_args(6,:));
inpar.n2 = py.sellmeir1.sellmeir1(wvs, 273.15 + 3, caf2_args(1,:),caf2_args(2,:),caf2_args(3,:),caf2_args(4,:),caf2_args(5,:),caf2_args(6,:));
inpar.n3 = py.sellmeir1.sellmeir1(wvs, 273.15 + 3, znse_args(1,:),znse_args(2,:),znse_args(3,:),znse_args(4,:),znse_args(5,:),znse_args(6,:));


figure();
if strcmpi(cType, 'adc')
    
    [inpar.tilt_valsCopy, inpar.clocking, inpar.I] = simple_ADCOPT(inpar);
    
    wphz = simpleADC(inpar);
    
    for ch = 1:inpar.numWavelengths
                 
        phzW = phz - wphz(:,:,ch);%phz + wphz(ch);
    
        %Epup(:,:,ch) = exp(1i*phz*inpar.lambda0/inpar.lambdas(ch)).*inpar.PUPIL;
        Epupw(:,:,ch) = exp(1i*wphz(:,:,ch)).*PUPIL;
        Epup_Wedge(:,:,ch) = exp(1i*phzW).*PUPIL;
        
        disp('here')
        disp(max(max(wphz(:,:,ch)/2/pi))-min(min(wphz(:,:,ch)/2/pi)));

    
        subplot(1,inpar.numWavelengths,ch);
        imagesc(inpar.xvalsPP,inpar.yvalsPP,angle(Epupw(:,:,ch)));
        axis image; 
        axis xy;
        axis([-1 1 -1 1]);
        title(['Phase at ',num2str(inpar.lambdas(ch)*1e9),'nm']);
        colorbar; 
        colormap(hsv);
    end
    
elseif strcmpi(cType, 'wedge')
    
    [inpar.tilt_valsCopy, inpar.clocking, inpar.I] = simple_PRISMOPT_WIP(inpar);
    
    wphz = simpleADC(inpar);
    
    for ch = 1:inpar.numWavelengths
        
        phzW = phz - wphz(:,:,ch);%phz + wphz(ch);
    
        %Epup(:,:,ch) = exp(1i*phz*inpar.lambda0/inpar.lambdas(ch)).*inpar.PUPIL;
        Epupw(:,:,ch) = exp(1i*wphz(:,:,ch)).*PUPIL;
        Epup_Wedge(:,:,ch) = exp(1i*phzW).*PUPIL;
%         inpar.pyyo = -1*pyyCopy(ch);
%     
%         %Requisite difference in phase computed here
%         wphz = simplewedge_max(inpar, inpar.lambdas(ch));
%     
%         phzW = (phz*inpar.lambda0/inpar.lambdas(ch)) - wphz;
%     
%         %Epup(:,:,ch) = exp(1i*phz*inpar.lambda0/inpar.lambdas(ch)).*inpar.PUPIL;
%         Epupw(:,:,ch) = exp(1i*wphz).*PUPIL;
%         Epup_Wedge(:,:,ch) = exp(1i*phzW).*PUPIL;
    
        subplot(1,inpar.numWavelengths,ch);
        imagesc(inpar.xvalsPP,inpar.yvalsPP,angle(Epupw(:,:,ch))/2/pi);
        axis image; 
        axis xy;
        axis([-1 1 -1 1]);
        title(['Phase at ',num2str(inpar.lambdas(ch)*1e9),'nm']);
        colorbar; 
        colormap(hsv);
    end
end

%% Get PSF after wedge/ADC applied
[PSFwBB, PSFw] = getPSF_mft(Epup_Wedge, inpar.lambdas, foc, coordsPP, coordsFP);

% figure()
% for ch = 1:inpar.numWavelengths
%     
%     subplot(1,inpar.numWavelengths,ch);
%     imagesc(inpar.xvalsPP/inpar.lambdaOverD,inpar.yvalsPP/inpar.lambdaOverD,PSFw(:,:,ch));
%     axis image;
%     axis xy;
%     axis([-2 2 -2 2]);
%     title(['PSF w/ Vortex at ' num2str(inpar.lambdas(ch)*1e9) 'nm']);
%     colorbar;%caxis([-3 0])
%     colormap(parula(256));
% end
% drawnow;

figure()
imagesc(inpar.xvalsFP,inpar.yvalsFP,PSFwBB);
axis image; 
axis([-3 3 -3 3]);
title('Broadband PSF w/ Vortex');
colorbar;
colormap(parula(256));
drawnow;

%% Calculate CORRECTED Throughputs
%-- Get null depth (using overlap integral)
eta_onAx = nan(1,inpar.numWavelengths);
for ch = 1:inpar.numWavelengths
    eta_onAx(ch) = (abs(sum(sum(PSFw(:,:,ch).*fibmode(:,:,ch)))).^2)/totalPower0;
end

%-- Compute 2D coupling map
eta_maps = zeros(Nxi,Nxi,inpar.numWavelengths);
for ch = 1:inpar.numWavelengths
    % Compute the monochromatic coupling map (spatial units of lambda/D)
    eta_maps(:,:,ch) = generateCouplingMap( fibmode(:,:,ch), PSFw(:,:,ch), totalPower0, 5*lambda0Fnum_samp, coordsFP);

end

disp('Key Coupling Points:')
for ch = 1:inpar.numWavelengths
    fprintf('lambda = %f nm,    on-axis coup = %e,    max = %f %%\n',inpar.lambdas(ch)*1e9, eta_onAx(ch), max(eta_maps(:,:,ch),[],'all')*100);
end

%% Display CORRECTED coupling maps
%-- Linear scale
figure();
for ch = 1:inpar.numWavelengths
    subplot(1,inpar.numWavelengths,ch);
    imagesc(inpar.xvalsFP,inpar.yvalsFP,eta_maps(:,:,ch));
    axis image; 
    axis([-3 3 -3 3]);
    title(['\eta at ',num2str(inpar.lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(gray(256));
end
drawnow

%-- Log Scale
figure();
for ch = 1:inpar.numWavelengths
    subplot(1,inpar.numWavelengths,ch);
    imagesc(inpar.xvalsFP,inpar.yvalsFP,log10(eta_maps(:,:,ch)));
    axis image; 
    axis([-3 3 -3 3]);
    title(['log10(\eta) at ',num2str(inpar.lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(gray(256));
end
drawnow

%% Display CORRECTED broadband coupling map if applicable
% This is what the photodiode should see when broadband light is applied

if inpar.numWavelengths > 1
    % Compute BB coupling map as average of all coupling maps
      % (This assumes the input spectrum is is flat --> ie. equal power at
      % all wavelengths)
    eta_map_BB = mean(eta_maps,3);
    
    % Plot Linear Scale
    figure();
    imagesc(inpar.xvalsFP, inpar.yvalsFP, eta_map_BB);
    axis image;
    axis([-3 3 -3 3]);
    title('Broadband \eta');
    colorbar;
    colormap(gray(256));
    
    % Plot Log Scale
    figure();
    imagesc(inpar.xvalsFP, inpar.yvalsFP, log10(eta_map_BB));
    axis image;
    axis([-3 3 -3 3]);
    title('Broadband log10(\eta)');
    colorbar;
    colormap(gray(256));
    
    % Print key values
    disp('---')
    fprintf('Broadband Performance:     on-axis coup = %e,    max = %f %%\n', mean(eta_onAx), max(eta_map_BB,[],'all')*100);
end

%% Characterize the CORRECTED Null Location
%-- find the centroid of eta_maps(:,:,ch)
%-- find min in centroid
%-- find coordinates of this min in eta_maps
%-- plot the difference in these coordinates in x and y position from the
%   origin, with lambda corresponding to ch on the bottom and offset on the
%   left
% 
%-- Finds the radius of the centroid to crop the image to depending on the maximum eta in 1 slice
%   of eta_maps(:,:,ch)
Xshift = zeros(inpar.numWavelengths,1);
Yshift = zeros(inpar.numWavelengths,1);
etaw = zeros(inpar.numWavelengths, 1);

for ch = 1:inpar.numWavelengths 
    map = eta_maps(:,:,ch); %one slice of the eta_maps cube
    map_max = max(map,[],'all'); %the maximum value in cmap
    [max_ind(1),max_ind(2)] = find(map_max==map,1); %linear coordinates of max value
    max_rho = sqrt(((Nxi/2+1)-max_ind(1))^2 + ((Nxi/2+1)-max_ind(2))^2);
    
    crp = 1.5*max_rho; %The length of one side of the cube to crop the image to

    cmap = map(end/2+1-floor(crp/2):end/2+1+floor(crp/2),end/2+1-floor(crp/2):end/2+1+floor(crp/2)); %the centroid
    cmap_min = min(cmap,[],'all'); %minimum value in the centroid
    [min_ind(1),min_ind(2)] = find(cmap_min==cmap); %indices of minimum value in the centroid
    min_ind = min_ind + (Nxi/2-floor(crp/2)); %adjust min values to reflect position in map, not cmap
   
    Xshift(ch) = Nxi/2+1-min_ind(2); %x value is wavelength, y value is offset
    Yshift(ch) = Nxi/2+1-min_ind(1);%^
    
    etaw(ch) = cmap_min;
end

%change lambdaOverD here 
inpar.lambdaOverD = lambda0Fnum_samp; %inpar.N/inpar.apRad/2; % lam/D in units of samples in the image plane

%-- Null shift plots for X, Y, and the trendlines that result
figure();
subplot(2,2,1);
plot(inpar.lambdas/inpar.lambda0,Xshift/inpar.lambdaOverD, '-o', 'Color', 'r');
title('Xshift');
subplot(2,2,2);
plot(inpar.lambdas/inpar.lambda0,Yshift/inpar.lambdaOverD, '-o', 'Color', 'b');
title('Yshift');

px = polyfit(inpar.lambdas/inpar.lambda0,Xshift'/inpar.lambdaOverD,1);
pxy = polyval(px,inpar.lambdas/inpar.lambda0);
subplot(2,2,3);
plot(inpar.lambdas/inpar.lambda0,pxy,'-o','Color','m')
title('X Offset Trend')
txt = ['p value: ' num2str(px)];
text(mean(inpar.lambdas/inpar.lambda0),mean(pxy),txt);

py_1 = polyfit(inpar.lambdas/inpar.lambda0,Yshift'/inpar.lambdaOverD,1);
pyy = polyval(py_1,inpar.lambdas/inpar.lambda0);
subplot(2,2,4);
plot(inpar.lambdas/inpar.lambda0,pyy,'-o','Color','g');
title('Y Offset Trend')
txt = ['p value: ' num2str(py_1)];
text(mean(inpar.lambdas/inpar.lambda0),mean(px),txt);

%-- Null value vs wavelength offset from central wavelength

figure();
%subplot(1,1,1);
semilogy(inpar.lambdas/inpar.lambda0,eta_onAx','-o','Color','b'); %lambdas/lambda0,,'-o','Color','r'
hold on;
% semilogy(inpar.lambdas/inpar.lambda0,etasCopy,'-o','Color','b');
semilogy(inpar.lambdas/inpar.lambda0,eta_onAxCopy','-o','Color','r');
ylim([1e-15 1]);
legend({'Wedge', 'No Wedge'}, 'Location', 'SouthEast');
title([pType ' with Wedge Null Value vs \lambda/\lambda0'])
xlabel('\lambda/\lambda0')
ylabel('\eta')
grid on
hold off;

%-- Actual positional offset of null for x and y overlayed
figure();
plot(inpar.lambdas/inpar.lambda0,Xshift/inpar.lambdaOverD, '-o', 'Color', 'r');
hold on
plot(inpar.lambdas/inpar.lambda0,Yshift/inpar.lambdaOverD, '-x', 'Color', 'b');
legend({'Xshift', 'Yshift'}, 'Location', 'SouthEast');
title('Null Movement')
xlabel('\lambda/\lambda_0')
ylabel('Null Shift [\lambda/D]')
grid on

%-- Overlay of trends in x and y null positional offset
figure();
plot(inpar.lambdas/inpar.lambda0,pxy,'-x','Color','r');
hold on
plot(inpar.lambdas/inpar.lambda0,pyy,'-o','Color','g');
plot(inpar.lambdas/inpar.lambda0,inpar.pyyCopy,'-o','Color','b');
legend({'Xshift', 'Yshift with ADC','Yshift without ADC'}, 'Location', 'SouthEast');
title([pType ' Induced VFN Null Movement'])
xlabel('\lambda/\lambda_0')
ylabel('Null Shift [\lambda/D]')
txt = ['Shift Slope no Wedge p = ' num2str(py_1Copy(1)) newline 'Shift Slope Wedge p = ' num2str(py_1(1))]; %['y trend p value: ' num2str(py_1) newline 'x trend p value: ' num2str(px)];
text(0.9,0.05,txt);
grid on

%% Compute the relative integration time
etamean = mean(eta_maps,3);

map_max = max(etamean,[],'all'); %the maximum value in cmap
[max_ind(1),max_ind(2)] = find(map_max==etamean,1); %linear coordinates of max value
max_rho = sqrt(((Nxi/2+1)-max_ind(1))^2 + ((Nxi/2+1)-max_ind(2))^2);
    
crp = 1.5*max_rho; %The length of one side of the cube to crop the image to

cmapmean = etamean(end/2+1-floor(crp/2):end/2+1+floor(crp/2),end/2+1-floor(crp/2):end/2+1+floor(crp/2)); %the centroid
eta_s = min(cmapmean,[],'all'); %minimum value in the centroid
[min_ind(1),min_ind(2)] = find(eta_s==cmapmean); %indices of minimum value in the centroid
min_ind = min_ind + (Nxi/2-floor(crp/2)); %adjust min values to reflect position in map, not cmap
        
eta_sX = min_ind(2);
eta_sY = min_ind(1);

[radProf2, ~] = VFN_An_radAverage(map,[eta_sX, eta_sY]);
[radMx, ~] = max(radProf2);

%Compute the relative integration time here
reltints = eta_s./(radMx.^2);
fprintf(['Relative Integration Time: ' num2str(reltints) '\n']);

