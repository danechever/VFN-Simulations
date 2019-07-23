%{
Test code used for developing the basic format of the SPIE 2019 KPIC VFN
simulations. This is based off the main_pupilBasedVFN_polychromatic code.

Goals:___
- Develop basic formatting and method for main SPIE figure
- Simulate VFN on-sky performance given known WF residuals
- Show a live image of this performance vs. time.

Notes:___
* Beg. of code is mostly untouched. Changes start in "Generate SPIE fig" section
* Wavefront residuals are plotted in [nm] at central lambda
* Eta_s is calculated using min in central region:
    * Multi-wavelength eta matrix gets cropped to central 1/32nd of frame
    * All wavelengths are averaged to get broadband null 
    * Min of this wavelength-averaged sub-matrix is treated as null
    --> Thus, the calculated null puts the fiber at the best location for a
        **broadband** null assuming a flat spectrum. 
    --> This has implications to be considered when the null splits.
    * The wavelength eta_s's (for the final plot) are calculated with the fiber
        at the location that produced the best broadband null, NOT the ideal
        position for each individual wavelength.
* Y-axis of eta_s plots are rescaled by a user-provided value to avoid having
    the extra exponent at the top of the figure.
* Eta_p is calculated via azimuthal averaging
    * VFN_An_radAverage (from VFN-Lab repo) is used for this averaging
    * The centering of the average is based on the **broadband** null location
        calculated as described above. 
    * The eta_p reported in the bottom left graph is the average accross the
        full band based on the location that produced the best **average** (ie.
        broadband) planet coupling. This is the planet location with the fiber
        at the ideal null position.
    * The wavelength-dependent couplings (for the final plot) are calculated
        with the fiber located at the ideal null position and planet located at
        the ideal **average/broadband** location as described in last bullet but
        at each wavelength.
* A circle is superimposed on the eta map to show where the ideal planet
    location is based on the Eta_p calculation described above. The average at
    this circle's location is the value reported as eta_p everywhere else.
%}

clear; %close all; 
addpath('VFNlib');
addpath(genpath('C:\Users\danie\Documents\MATLAB\VFN-Lab\AnalysisCode\'))

%% Input parameters 

% Define smapling info
N = 2^11; % Size of computational grid (NxN samples) 
apRad = 128; % Aperture radius in samples 

% Define wavelength info
lambda0 = 780e-9; %central wavelength
fracBW = 0.18; %\Delta\lambda/\lambda
numWavelengths = 5;% number of discrete wavelengths 
lambdas = getWavelengthVec(lambda0,fracBW,numWavelengths);% array of wavelengths (meters)

% Define charge of the vortex mask at each wavelength
charge = 2*ones(1,numWavelengths); % achromatic
%charge = 2*lambda0./lambdas; % simple scalar model

% Define wavefront error at the central wavelength
nolls = 4:8;
coeffs = 0.5*[0, 0.025, 0.0, 0.025, 0.0];

% Give offsets for the vortex mask
offsetX = 0*apRad;%0.0952*apRad;
offsetY = 0*apRad;%0.0524*apRad; 

% Flag to plot coupling maps in log scale
islogcoup = true;

%% Generate the coordinate system

coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
xvals = coords.xvals;% Helpful for plotting
yvals = coords.yvals;

%% Create array with pupil function

PUPIL = makeCircularPupil( apRad, N );
[normI, totalPower0] = getNormalization(PUPIL);% Normalization factors
lambdaOverD = N/apRad/2; % lam/D in units of samples in the image plane

% figure(1)
% imagesc(xvals/apRad,yvals/apRad,PUPIL);
% axis image; 
% axis([-1 1 -1 1]);
% title('Pupil');
% colorbar; 
% colormap(parula(256));
% drawnow;

%% Define pupil field

phz = generateZernike_fromList( nolls, coeffs, PUPIL, apRad, coords); 

figure(2);
for ch = 1:numWavelengths
    Epup(:,:,ch) = exp(1i*phz*lambda0/lambdas(ch)).*PUPIL;
    
    subplot(1,numWavelengths,ch);
    imagesc(xvals/apRad,yvals/apRad,angle(Epup(:,:,ch)));
    axis image; 
    axis([-1 1 -1 1]);
    title(['Phase at ',num2str(lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(parula(256));
   
end
drawnow;

%% Get PSF without vortex mask

% Get broadband PSF
iPSF_BB = getPSF(Epup,lambda0,lambdas,normI,coords);

figure(3)
imagesc(xvals/lambdaOverD,yvals/lambdaOverD,iPSF_BB);
axis image; 
axis([-3 3 -3 3]);
title('broadband PSF w/o vortex');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Make vortex mask 

EPM = generateVortexMask( charge, coords, [offsetX offsetY] );

central_band_index = ceil(numWavelengths/2);

figure(4)
imagesc(xvals/apRad,yvals/apRad,angle(EPM(:,:,central_band_index).*Epup(:,:,central_band_index)));
axis image; 
axis([-1 1 -1 1]);
title('Pupil phase at \lambda_0');
colormap(hsv(256));
colorbar; 
drawnow;

%% Get PSF with vortex mask

iPSFv_BB = getPSF(Epup.*EPM,lambda0,lambdas,normI,coords);

figure(5)
imagesc(xvals/lambdaOverD,yvals/lambdaOverD,iPSFv_BB);
axis image; 
axis([-3 3 -3 3]);
title('broadband PSF w/ vortex');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Generate coupling maps for each wavelength

% Parameters for Thorlabs SM600
    % link: https://www.thorlabs.com/NewGroupPage9_PF.cfm?ObjectGroup_ID=949
fiber_props.core_rad = 11e-6/2;% Core radius [um]
fiber_props.n_core = 1.4606;% core index (interpolated from linear fit to 3 points)
fiber_props.n_clad = 1.4571;% cladding index (interpolated from linear fit to 3 points)
Fnum = 5; % focal ratio of the beam at the fiber
fiber_props.type = 'bessel';

eta_maps = generateCouplingMap_polychromatic( Epup.*EPM, fiber_props, lambda0, Fnum, lambdas, totalPower0, lambdaOverD, 3*lambdaOverD, coords);

figure(6);
for ch = 1:numWavelengths
    subplot(1,numWavelengths,ch);
    if islogcoup
        imagesc(xvals/lambdaOverD,yvals/lambdaOverD,log10(eta_maps(:,:,ch)));
    else
        imagesc(xvals/lambdaOverD,yvals/lambdaOverD,eta_maps(:,:,ch));
    end
    axis image;
    axis([-2 2 -2 2]);
    title(['\eta at ',num2str(lambdas(ch)*1e9),'nm']);
    colorbar;
    if islogcoup; caxis([-4 0]); end
    colormap(gray(256));
end

%% Generate SPIE fig

figure(7);

%-- Get central wavelength index
lam0 = ceil(numWavelengths/2);

%-- Axis limits
% PSF plot
axlimPSF = [-3 3 -3 3];     % x-y axis limits
% Eta Map
axlim2DETA = axlimPSF;      % x-y axis limits
% Eta_p plots
axYlimEtaP = [0 12];        % y axis limit
% Eta_s plots
axYlimEtaS = [0 0.25];      % y axis limit
eta_sYSCL  = 1e3;           % scaling for y axis

%-- Plot instantaneous WFE/R (at central lambda)
% Convert WFE from rads to nm
wfr = phz/2/pi*lambda0*1e9;
% Plot
subplot(2,3,1); 
imagesc(xvals/apRad,yvals/apRad,wfr);
axis image; 
axis([-1 1 -1 1]);
title(['Wavefront Residuals at ',num2str(lambdas(lam0)*1e9),'nm']);
xlabel('x/D', 'FontWeight', 'bold')
ylabel('y/D', 'FontWeight', 'bold')
wfrCBAR = colorbar; 
wfrCBAR.Label.String = 'Residuals [nm]';
wfrCBAR.Label.FontWeight = 'bold';
colormap(parula(256));

%-- Plot instantaneous Donut PSF (polychromatic)
subplot(2,3,2);
imagesc(xvals/lambdaOverD,yvals/lambdaOverD,iPSFv_BB);
axis image; 
axis(axlimPSF);
title('Broadband PSF w/ vortex');
xlabel('\lambda/D', 'FontWeight', 'bold')
ylabel('\lambda/D', 'FontWeight', 'bold')
colorbar;%caxis([-3 0])
colormap(parula(256));

%-- Plot eta_s vs. time [eta_s first since it's needed for eta_p and eta map]
% Get eta_s as min in center of all eta_maps averaged together
    % Average eta_maps since fiber can only be in one location so need to make
    % sure we take the eta_s at this same locaiton always.
crp = 68;   %1/2 of crop window for min analysis
eta_cent = eta_maps(round(end/2-end/crp):round(end/2+end/crp),round(end/2-end/crp):round(end/2+end/crp),:);
eta_cent = mean(eta_cent,3);
eta_s = min(eta_cent,[],[1 2]);    % Get min of average in central region as broadband null
% Get indices for null location
[eta_sInd(1), eta_sInd(2)] = find(eta_s==eta_cent);
eta_sInd = eta_sInd + round(N/2-N/crp-1);   % Account for cropping from before
% Rescale broadband null value for plot
eta_s = eta_s*eta_sYSCL;   
% Plot
subplot(2,3,5);
scatter(1, eta_s)
ylim(axYlimEtaS)
ylabel('\eta_s [\times 10^{-3}]', 'FontWeight', 'bold')
xlabel('WF Sample', 'FontWeight', 'bold')
title('Broadband Star Coupling')

%-- Plot eta_p vs. time [eta_p before eta map since needed for circle plot]
% Iterate through wavelengths to get azimuthal average of coupling
%eta_sI = nan(numWavelengths, 2);
eta_pAvg = nan(numWavelengths, N/2);
for i = 1:numWavelengths
    % Azimuthally average
    tmpAvg = VFN_An_radAverage(eta_maps(:,:,i), eta_sInd);
    % Only populate part of matrix that is valid from tmpAvg
    eta_pAvg(i,1:length(tmpAvg)) = tmpAvg;
end
% Average along wavelength dimension to find best average planet coupling
[eta_p, eta_pInd] = max(mean(eta_pAvg, 1));
eta_p = eta_p*100;  % convert to [%]
% Plot
subplot(2,3,4);
scatter(1, eta_p)
ylim(axYlimEtaP)
ylabel('\eta_p [%]', 'FontWeight', 'bold')
xlabel('WF Sample', 'FontWeight', 'bold')
title('Broadband Planet Coupling')

%-- Plot instantaneous eta map (at central lambda)
subplot(2,3,3);
imagesc(xvals/lambdaOverD,yvals/lambdaOverD,eta_maps(:,:,lam0));
% Show circle where coupling was averaged. 
    % fliplr needed to account for non-conventional [col row] indexing of viscircles
    % subtract -2 and max(xvals)... to center on image axes.
    % subtract -1 from eta_pInd since index 1 but axes are index 0
viscircles(fliplr((eta_sInd-2)/lambdaOverD-max(xvals)/lambdaOverD),  (eta_pInd-1)/lambdaOverD);
axis image;
axis(axlim2DETA);
title(['\eta at ',num2str(lambdas(lam0)*1e9),'nm']);
xlabel('\lambda/D', 'FontWeight', 'bold')
ylabel('\lambda/D', 'FontWeight', 'bold')
colorbar;
colormap(gray(256));

%-- Plot eta_p and eta_s vs. wavelength
% Eta_s vs wavelength (at eta_sInd location)
eta_sL = squeeze(eta_maps(eta_sInd(1), eta_sInd(2), :));
% Eta_p vs wavelength (average w/ centering based on eta_sInd)
eta_pL = eta_pAvg(:,eta_pInd)*100;
subplot(2,3,6);
title('Coupling Fractions Across Band')
yyaxis left
plot(lambdas*1e9, eta_sL*eta_sYSCL)
hold on
scatter(lambdas*1e9, eta_sL*eta_sYSCL)
hold off
ylim(axYlimEtaS)
ylabel(['\eta_s [\times 10^{-' num2str(log10(eta_sYSCL)) '}]'], 'FontWeight', 'bold')
xlabel('Wavelength [nm]', 'FontWeight', 'bold')
yyaxis right
plot(lambdas*1e9, eta_pL)
hold on
scatter(lambdas*1e9, eta_pL)
hold off
ylim(axYlimEtaP)
ylabel('\eta_p [%]', 'FontWeight', 'bold')


%% TO-DO LIST!
%{
1) Rescale WF plot to nm
    *** DONE
    - NOTE: phz matrix is same as what is plotted in WFR currently
    - phz/2*pi = waves
    - Generate WF in nm (using my code). Convert to waves. get rms. Generate in
    garys code using rms waves. Get back to nm.
2) get eta_p via averaging instead of max
    *** DONE
2a) Add circle on coupling plot showing where average coupling was done
    *** DONE
3) modify code to loop through multiple frame calculations
    - have each plot in figure be its own matrix
    - save each frame as an itr axis in matrices
4) use rand() to generate random wavefronts and create time-dependent plots
5) create fully dynamic plot (ie. gif)
%}