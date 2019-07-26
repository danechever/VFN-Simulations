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
    * Multi-wavelength eta matrix gets cropped to central subsection of frame
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
    * [Updated]: User can also select to have axes dynamically scaled by code to
    go from 0 to 1.2*max(eta_sL).
* Eta_p is calculated via azimuthal averaging
    * VFN_An_radAverage (from VFN-Lab repo) is used for this averaging
    * The centering of the average is based on the **broadband** null location
        calculated as described above. 
        ** [Updated]: Centering is based on center of frame since this is where
        the star is actually located and hence the planet distance is w.r.t here
    * The eta_p reported in the bottom left graph is the average accross the
        full band based on the location that produced the best **average** (ie.
        broadband) planet coupling in the first frame.
    * The wavelength-dependent couplings (for the final plot) are calculated
        w.r.t the real star location as described in last bullet but at each 
        wavelength.
    * [Updated]: Averaging is now not just at a single separation, it is over a
        user-defined range of radii. This allows for uncertainty in the planet
        location.
* A circle is superimposed on the eta map to show where the ideal planet
    location is based on the Eta_p calculation described above. The average at
    this circle's location is the value reported as eta_p everywhere else.
    ** [Updated]: Two circles are now used to show the inner and outer range of
    the averaged region.
* Keck pupil is used instead of simple circle.
* Wavefront residuals are modulated accross samples using sine function to have
    continuously variable WF.
* Magnitude of each WF residual across samples is plotted at the end.
* "subtightplot" library is used to provide control of subplot location.
* Full figure is only plotted once and then only the data is updated on each
    iteration of wavefront samples.
* Gif is optionally output at the end.

** NOTE: this version is very RAM intensive since all matrices are preallocated.
    - This provides full matrices at end of run which is useful for debugging
        since you only need to run the actual calculations once and then you can
        modify and run the figure generation separately.
    - Since so much RAM is used (max allowable by MATLAB), the overall runtime
        actually increases since we're memory-limited.
    * Faster, more efficient way to do this is to not save full matrices, just
        the individual frames each time and save the frame at each itr.
%}

clear; %close all; 
addpath('VFNlib');
addpath(genpath('C:\Users\danie\Documents\MATLAB\VFN-Lab\AnalysisCode\'))

%% Input parameters 

% Define smapling info
N = 2^11; % Size of computational grid (NxN samples) 
apRad = 128; % Aperture radius in samples 

% Define wavelength info
lambda0 = 2200e-9; %central wavelength
fracBW = 0.20; %\Delta\lambda/\lambda
numWavelengths = 5;% number of discrete wavelengths 
lambdas = getWavelengthVec(lambda0,fracBW,numWavelengths);% array of wavelengths (meters)

% Define charge of the vortex mask at each wavelength
charge = 2*ones(1,numWavelengths); % achromatic
%charge = 2*lambda0./lambdas; % simple scalar model

% Define wavefront error at the central wavelength
nolls = [5, 6];
coeffs = 1*[0.025, 0.01];
% Create matrix-version of coeffs to simulate changing WFE
wfSamps = 50;                      % Number of WF samples to synthesize
wfModPhzEnd = 4*pi;
wfMod = nan(wfSamps, length(nolls));
for i = 1:length(nolls)
    phzShift = rand*2*pi;
    wfMod(:, i) = sin(linspace(0+phzShift, wfModPhzEnd+phzShift, wfSamps))';    % smoothly modulate WFE amplitude
end
%wfMod = repmat(wfMod, 1, length(coeffs));
coeffs = repmat(coeffs, wfSamps, 1);
coeffs = coeffs.*wfMod;             % Apply modulation

% Give offsets for the vortex mask
offsetX = 0*apRad;%0.0952*apRad;
offsetY = 0*apRad;%0.0524*apRad; 

% Flag to plot coupling maps in log scale
islogcoup = true;

% Flag to plot intermediate figure (before SPIE figure)
isPlotSimp = true;

%-- Saving parameters
% Flag to save data
isSaveGif = false;
% Save folder
svfld = 'C:\Users\danie\Desktop\';
svnm  = 'PhaseShiftedGif.gif';

%% Generate the coordinate system

coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
xvals = coords.xvals;% Helpful for plotting
yvals = coords.yvals;

%% Create array with pupil function

%PUPIL = makeCircularPupil( apRad, N );
[PUPIL, pupCircDiam] = makeKeckPupil( 2*apRad, N );
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

%% Make vortex mask (once before loop since does not change)
EPM = generateVortexMask( charge, coords, [offsetX offsetY] );

%% Loop through WF samples

%-- Preallocate data matrices
phzs = nan(wfSamps, N, N);                      % WF residuals
Epups = nan(wfSamps, N, N, numWavelengths);     % Pupil fields
iPSF_BBs = nan(wfSamps, N, N);                  % Broadband PSFs (no vort)
iPSFv_BBs = nan(wfSamps, N, N);                 % Broadband PSFs (vort)
eta_mapss = nan(wfSamps, N, N, numWavelengths); % Eta_maps

for i = 1:wfSamps
fprintf('wfSamp %03d of %03d\n', i, wfSamps);
    
%% Define pupil field
phzs(i,:,:) = generateZernike_fromList( nolls, coeffs(i,:), PUPIL, pupCircDiam/2, coords); 

if isPlotSimp
    figure(2);
end
for ch = 1:numWavelengths
    Epups(i,:,:,ch) = exp(1i*squeeze(phzs(i,:,:))*lambda0/lambdas(ch)).*PUPIL;
    
    if isPlotSimp
        subplot(1,numWavelengths,ch);
        imagesc(xvals/apRad,yvals/apRad,angle(squeeze(Epups(i,:,:,ch))));
        axis image; 
        axis([-1 1 -1 1]);
        title(['Phase at ',num2str(lambdas(ch)*1e9),'nm']);
        colorbar; 
        colormap(parula(256));
    end
   
end
drawnow;

%% Get PSF without vortex mask

% Get broadband PSF
iPSF_BBs(i,:,:) = getPSF(squeeze(Epups(i,:,:,:)),lambda0,lambdas,normI,coords);

if isPlotSimp
    figure(3)
    imagesc(xvals/lambdaOverD,yvals/lambdaOverD,squeeze(iPSF_BBs(i,:,:)));
    axis image; 
    axis([-3 3 -3 3]);
    title('broadband PSF w/o vortex');
    colorbar;%caxis([-3 0])
    colormap(parula(256));
    drawnow;
end 

%% Plot vortex mask
central_band_index = ceil(numWavelengths/2);

if isPlotSimp
    figure(4)
    imagesc(xvals/apRad,yvals/apRad,angle(EPM(:,:,central_band_index).*squeeze(Epups(i,:,:,central_band_index))));
    axis image; 
    axis([-1 1 -1 1]);
    title('Pupil phase at \lambda_0');
    colormap(hsv(256));
    colorbar; 
    drawnow;
end

%% Get PSF with vortex mask

iPSFv_BBs(i,:,:) = getPSF(squeeze(Epups(i,:,:,:)).*EPM,lambda0,lambdas,normI,coords);

if isPlotSimp
    figure(5)
    imagesc(xvals/lambdaOverD,yvals/lambdaOverD,squeeze(iPSFv_BBs(i,:,:)));
    axis image; 
    axis([-3 3 -3 3]);
    title('broadband PSF w/ vortex');
    colorbar;%caxis([-3 0])
    colormap(parula(256));
    drawnow;
end

%% Generate coupling maps for each wavelength

% Parameters for Thorlabs SM600
    % link: https://www.thorlabs.com/NewGroupPage9_PF.cfm?ObjectGroup_ID=949
fiber_props.core_rad = 11e-6/2;% Core radius [um]
fiber_props.n_core = 1.4606;% core index (interpolated from linear fit to 3 points)
fiber_props.n_clad = 1.4571;% cladding index (interpolated from linear fit to 3 points)
Fnum = 5; % focal ratio of the beam at the fiber
fiber_props.type = 'bessel';

eta_mapss(i,:,:,:) = generateCouplingMap_polychromatic( squeeze(Epups(i,:,:,:)).*EPM, fiber_props, lambda0, Fnum, lambdas, totalPower0, lambdaOverD, 3*lambdaOverD, coords);

if isPlotSimp
    figure(6);
    for ch = 1:numWavelengths
        subplot(1,numWavelengths,ch);
        if islogcoup
            imagesc(xvals/lambdaOverD,yvals/lambdaOverD,log10(squeeze(eta_mapss(i,:,:,ch))));
        else
            imagesc(xvals/lambdaOverD,yvals/lambdaOverD,squeeze(eta_mapss(i,:,:,ch)));
        end
        axis image;
        axis([-2 2 -2 2]);
        title(['\eta at ',num2str(lambdas(ch)*1e9),'nm']);
        colorbar;
        if islogcoup; caxis([-4 0]); end
        colormap(gray(256));
    end
end
drawnow

end

%% Generate SPIE figure (first sample only)

% Use subtightplot for control over plot spacing
subplot = @(m,n,p) subtightplot (m, n, p, [0.12 0.06], [0.1 0.05], [0.05 0.05]);

% Define font and linewidth sizes
fontszAx = 12;
fontblAx = 'normal';
fontszTi = 14;

% Preallocate additional matrices
eta_ss = nan(wfSamps, 1);
eta_ps = nan(wfSamps, 1);
eta_sLs = nan(wfSamps, numWavelengths);
eta_pLs = nan(wfSamps, numWavelengths);
eta_pAvgs = nan(wfSamps, numWavelengths, N/2-1);

gifFig = figure(7);
set(gifFig, 'Units','normalized', 'OuterPosition', [0.05 0.05 0.9 0.9], 'Color', 'w');

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
isaxYlimEtaS = false;       % Flag to mark whether to hardcode y-axis
axYlimEtaS = [0 0.25];      % y axis limit
eta_sYSCL  = 1e3;           % scaling for y axis

%-- Width of planet coupling average region
prg = 1;        % Average over region +/-prg pixels wide

%-- Plot instantaneous WFE/R (at central lambda)
% Convert WFE from rads to nm
wfr = phzs/2/pi*lambda0*1e9;
% Plot
axWFR = subplot(2,3,1); 
datWFR = imagesc(xvals/apRad,yvals/apRad,squeeze(wfr(1,:,:)));
axis image; 
axis([-1 1 -1 1]);
set(axWFR, 'TickDir', 'out');
title(['Wavefront Residuals at ',num2str(lambdas(lam0)*1e9),'nm'], 'FontSize', fontszTi);
xlabel('x/D', 'FontWeight', fontblAx, 'FontSize', fontszAx)
ylabel('y/D', 'FontWeight', fontblAx, 'FontSize', fontszAx)
wfrCBAR = colorbar; 
wfrCBAR.Label.String = 'Residuals [nm]';
wfrCBAR.Label.FontWeight = fontblAx;
wfrCBAR.Label.FontSize = fontszAx;
colormap(parula(256));
caxis([min(wfr, [], 'all') max(wfr, [], 'all')])

%-- Plot instantaneous Donut PSF (polychromatic)
axPSF = subplot(2,3,2);
datPSF = imagesc(xvals/lambdaOverD,yvals/lambdaOverD,squeeze(iPSFv_BBs(1,:,:)));
axis image; 
set(axPSF, 'TickDir', 'out');
axis(axlimPSF);
title('Broadband PSF w/ vortex', 'FontSize', fontszTi);
xlabel('\lambda/D', 'FontWeight', fontblAx, 'FontSize', fontszAx)
ylabel('\lambda/D', 'FontWeight', fontblAx, 'FontSize', fontszAx)
colorbar;%caxis([-3 0])
colormap(parula(256));
caxis([min(iPSFv_BBs, [], 'all') max(iPSFv_BBs, [], 'all')])

%-- Plot eta_s vs. time [eta_s first since it's needed for eta_p and eta map]
% Get eta_s as min in center of all eta_maps averaged together
    % Average eta_maps since fiber can only be in one location so need to make
    % sure we take the eta_s at this same locaiton always.
crp = 68;   %1/2 of crop window for min analysis
eta_cent = squeeze(eta_mapss(1,round(end/2-end/crp):round(end/2+end/crp),round(end/2-end/crp):round(end/2+end/crp),:));
eta_cent = mean(eta_cent,3);
eta_s = min(eta_cent,[],[1 2]);    % Get min of average in central region as broadband null
% Get indices for null location
[eta_sInd(1), eta_sInd(2)] = find(eta_s==eta_cent);
eta_sInd = eta_sInd + round(N/2-N/crp-1);   % Account for cropping from before
% Eta_s vs wavelength (at eta_sInd location)
eta_sLs = squeeze(eta_mapss(:, eta_sInd(1), eta_sInd(2), :))*eta_sYSCL;
% Rescale broadband null value for plot
eta_ss(1) = mean(eta_sLs(1,:), 2);   
% Plot
axNull = subplot(2,3,5);
datNull = scatter(1:wfSamps, eta_ss);
if ~isaxYlimEtaS
    axYlimEtaS(2) = 1.2*max(eta_sLs(:));
end
ylim(axYlimEtaS)
xlim([0 wfSamps])
box on
ylabel('\eta_s [\times 10^{-3}]', 'FontWeight', fontblAx, 'FontSize', fontszAx)
xlabel('WF Sample', 'FontWeight', fontblAx, 'FontSize', fontszAx)
title('Broadband Star Coupling', 'FontSize', fontszTi)

%-- Plot eta_p vs. time [eta_p before eta map since needed for circle plot]
% Iterate through wavelengths to get azimuthal average of coupling
%eta_sI = nan(numWavelengths, 2);
for i = 1:numWavelengths
    % Azimuthally average
    [tmpAvg, qvec] = VFN_An_radAverage(squeeze(eta_mapss(1,:,:,i)), [N/2+1, N/2+1]);
    % Only populate part of matrix that is valid from tmpAvg
    eta_pAvgs(1,i,:) = tmpAvg;
end
% Average along wavelength dimension to find best average planet coupling
[~, eta_pInd] = max(mean(squeeze(eta_pAvgs(1,:,:)), 1));
eta_ps(1) = mean(eta_pAvgs(1,:,eta_pInd-prg:eta_pInd+prg), [2 3])*100;  % convert to [%]
% Plot
axPeak = subplot(2,3,4);
datPeak = scatter(1:wfSamps, eta_ps);
ylim(axYlimEtaP)
xlim([0 wfSamps])
box on
ylabel('\eta_p [%]', 'FontWeight', fontblAx, 'FontSize', fontszAx)
xlabel('WF Sample', 'FontWeight', fontblAx, 'FontSize', fontszAx)
title([sprintf('Broadband Planet Coupling\nfrom %3.1f', qvec(eta_pInd-prg)/lambdaOverD) '\lambda/D to ' num2str(qvec(eta_pInd+prg)/lambdaOverD,'%3.1f') '\lambda/D'], 'FontSize', fontszTi)

%-- Plot instantaneous eta map (at central lambda)
axEMap = subplot(2,3,3);
datEMap = imagesc(xvals/lambdaOverD,yvals/lambdaOverD,squeeze(eta_mapss(1,:,:,lam0)));
% Show circle where coupling was averaged. 
    % subtract -2 and xvals(end)... to center on image axes.
crc1 = viscircles(([N/2-1, N/2-1]/lambdaOverD-xvals(end)/lambdaOverD),  (eta_pInd-prg)/lambdaOverD, 'LineWidth', 1.0);
set(crc1.Children(2), 'Visible', 'off');   % Hide outline of circle
crc2 = viscircles(([N/2-1, N/2-1]/lambdaOverD-xvals(end)/lambdaOverD),  (eta_pInd+prg)/lambdaOverD, 'LineWidth', 1.0);
set(crc2.Children(2), 'Visible', 'off');   % Hide outline of circle
axis image;
axis(axlim2DETA);
set(axEMap, 'TickDir', 'out');
title(['\eta at ',num2str(lambdas(lam0)*1e9),'nm'], 'FontSize', fontszTi);
xlabel('\lambda/D', 'FontWeight', fontblAx, 'FontSize', fontszAx)
ylabel('\lambda/D', 'FontWeight', fontblAx, 'FontSize', fontszAx)
colorbar;
colormap(gray(256));
caxis([min(eta_mapss, [], 'all') max(eta_mapss, [], 'all')])

%-- Plot eta_p and eta_s vs. wavelength
% Eta_p vs wavelength (average w/ centering based on eta_sInd)
eta_pLs(1,:) = mean(eta_pAvgs(1,:,eta_pInd-1:eta_pInd+1),3)*100;
% Plot
axCoup = subplot(2,3,6);
title('Coupling Fractions Across Band', 'FontSize', fontszTi)
yyaxis left
datCoupSLin = plot(lambdas*1e9, eta_sLs(1,:), '-o');
%hold on
%datCoupSSca = scatter(lambdas*1e9, eta_sLs(1,:));
%hold off
ylim(axYlimEtaS)
ylabel(['\eta_s [\times 10^{-' num2str(log10(eta_sYSCL)) '}]'], 'FontWeight', fontblAx, 'FontSize', fontszAx)
xlabel('Wavelength [nm]', 'FontWeight', fontblAx, 'FontSize', fontszAx)
yyaxis right
datCoupPLin = plot(lambdas*1e9, eta_pLs(1,:), '-o');
%hold on
%datCoupPSca = scatter(lambdas*1e9, eta_pLs(1,:));
%hold off
ylim(axYlimEtaP)
ylabel('\eta_p [%]', 'FontWeight', fontblAx, 'FontSize', fontszAx)

%% Iterate through remaining samples, updating/saving figure each time

if isSaveGif
    % Save initial frame
    frame = getframe(gifFig); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    imwrite(imind,cm,[svfld svnm],'gif', 'Loopcount',inf);
end

% Average eta_mapss along wavelength dim. once before loop
eta_mapssMN = mean(eta_mapss, 4);

for i = 2:wfSamps
%-- Update WFE/R 
%imagesc(xvals/apRad, yvals/apRad, squeeze(wfr(i,:,:)), 'Parent', axWFR)
set(datWFR, 'CData', squeeze(wfr(i,:,:)));

%-- Update Donut PSF
%imagesc(xvals/lambdaOverD,yvals/lambdaOverD,squeeze(iPSFv_BBs(i,:,:)), 'Parent', axPSF);
set(datPSF, 'CData', squeeze(iPSFv_BBs(i,:,:)));

%-- Update eta_s
    % Note, keep fiber at same location to not "correct" for T/T
%eta_sLs(i,:) = squeeze(eta_mapss(i, eta_sInd(1), eta_sInd(2),:));
eta_ss(i) = mean(eta_sLs(i,:));
%scatter(1:wfSamps,eta_ss, 'Parent', axNull)
set(datNull, 'YData', eta_ss);

%-- Update eta_p
for ii = 1:numWavelengths
    % Azimuthally average
    tmpAvg = VFN_An_radAverage(squeeze(eta_mapss(i,:,:,ii)), [N/2+1, N/2+1]);
    % Only populate part of matrix that is valid from tmpAvg
    eta_pAvgs(i,ii,:) = tmpAvg;
end
eta_pLs(i,:) = mean(eta_pAvgs(i,:,eta_pInd-prg:eta_pInd+prg),3)*100;
eta_ps(i) = mean(eta_pLs(i,:));
%scatter(1:wfSamps, eta_ps, 'Parent', axPeak)
set(datPeak, 'YData', eta_ps);

%-- Update eta map
%imagesc(xvals/lambdaOverD,yvals/lambdaOverD,squeeze(eta_mapss(i,:,:,lam0)));
%viscircles(fliplr((eta_sInd-2)/lambdaOverD-max(xvals)/lambdaOverD),  (eta_pInd-1)/lambdaOverD);
set(datEMap, 'CData', squeeze(eta_mapss(i,:,:,lam0)));

%-- Update eta_p and eta_s vs. wavelength
%yyaxis(axCoup, 'left')
%plot(lambdas*1e9, eta_sLs(i,:), 'Parent', axCoup)
%hold on
%scatter(lambdas*1e9, eta_sLs(i,:), 'Parent', axCoup)
set(datCoupSLin, 'YData', eta_sLs(i,:));
%set(datCoupSSca, 'YData', eta_sLs(i,:));
%hold off
%yyaxis(axCoup, 'right')
%plot(lambdas*1e9, eta_pLs(i,:), 'Parent', axCoup)
%hold on
%scatter(lambdas*1e9, eta_pLs(i,:), 'Parent', axCoup)
%hold off
set(datCoupPLin, 'YData', eta_pLs(i,:));
%set(datCoupPSca, 'YData', eta_pLs(i,:));
drawnow

if isSaveGif
    % Capture the plot as an image 
    frame = getframe(gifFig); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    imwrite(imind,cm,[svfld svnm],'gif','WriteMode','append'); 
end
%pause(1)
end

%% Display wavefront coefficients
zrnFig = figure(8);
mrks = {'-o', '-s', '-d', '-^', '-*', '-x', '-p'};
legs = cell(length(nolls),1);
% plot(1:wfSamps, coeffs, mrks)
for i = 1:length(nolls)
    plot(coeffs(:,i), char(mrks(i)), 'MarkerSize', 4, 'LineWidth', 2);
    if i == 1
        hold on
    end
    legs(i) = {num2str(nolls(i))};
end
hold off

legH = legend(legs);
title(legH, 'Noll Index')

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
    *** DONE
    - have each plot in figure be its own matrix
    - save each frame as an itr axis in matrices
4) use rand() to generate random wavefronts and create time-dependent plots
5) create fully dynamic plot (ie. gif)
    *** DONE
6) Account for tip/tilt errors:
    *** DONE: keep fiber at same location from first null
    - By finding best null in each frame individually, we are basically
    correcting the fiber positioning infinitly well. In reality, tip/tilt errors
    will result in the fiber positioning being slightly incorrect. We need to
    account for this accordingly.
7) Calculate eta_p at specific location
    *** IGNORED SINCE GOING OFF CENTER OF FRAME
    - Add ability to calculate eta_p as average at a specific lam/D from the null
    point (or from wherever the fiber is centered).
    - This is because as is, we will state where the ideal planet coupling
    occurs but won't show what the planet coupling for a fixed planet source
    would be.
8) Average planet coupling over range of radii
    *** DONE
9) Center planet coupling calc. on frame, not null.
    *** DONE
10) Pass circum_diam/2 to generateZernike... INSTEAD of apRad
    *** DONE
11) Make figure and extract frame at each itr of first loop
    *** DONE
%}