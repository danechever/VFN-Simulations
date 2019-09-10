%{
Code for simulating the KPIC-VFN on-sky performance given known AO and other 
performance parameters. This script produces the wavefronts only.

This is based off the KPICVFN_OnSky_Simulator code.
    - This is the analysis section of that code.

Goals:___
- Simulate on-sky performance given synthetic or real WF residual data
- Account for atmospheric dispersion and ADC performance

Notes:___
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
* Synthetic WF residuals are modulated accross samples using sine function to 
    have continuously variable WF.
    ** [Updated] This is now optional.
* Magnitude of each WF residual across samples is plotted at the end.
    ** [Updated] This is only done when a synthetic WF is used
    ** [Updated] This is done for both synthetic and real wavefronts and is done
        at the beginning (for synthetic) or live during the loop (for real)
* "subtightplot" library is used to provide control of subplot location.
    * This is NOT a native matlab library; user must download from file exchange
* Full figure is only plotted once and then only the data is updated on each
    iteration of wavefront samples.
* Gif is optionally output at the end.
* For real WF data:
    - User must provide fits files with [xdim, ydim, sampleDim] format
    - User also provides a course mask for this data
        - That mask is only used for identifying which pixels/actuators are
        actually within the real beam.
    - The WF data is rotated and rescaled to match the sampling of our
        simulation. This is done using interpolation. I've tried to do these
        steps in such a way that the data is minimally affected.
    - Once the WF is oriented and scaled correctly, our simulated keck pupil
        (which includes central obscuration, spiders, fine outer pattern, etc.)
        is applied to the WF data. 
    ** [Updated] Fudge-factors have been added for user to tune the centering,
        orientation, and size of the simulated pupil within the WF data to
        ensure it is sampling the same region as the real pupil mask.
    ** [Updated] the raw WF (described above) can be used directly or a zernike
        decomp can be done and this reconstructed wavefront can be applied
        instead.
* Plotting is now done much more efficiently:
    - As before, main figure is only produced once and then updated.
    - However, instead of saving all frames of all the data, only the current
        frame's is saved and this gets overwritten on each WF sample iteration.
    - This vastly decreases RAM use and increases calculation speed.
* Code now handles tip/tilt residuals as well as ADC residulas
    - User can input a tip/tilt residual fits file to process real data. This
        file must be in [arcsec] and have dimensionality [wfSamps, ti[/tilt]
    - Synthetic residuals can also be provided in [mas] directly. These can be
        modulated in the code or kept constant.
    - ADC residuals are provided at each sample and wavelength.
%}

clear; close all; 
addpath('C:\Users\danie\Documents\MATLAB\VFN-Simulations\VFNlib');
addpath(genpath('C:\Users\danie\Documents\MATLAB\VFN-Lab\AnalysisCode\'))

%% Input parameters 

% Define smapling info
N = 2^12; % Size of computational grid (NxN samples) 
apRad = 128*2; % Aperture radius in samples 

% Define wavelength info
lambda0 = 2200e-9; %central wavelength
fracBW = 0.20; %\Delta\lambda/\lambda
numWavelengths = 5;% number of discrete wavelengths 
lambdas = getWavelengthVec(lambda0,fracBW,numWavelengths);% array of wavelengths (meters)

%--- Define wavefront residuals at the central wavelength
isrealWF = true;      % Flag to mark if real or synthetic WFs should be used
% Number of zernike modes used in reconstruction/analysis
recNMds = 36;
if isrealWF
    % Real WF data parameters
    %wfPTH = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\SPIE_2019\telemetry_20190420\';
    wfPTH = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\SPIE_2019\telemetry_20190617\';
        % NOTE: June data only has 1 minute's worth (59,000 samples)
    wfRSN = 'residualWF.fits';          % Filename for residuals
    wfMKN = 'DMPupil.fits';             % Filename for pupil mask on residuals
    % NOTE::: Raw data assumed to be in [Volts]
    % Define sample start, end, and step values
    wfstrt = 1;
    wfSamps = 118;
    wfStepSz = 500;
    % Scaling factor for data (nm/V)
    wfSCL = 600;     
    % Zernikes to plot in figure
    nolls = [2, 3, 5, 6, 7, 8];
    % Flag to choose if WF reconstruction or raw WF should be used
    israwWF = false;
    %%%% NOTE:::: DEAL WITH MANUAL CORRECTIONS IN WF SECTION BELOW!!!
else
    % Synthetic wavefront parameters
    nolls = [5];
    coeffs = 0.0*[0.01];    % Waves RMS at given nolls
    wfSamps = 10;                      % Number of WF samples to synthesize
    isWFPhzMod = false;                % Flag to chose if phase should be modulated
    wfModPhzEnd = 2*pi;                % Temporal modulation phase end (without delay)
end

%-- Define Tip/Tilt residuals  (ADC PARAMS ARE DEFINED SEPARATELY)
% NOTE::: To disable TT residuals, use isrealTT = false and set ttres = [0 0];
isrealTT = true;      % Flag to mark if real or synthetic TTs should be used
% Scaling factor for data (arcsec to waves PV at lambda0)
keckD = 10.949;                 % Real-world keck pupil diameter [m] - 10.949 = circumscribed
ttSCL = keckD/206265/lambda0;       % (D in [m])/(arcsec/rad)/(lambda in [m])  = arcsec2wavesPV
if isrealTT
    % Real TT data parameters
    %ttPTH = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\SPIE_2019\telemetry_20190420\';
    ttPTH = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\SPIE_2019\telemetry_20190617\';
    ttRSN = 'residualTT.fits';          % Filename for residuals
    % NOTE::: raw data assumed to be in [arcsec]
    % Define sample start values
    if exist('wfstrt', 'var')
        % By default, use wfstrt when defined
        ttstrt = wfstrt;
    else
        % NOTE::: If wfstrt is not defined, user must define a ttsrt manually
        ttstrt = 100;
    end 
    if exist('wfStepSz', 'var')
        % By default, use wfStepSz when defined
        ttStepSz = wfStepSz;
    else
        % NOTE::: If wfStepSz is not defined, user must define a ttStepSz manually
        ttStepSz = 1;
    end
    % NOTE::: no ttSamps since should sample exactly as many TTs as WFs     
else
    % Synthetic TT parameters
    ttres = 0*[41.4452, 0];    % Peak TT error in [mas] at central wavelength
    isTTPhzMod = true;     % Flag to chose if phase should be modulated
    ttModPhzEnd = 2*pi;     % Temporal modulation phase end (without delay) 
end

%-- Define ADC residuals:   in [mas]
% Flag to mark if real adc residuals should be used
isrealADC = true;
if isrealADC
    %-- Real ADC params
    adcPTH = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\SPIE_2019\ADC_residuals\';
    adcRSN = 'corrected_dispersion_K_z30.csv';
    %adcRSN = 'uncorr_dispersion_K_z30.csv';
    % NOTE::: files should be in csv, col1=wavelengths, col2=lambdas in [mas]. 
    %   Values are interpolated from this data. 
else
    %-- User-provided manual ADC values 
    % NOTE::: To disable ADC residuals, set adcres = to all 0's
    %adcres = 0*ones(2,numWavelengths);         % Disable ADC residuals
    % NOTE::: To define constant ADC residuals, provide 2D matrix here
    adcres = 0*[0 0 0 0 0;...         % Tip residuals at each wavelength in [mas]
             -10 -5 0 5 10];          % Tilt residuals at each wavelength in [mas]
    % Make into 3D matrix with dimensionaly: [wfSamps, tip/tilt, wavelength]
    adcres = permute(repmat(adcres, 1, 1, wfSamps),[3 1 2]);
    % NOTE::: To provide time-varyin ADC:
    % Optionally, provide 3D matrix directly with time-varying ADC residuals
        % Dimensionality still needs to match [wfSamps, tip/tilt, wavelength]
end

% Flag to plot intermediate figures (before SPIE figure)
isPlotSimp = true;

% Flag to display zernike figure
isPlotResCof = true;

%-- Saving parameters
% Flag to save gif
isSaveGif = false;
% Flag to save fits
isSaveFit = false;
% Save folder
svfld = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\SPIE_2019\KPIC_OnSkySims\res12_apRad256_AllOn_118Samps_JuneRunBLAH\';
if isSaveFit
    % Create foler if it doesn't already exist
    mkdir(svfld)
end
% Save name for gif
svnmGif = 'WFR.gif';
% Delay time for gif (time between frames in [s])
gifDelay = 0.3; 

%% Generate the coordinate system

coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
xvals = coords.xvals;% Helpful for plotting
yvals = coords.yvals;

%% Create array with pupil function

%PUPIL = makeCircularPupil( apRad, N );
[PUPIL, pupCircDiam] = makeKeckPupil( 2*apRad, N );
[normI, totalPower0] = getNormalization(PUPIL);% Normalization factors
lambdaOverD = N/(pupCircDiam/2)/2; % lam/D in units of samples in the image plane

% figure(1)
% imagesc(xvals/(pupCircDiam/2),yvals/(pupCircDiam/2),PUPIL);
% axis image; 
% axis([-1 1 -1 1]);
% title('Pupil');
% colorbar; 
% colormap(parula(256));
% drawnow;

%% Deal with WF extraction (real data) or synthesis (synthetic data)
if isrealWF
    %*** Real data
    
    %-- Read data 
    wfres = fitsread([wfPTH, wfRSN]);   % Full residuals file
    wfmsk = fitsread([wfPTH, wfMKN]);   % Mask file
    % Extract specific samples
    wfres = wfres(:,:,wfstrt:wfStepSz:wfstrt+wfStepSz*(wfSamps-1));
    % Rescale the WFR to waves at lambda0
    wfres = wfres*wfSCL*1e-9/lambda0;    %x1e-9 to go to [m] since lambda0 is in [m]
    % Copy raw mask for plotting later
    wfmskRAW = wfmsk;
    % Extract in-pupil pixels (raw file is 0=out, 1=in pupil, 2=in DM)
    wfmsk(wfmsk~=1) = 0;
    
    %-- Mask corrections (apply manually on case-by-case basis)
    % Fill in central pixels
    wfmsk(10:12, 9:10) = 1;
    % Correct shifted pixels (1 extra at top and one missing at bottom)
    wfmsk(2,13) = 0; wfmsk(19,13) = 1;
    % Fill in closest pixels (causing simulated mask to be severely undersized)
    wfmsk(15,17) = 1; wfmsk(17, 16) = 1; wfmsk(4, 15) = 1;
    
    %-- Aditional nobs to turn for mask features:
    % rot: angle [degrees] by which to rotate data to match our pupil
    % sclfdg: Fudge factor for mask resizing to match simulated pupil to wfmsk
    % cntfdg: Fudge factor for mask centering to match sim pup to wfmsk
    
    %-- Pre-process WF before loop
    % Mask WF using charlotte's pupil
    wfmsk = repmat(wfmsk, [1 1 wfSamps]);
    wfres = wfres.*wfmsk;
    % Pad sides of WF so that up-sampling produces even-sized grid
    wfres = padarray(wfres, [1 1 0]);
    
else
    %*** Synthetic data
    % Create matrix-version of coeffs to simulate changing WFE
    wfMod = ones(wfSamps, length(nolls));
    if isWFPhzMod
        for i = 1:length(nolls)
            % Create random phase delay
            phzShift = rand*2*pi;
            wfMod(:, i) = sin(linspace(0+phzShift, wfModPhzEnd+phzShift, wfSamps))';    % smoothly modulate WFE amplitude
        end
    end
    %wfMod = repmat(wfMod, 1, length(coeffs));
    coeffs = repmat(coeffs, wfSamps, 1);
    coeffs = coeffs.*wfMod;             % Apply modulation
end

%% Display WF Mask (for real data) or Zernikes (for synthetic data)
if isrealWF
    if isPlotSimp
        % Show raw mask
        figure(101);
        imagesc(wfmskRAW);
        axis image;
        title('Original WF Mask - from PyWFS')
        % Show corrected mask
        figure(102);
        imagesc(wfmsk(:,:,1))
        axis image;
        title('User-Corrected WF Mask - from PyWFS')
    end
    
    % NOTE: Zernikes for real data are calculated during main loop
    
else
    if isPlotResCof
        % Show synthetic coefficients
        zrnFig = figure(100);
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
        title('Zernike Decomposition')
        xlabel('WF Sample')
        ylabel('RMS [waves at \lambda_0]')
        legH = legend(legs);
        title(legH, 'Noll Index')
    end
end    

%% Deal with TT extraction (real data) or synthesis (synthetic data)
if isrealTT
    %*** Real data
    
    %-- Read data 
    ttres = fitsread([ttPTH, ttRSN]);   % Full residuals file
    % Extract specific samples
    ttres = ttres(ttstrt:ttStepSz:ttstrt+ttStepSz*(wfSamps-1),:);
    % Rescale the TTR to waves PV at lambda0
    ttres = ttres*ttSCL;    
       
else
    %*** Synthetic data
    
    %-- Create matrix-version of coeffs to simulate changing TT
    ttMod = ones(wfSamps, length(ttres));
    if isTTPhzMod 
        % Modulate phase as requested by user
        for i = 1:length(ttres)
            % Create random phase delay
            phzShift = rand*2*pi;
            ttMod(:, i) = sin(linspace(0+phzShift, ttModPhzEnd+phzShift, wfSamps))';    % smoothly modulate TT amplitude
        end
    end
    ttres = repmat(ttres, wfSamps, 1);
    ttres = ttres.*ttMod;             % Apply modulation
    %-- Convert from mas to [waves PV at lambda0]
    ttres = ttres*1e-3*ttSCL;        % *1e-3 to go from mas to arcsec for ttSCL
end

%% Display TT residuals (no atmospheric dispersion yet)
if isPlotResCof
    ttFig = figure(104);
    % Set axes colors to black [0 0 0]
    set(ttFig,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
    yyaxis right
    % Plot right axis [waves RMS]
    ttPLT = plot(1:wfSamps, ttres(:,1), 'r-o', 1:wfSamps, ttres(:,2), 'b-o', 'MarkerSize', 4, 'LineWidth', 2);
    title('Tip Tilt Residuals - Before Atm. Disp.')
    xlabel('WF Sample')
    legend('Tip', 'Tilt');
    ylabel('RMS [waves at \lambda_0]')
    % Save ylimits for rescaling on left side
    ttYlimL = ylim;
    % Add left label [mas]
    yyaxis left
    ylabel('RMS [mas]')
    % Rescale using mean of ttSCL to account for difference in Tip/Tilt PV2RMS
    ylim(ttYlimL*1e3/ttSCL);
end

%% Deal with atmospheric dispersion
if isrealADC
    %-- Load real data
    % Read file
    adcraw = csvread([adcPTH adcRSN], 1, 0);    % ignore first row with names
    % Interpolate dispersion for given wavelengths
    adcres = interp1(adcraw(:,1),adcraw(:,2),lambdas*1e6, 'linear', 'extrap');
    % Pad to inlcude 0's for non-dispersion axis
    adcres = padarray(adcres, 1, 0, 'pre');
    % Make into 3D matrix with dimensionaly: [wfSamps, tip/tilt, wavelength]
    adcres = permute(repmat(adcres, 1, 1, wfSamps),[3 1 2]);
end

% Convert from [mas] to [waves at lambda0]
adcres = adcres*1e-3*ttSCL;

% Convert ttres to 3D to include wavelengths. New dim.: [wfSamps, tip/tilt, wavelength];
ttres = repmat(ttres, 1, 1, numWavelengths);

% Add in ADC residuals as additional TT errors directly
ttres = ttres+adcres;

%% Save fits with ttres for star centering in figure-generating script
if isSaveFit
    fitswrite(ttres, [svfld 'ttres.fits']);
end

%% Generate unit zernike modes (once before loop since does not change)
    % Generating unit modes once improves runtime by 25% on top of fiber mode
    % improvement for N=2048 and wfSamps=10

% Define zernikes on which to decompose
if (isrealWF && israwWF)
    %-- Reconstruction not requested; use raw WF
    % Set loop to only decompose requested modes (nolls)
    jvls = nolls;
else
    %-- Use reconstruction for WF instead of raw values
    % Set loop to decompose on all recNMds
    jvls = 1:recNMds;
end

% Preallocate zernike matrix
unitzrn = nan(N, N, length(jvls));
% Create basis
for j = jvls
    % Generate zernike (on keck pupil) with 1 wave RMS error
    tmpzrn = generateZernike_fromList(j, 1, PUPIL, pupCircDiam/2, coords);
    % Rescale zernike from radians to waves
    unitzrn(:,:,j) = tmpzrn/(2*pi);
end

%% Loop through WF samples
if isrealWF
    %-- Preallocate realWF zernike coeffs FOR PLOTTING
        % NOTE: WF reconstruction is still done over 36 zernikes
    coeffs = nan(wfSamps, length(nolls));
end

%-- Define font and linewidth sizes
fontszAx = 14;
fontblAx = 'normal';
fontszTi = 14;

%-- Axis limits
% realWF Zernike plot
zrnSCL = 1e2;               % scaling factor for y axis
axYlimZrn = [-10 10];          % y axis limits (in units of 1e3)

%-- Coloraxes
% WFR map
if isrealWF
    % Real WF thus get actual value directly from full matrix
        % slightly inaccurate since different mask and interp but close enough
    %caxlimWFR = [min(wfres(:)*lambda0*1e9) max(wfres(:)*lambda0*1e9)];
    caxlimWFR = [-400 400];
else
    % Approximate PV based on the RMS values from the coeffs
    caxlimWFR = [3*min(coeffs(:)*lambda0*1e9) 3*max(coeffs(:)*lambda0*1e9)];
end

%-- Crop region for fits files
wfrCrp = [N/2-round(pupCircDiam/2)-20:N/2+round(pupCircDiam/2)+20];

for i = 1:wfSamps
fprintf('wfSamp %06d of %06d\n', i, wfSamps);
    
%% Define pupil field
%*** DEAL WITH WAVEFRONT RESIDUALS
if isrealWF  
    
    %-- FORMAT REAL WAVEFRONT FOR SIMULATION
    % Extract frame
    wfTMP = wfres(:,:,i);
    % Up-sample to decrease error sensitivity in centering and rotation
        % Want enough pts that rot. interp. via nearest-neighbor looks good.
    wfTMP = interp2(wfTMP, 5, 'nearest');
    % Recenter and crop (assuming roughly symmetric)
    [rws, cls] = find(wfTMP);  % find min and max row/col indices
    cntfdg = [0 -10];   % centering fudge factor  [row col] ** must be even
    rw1 = min(rws); rw2 = max(rws);
    cl1 = min(cls); cl2 = max(cls);
    wfTMP = wfTMP(rw1-1+cntfdg(1)/2:rw2+cntfdg(1)/2, cl1-1+cntfdg(2)/2:cl2+cntfdg(2)/2);
    % Rotate to fix orientation
    rot = 27;       % 27 provides decent alignment for furthest segments
    wfTMP = imrotate(wfTMP, rot);
    % Find inscribed circle size
    [rws, cls] = find(~wfTMP);
    crc = 2*min(sqrt((rws-size(wfTMP,1)/2-1).^2 + (cls-size(wfTMP,2)/2-1).^2));
    % OPTIONAL: show inscribed circle
    %figure; imagesc(wfTMP)
    %viscircles([size(wfTMP,2)/2-1, size(wfTMP,1)/2-1], crc/2);
    %title('Inscribed Circ B4 Resize')
    % Resize so that simulated pupil closely matches real pupil (wfmsk)
    sclfdg = -23;    % fudge factor for pupil resizing ** (-)=grow mask; (+)=shrink mask
	x = 0:(pupCircDiam+sclfdg)/crc:(size(wfTMP,2)-1)*(pupCircDiam+sclfdg)/crc;
    y = 0:(pupCircDiam+sclfdg)/crc:(size(wfTMP,1)-1)*(pupCircDiam+sclfdg)/crc;
    xi = 0:(size(wfTMP,2)-1)*(pupCircDiam+sclfdg)/crc;
    yi = 0:(size(wfTMP,1)-1)*(pupCircDiam+sclfdg)/crc;
    [xt, yt] = meshgrid(x,y); %create vector arrays
    [xit, yit] = meshgrid(xi,yi); %create vector arrays
    wfTMP = interp2(xt,yt,wfTMP,xit,yit);
    % Pad to match Simulated pupil diameter
    pd = size(wfTMP) - [N N];
    if pd(1) < 0
        % Need to add rows 
        wfTMP = padarray(wfTMP, [-round(pd(1)/2) 0]);
    end
    if pd(2) < 0
        % Need to add columns
        wfTMP = padarray(wfTMP, [0 -round(pd(2)/2)]);
    end
    % Final crop to ensure proper size 
        % Accounts for initial oversize or padding which adds 1 extra row/col
    wfTMP = wfTMP(1:N, 1:N);
    % Apply pupil!  (Not actually, pupil will be applied at Epup calc)
    % wfTMP = wfTMP.*PUPIL;
    % OPTIONAL: show final pupil with WF
    %figure; imagesc(wfTMP.*PUPIL); title('Final Raw WF - Simulated Pupil')
    
    %-- DECOMPOSE WAVEFRONT INTO ZERNIKES;
    recCof = nan(recNMds,1);
    if ~(isrealWF && israwWF)
        % Preallocate reconstructed WF matrix
        recWF = zeros(N);
    end
    for j = jvls
        % Extract current mode
        tmpzrn = unitzrn(:,:,j);
        % Left-divide to get coefficient
        recCof(j) = tmpzrn(logical(PUPIL))\wfTMP(logical(PUPIL));
        if ~(isrealWF && israwWF)
            recWF = recWF + tmpzrn*recCof(j);
        end
    end
    coeffs(i,:) = recCof(nolls);
    
    if isPlotResCof
        %-- DISPLAY REQUESTED ZERNIKE COEFFICIENTS (nolls)
        if i == 1
            % Show requested coefficients
            zrnFig = figure(100);
            legs = cell(length(nolls),1);
            for j = 1:length(nolls)
                cofPlt(j) = plot(coeffs(:,j)*zrnSCL, '-o', 'MarkerSize', 4, 'LineWidth', 2);
                if j == 1
                    hold on
                end
                legs(j) = {num2str(nolls(j))};
            end
            hold off
            xlim([0 wfSamps]);
            ylim(axYlimZrn);
            title('Zernike Decomposition (Approx in Pupil)', 'FontSize', fontszTi);
            xlabel('WF Sample')
            ylabel('RMS [% waves at \lambda_0]')% \times 10^{-2}]')
            legH = legend(legs);
            title(legH, 'Noll Index')
        else
            for j = 1:length(nolls)
                set(cofPlt(j), 'YData', coeffs(:,j)*zrnSCL);
            end
        end
    end
    
    %-- SET WF FOR USE AND CONVERT FROM WAVES TO RADIANS
    if (isrealWF && israwWF)
        wfphz = wfTMP*2*pi;
    else
        % Renormalize one last time to ensure RMS values match.
        rmsWF = sqrt(mean(wfTMP(logical(PUPIL)).^2));
        recWF = recWF*rmsWF/sqrt(mean(recWF(logical(PUPIL)).^2));
        %recWF = recWF*(recWF(logical(PUPIL))\wfTMP(logical(PUPIL)));
        wfphz = recWF*2*pi;
    end
    
    %-- DISPLAY PUPIL PROJECTIONS
    if isPlotSimp && (i == 1)
        crcPup = makeCircularPupil(pupCircDiam/2, N);
        wfPLT = wfTMP; wfPLT(wfPLT==0) = nan;
        figure(103); 
        imagesc(xvals/(pupCircDiam/2), yvals/(pupCircDiam/2), wfPLT+0.05*PUPIL+0.1*crcPup); 
        axis image
        axis([-1.2 1.2 -1.2 1.2])
        title('Pupil Projections')
    end
    
else
    %-- Apply synthetic wavefronts
    wfphz = generateZernike_fromList( nolls, coeffs(i,:), PUPIL, pupCircDiam/2, coords); 
end

for ch = 1:numWavelengths
    %*** DEAL WITH TT AND ADC RESIDUALS (in loop to do by wavelength)
    % Check if there are any non-zero values in ttres in case user has disabled 
        % TT/ADC residuals via setting all to 0
    if find(ttres(:,:,ch))
        %-- Non-zero value found thus apply TT/ADC residuals
        % Build TT/ADC wavefronts
        ttphz = 2*pi*ttres(i,1,ch)*lambdaOverD*coords.X/N;  % Create X tilt
        ttphz = ttphz + 2*pi*ttres(i,2,ch)*lambdaOverD*coords.Y/N;  % Add Y tilt
    else
        ttphz = zeros(N, N);
    end
    
    %*** COMBINE ALL PHASE ERRORS 
    pupphz(:,:,ch) = (ttphz + wfphz)*lambda0/lambdas(ch);   % also rescale to lambda
end

% Convert WFE from rads to nm
wfr = wfphz/2/pi*lambda0*1e9.*PUPIL;      %.*PUPIL Optional to show with pupil
% Set out-of-pupil values to NaN for plotting
wfr(wfr==0) = nan;

%% Plot if desired
if isPlotSimp
    %-- Display WFR part of error only (top left plot in SPIE figure)
    if i == 1
        %-- First iteration: make the figure
        fig = figure(7);
        set(fig, 'units', 'inches', 'Position', [0 0 7 6], 'color', 'w');
        % Plot
        datWFR = imagesc(xvals/pupCircDiam,yvals/pupCircDiam,wfr);  
        axis image; 
        axis([-0.5 0.5 -0.5 0.5]);
        set(datWFR.Parent, 'TickDir', 'out');
        title(['Wavefront Residuals at ',num2str(lambda0*1e9),'nm'], 'FontSize', fontszTi);
        xlabel('x/D', 'FontWeight', fontblAx, 'FontSize', fontszAx)
        ylabel('y/D', 'FontWeight', fontblAx, 'FontSize', fontszAx)
        wfrCBAR = colorbar; 
        wfrCBAR.Label.String = 'Residuals [nm]';
        wfrCBAR.Label.FontWeight = fontblAx;
        wfrCBAR.Label.FontSize = fontszAx;
        %set(axWFR, 'Colormap', parula(256));
        colormap(datWFR.Parent, parula(256));
        caxis(caxlimWFR);
    else
        %-- Other iterations: update figure
        set(datWFR, 'CData', wfr);           
    end

    drawnow;
end

%% Save Results
if isSaveFit
    %-- Save WFR part of data
    fitswrite(wfr(wfrCrp, wfrCrp), [svfld sprintf('wfr%06d.fits',i)]);
    %-- Save combined pupil phase
    fitswrite(pupphz(wfrCrp, wfrCrp,:), [svfld sprintf('pupphz%06d.fits',i)]);
end

if isSaveGif && isPlotSimp
    % Save initial frame
    frame = getframe(fig); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if i == 1
        imwrite(imind,cm,[svfld svnmGif],'gif', 'Loopcount',inf, 'DelayTime', gifDelay);
    else
        imwrite(imind,cm,[svfld svnmGif],'gif','WriteMode','append','DelayTime', gifDelay); 
    end
end

end

if isSaveGif && isPlotSimp
    saveas(fig, [svfld 'WFR_FinalFrame.png'])
end