%{
Code for analyzing the WFE reported by the PyWFS on KPIC. It takes the raw WF
provided by the PyWFS, cleans it up, and outputs the: RMS WFE, zernike
coefficients, and cleaned-up but still pixelated WF.

*** Based heavily off KPICVFN_OnSky_WFGenerator.m

*** This script is meant to be used with KPICVFN_OnSky_NullCalculator

Goals:___
- Process raw WF residuals from PyWFS to provide useful data for analysis
- This data will be utilized for determining a simple equation relating RMS WFE
    to null depth.

Notes:___
* Does not do anything with T/T or ADC residuals at the moment.
* Alginment of PyWFS data with Keck PUPIL is done differently than before
    - PUPIL is rotated instead of PyWFS data
    - Resampling is done slightly differently
    - These mods are done once before the loop to minimize runtime.
    * mods are explained in more detail in pertinent section
* The PUPIL and most other matrices are cropped to the central region only
    - This should reduce memory requirements as well as runtime
    - Does not affect results here since no FFT are performed
* There is an option to decompose the WF into Zernikes using Charlotte's pupil
* Did a runtime profiler to check performance:
    - Took 6.5min to analyze 590 WF's
    - Costliest chunks of code were:
      - Creating Keck PUPIL: 57 s
      - Creating gif (getframe()): 73s
* When handling lots of data, RAM becomes an issue. Thus wfres is saved to HD
    - To reduce RAM, wfres is broken into chunks and processed in pieces
    - Each piece is pre-processed then saved as a fits which is then reloaded
        independently in the main loop
    * THUS, some cubes will always need to be saved and a save path should be
        provided
[01/2021]
* Modified to work on linux
    - Use "filesep" and relative pathing (VFN-Lab and falco-matlab should
        be in same higher-level directory)

Output:___
	- Zernike Coeffs (recCof but full vector at every instance) [2D matrix: row=time instance, col=zernike coeff]
		This is the zernikes projected on the Keck PUPIL, NOT charlotte's pupil
		* These are before the rmsWF renormalization!! Can renormalize later using rms_WF if needed *
		'zern_coeff%06d.fits'
	- Zernike Coeffs (on CHARLOTTE's pupil) [2D matrix: row=time, col=zern coeff]
		'CH_zern_coeff.fits'
	- RMS WFE at every instance (rmsWF) [2D matrix: row=time, col1=RMS of original (ln 256), col2=RMS of reformatted(ln 575)]
		rmsWF(:,1) is the rms within Charlotte's pupil
		rmsWF(:,2) is the rms within our Keck pupil
		'rms_WF%06d.fits'
	- raw wavefront reformatted (recentered/rescaled/sampled to Keck PUPIL size) [3d cube: rowXcol=2D wavefront, pane=time instance]
		'wfr_reform.fits'
		Header of this file has run parameters
	- reconstructed WF (reconstruction from zernikes, with RMS normalization) [3D cube: rowXcol=2D wavefront, pane=time instance]
		'wfr_recons%06d.fits'
	- pupils used (reformatted to final shape) [3D cube: rowXcol=2D wavefront, pane1=wfmskRaw pane2=wfmsk pane3=PUPIL]
		These are the all the pupils in a comparable form. Thus charlotte's pupil and my corrections to it are here with 
		  the same processing as the WF residuals so that they are at the same scale as PUPIL. This is also PUPIL after
		  rotation and cropping so that it's exactly as projected onto the final data.
		'wfr_pupils.fits'
%}

clear; close all; 
addpath(['..' filesep 'VFNlib']);
addpath(genpath(['..' filesep '..' filesep 'VFN-Lab' filesep 'AnalysisCode']))
addpath(genpath(['..' filesep '..' filesep 'falco-matlab']))

%% Input parameters 

% Define smapling info
N = 2^9; % Size of computational grid (NxN samples) 
apRad = N/2-4; % Aperture radius in samples 

% Define wavelength info
lambda0 = 2200e-9; %central wavelength

%--- Define wavefront residuals at the central wavelength
% Number of zernike modes used in reconstruction/analysis
recNMds = 36;
% Real WF data parameters
%wfPTH = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\SPIE_2019\telemetry_20190420\';
wfPTH = '/media/Data_Drive/VFN/KPIC_PyWFS_Data/telemetry_20190617/';
    % NOTE: June data only has 1 minute's worth (59,000 samples)
wfRSN = 'residualWF.fits';          % Filename for residuals
wfMKN = 'DMPupil.fits';             % Filename for pupil mask on residuals
% NOTE::: Raw data assumed to be in [Volts]
% Define sample start, end, and step values
wfstrt = 1;
wfSamps = 5900;
wfStepSz = 10;
% Scaling factor for data (nm/V)
wfSCL = 600;     
% Zernikes to plot in figure
nolls = 1:recNMds;
% Keck PUPIL Rotation value
PUP_rot = 30;
%%%% NOTE:::: DEAL WITH MANUAL CORRECTIONS IN WF SECTION BELOW!!!

%-- Control Flags
% Flag to decompose into Zernikes using Charlotte's pupil
    % This decomp is in addition to the normal decomp and processing
isDecompCH = true;
% Flag to plot intermediate figures (before SPIE figure)
isPlotSimp = true;
% Flag to display zernike figure
isPlotResCof = true;
% Flag to make figures visible
visflg = 'off';     % should be 'on' or 'off'
% Flag to make all figures visible at the end
isFigsVisEnd = true;


%-- Saving parameters
% Flag to save gif
isSaveGif = true;
% Flag to save fits of WF
    % NOTE: The wfr_reform will ALWAYS be saved to reduce RAM (see header)
isSaveFit = true;
% Flag to save rmsWF and Zernike Coeffs (including zerns on charlotte pupil)
isSaveRMS_COF = true;
% Checkpoints for save data
    % Number of frames to save per cube of reconstructed WF
    % Also number of time instances to save per cube of rmsWF and Zernike coeffs
SAVE_CHK = 590;
% Save folder
svfld = '/media/Data_Drive/VFN/KPIC_PyWFS_Data/telemetry_20190617/Processed/FullSet/';
if isSaveFit || isSaveGif || isSaveRMS_COF
    % Create foler if it doesn't already exist
    mkdir(svfld)
end
% Save name for gif
svnmGif = 'WFR.gif';
% Delay time for gif (time between frames in [s])
gifDelay = 0.1; 

%% Generate the coordinate system
fprintf('Generating Coords\n')
coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
xvals = coords.xvals;% Helpful for plotting
yvals = coords.yvals;

%% Create array with pupil function
fprintf('Generating PUPIL\n')
%PUPIL = makeCircularPupil( apRad, N );
[PUPIL, pupCircDiam] = makeKeckPupil( 2*apRad, N );

% Rotate pupil to match orientation of PyWFS data
PUPIL = imrotate(PUPIL, PUP_rot, 'crop');    % 'crop' maintains size and centering

% figure(1)
% imagesc(xvals/(pupCircDiam/2),yvals/(pupCircDiam/2),PUPIL);
% axis image; 
% axis([-1 1 -1 1]);
% title('Pupil');
% colorbar; 
% colormap(parula(256));
% drawnow;

%% Prepare plotting params

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
% Real WF thus get actual value directly from full matrix
    % slightly inaccurate since different mask and interp but close enough
%caxlimWFR = [min(wfres(:)*lambda0*1e9) max(wfres(:)*lambda0*1e9)];
caxlimWFR = [-400 400];

%% Load WF Data (real data)
fprintf('Loading WF data\n')
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

%% (OPTIONAL) calculate zernike's on charlotte's pupil
% This section is unoptimized in runtime/memory since the matrix is still very
    % small at this point so cost is low
if isDecompCH
    fprintf('Calculating Zernikes on Charlotte PUPIL\n')
    %-- Create copy of wfres to manipulate
    wfresCH = wfres;
    
    %-- Crop to center on Charlotte's pupil (does not remove anything in mask)
    wfTMP = wfmsk;
    [rws,cls] = find(wfTMP);
    rw1 = min(rws);rw2 = max(rws);
    cl1 = min(cls); cl2 = max(cls);
    cntfdg = [0 0];   % centering fudge factor  [row col] ** must be even
    cnt = [rw1+cntfdg(1)/2, rw2+cntfdg(1)/2,   cl1+cntfdg(2)/2, cl2+cntfdg(2)/2];
    wfresCH = wfresCH(cnt(1):cnt(2),cnt(3):cnt(4),:);
    wfmskCH = wfmsk(cnt(1):cnt(2),cnt(3):cnt(4));
    
    %-- MANUALLY pad sides to make matrices square and recenter on coords
    % Look at the resulting coordsCH.RHO matrix to try to also center the pupil
        % on the coordinate system this way as well. This is done by making sure
        % the edges of coordsCH.RHO align roughly with edge of this pupil. Edge
        % is defined by the apradCH value provided below
    wfresCH = padarray(wfresCH, [1 0 0]);
    wfmskCH = padarray(wfmskCH, [1 0 0]);
    wfresCH = padarray(wfresCH, [0 1 0], 'pre'); 
    wfmskCH = padarray(wfmskCH, [0 1 0], 'pre');
    apradCH = 9;
    
    %-- Apply mask to data
        % No repmat needed; 3D.*2D works as long as matrices are "compatible"
    wfresCH = wfresCH.*wfmskCH;
    
    %-- Create coords for these dimensions
    coordsCH = generateCoordinates(size(wfresCH,1));% Creates NxN arrays with coordinates 
    
    %-- Create unit Zernikes on Charlotte's pupil
    % Define zernikes on which to decompose
    % Set loop to decompose on all recNMds
    jvls = 1:recNMds;
    % Preallocate zernike matrix
    unitzrn = nan(size(wfresCH,1), size(wfresCH,2), length(jvls));
    % Create basis
    for j = jvls
        % Generate zernike (on keck pupil) with 1 WAVE RMS error
        unitzrn(:,:,j) = (1/(2*pi))*generateZernike_fromList(j, 1, wfmskCH, apradCH, coordsCH);
    end
    if any(isnan(unitzrn),'all')
        error('an element in unitzrn was left unallocated')
    end
    
    %-- Decompose
    % Preallocate result matrix
    coeffsCH = nan(size(wfresCH,3), length(nolls));
    % Iterate through WF frames decomposing
    for i = 1:size(wfresCH,3)
        wfTMP = wfresCH(:,:,i);
        for j = jvls
            % Extract current mode
            tmpzrn = unitzrn(:,:,j);
            % Left-divide to get coefficient
            coeffsCH(i,j) = tmpzrn(logical(wfmskCH))\wfTMP(logical(wfmskCH));
        end
    end
    if any(isnan(coeffsCH),'all')
        error('an element in coeffsCH was left unallocated')
    end
    
    %-- Save results
    if isSaveRMS_COF
        fitswrite(coeffsCH, [svfld 'CH_zern_coeff.fits']);
    end
    
    if isPlotResCof
    %-- DISPLAY REQUESTED ZERNIKE COEFFICIENTS (nolls)
        % Show requested coefficients
        zrnFigCH = figure(105);
        set(zrnFigCH, 'color', 'w', 'Visible', visflg);
        legs = cell(length(nolls),1);
        for j = 1:length(nolls)
            cofPltCH(j) = plot(coeffsCH(:,j)*zrnSCL, '-o', 'MarkerSize', 4, 'LineWidth', 2);
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
        legH = legend(legs, 'Location', 'eastoutside');
        title(legH, 'Noll Index')
    end
    
    %-- Cleanup workspace
    clear unitzrn coordsCH wfresCH wfTMP rws cls rw1 cl1 rw2 cl2 cntfdg cnt 
    clear wfmskCH apradCH jvls legs
    
end

%% Process WF Data as usual
fprintf('Finishing WF formatting\n')

%-- Calculate RMS WFE of raw WF with charlotte's mask
wfmsk(wfmsk == 0) = nan;
% RMS along axes [1 2] using nanmean to omit 0's. Squeeze to get 1D vector
rmsWF_ch = squeeze(sqrt(nanmean((wfres.*wfmsk).^2,[1 2])));

%-- Set wfmsk nan's back to 0's
wfmsk(isnan(wfmsk)) = 0;

%-- Mask corrections (apply manually on case-by-case basis)
% Fill in central pixels
wfmsk(10:12, 9:10) = 1;
% Correct shifted pixels (1 extra at top and one missing at bottom)
wfmsk(2,13) = 0; wfmsk(19,13) = 1;
% Fill in closest pixels (causing simulated mask to be severely undersized)
wfmsk(15,17) = 1; wfmsk([5,17], 16) = 1; wfmsk(3:4, 15) = 1; wfmsk(12,19) = 1; wfmsk(14,18) = 1; wfmsk(19,[5,15])=1;

%-- Aditional nobs to turn for mask features:
% rot: angle [degrees] by which to rotate data to match our pupil
% sclfdg: Fudge factor for mask resizing to match simulated pupil to wfmsk
% cntfdg: Fudge factor for mask centering to match sim pup to wfmsk

%-- Pre-process WF before loop
% Mask WF using charlotte's pupil (with corrections by user)
wfres = wfres.*wfmsk;
% Pad sides of WF so that up-sampling produces even-sized grid
wfres = padarray(wfres, [1 1 0]);

%% Display WF Mask (for real data) 
if isPlotSimp
    % Show raw mask
    Rmskfig = figure(101);
    set(Rmskfig, 'color', 'w', 'Visible', visflg);
    imagesc(wfmskRAW);
    axis image;
    title('Original WF Mask - from PyWFS')
    % Show corrected mask
    Cmskfig = figure(102);
    set(Cmskfig, 'color', 'w', 'Visible', visflg);
    imagesc(wfmsk)
    axis image;
    title('User-Corrected WF Mask - from PyWFS')
end
% NOTE: Zernikes for real data are calculated during main loop

%% Generate unit zernike modes (once before loop since does not change)
fprintf('Generating unit Zernikes\n')

% Define zernikes on which to decompose
% Set loop to decompose on all recNMds
jvls = 1:recNMds;

% Preallocate zernike matrix
unitzrn = nan(N, N, length(jvls));
% Create basis
for j = jvls
    % Generate zernike (on keck pupil) with 1 WAVE RMS error
    unitzrn(:,:,j) = (1/(2*pi))*generateZernike_fromList(j, 1, PUPIL, pupCircDiam/2, coords);
end

if any(isnan(unitzrn),'all')
    error('an element in unitzrn was left unallocated')
end

%% Recenter + Resample PyWFS data to match simulated PUPIL
% This is based off what was done in _WFGenerator except:
% - Manipulations to reformat the pupils are done outside the for loop
%   - This should decrease runtime since we're not repeating calcs unnecessarily
% - The keck PUPIL is rotated to match the PyWFS pupil
%   - This should simplify things since we're rotating 1 2D matrix instead of
%       rotating the PyWFS data every time frame by frame 
% - Crop bounds are handled more carefully in this code

%- Add in non-Keck PUPIL masks to process these in same way as raw data
    % This is to accurately compare these masks to the final PUPIL mask
if isSaveFit
    % Append pupils
    wfmskRAW = padarray(wfmskRAW, [1 1 0]);
    wfmsk = padarray(wfmsk, [1 1 0]);
    wfres = cat(3,wfres,wfmskRAW,wfmsk);
end

% Optional variable save for when manually tuning the fudge factors below
% wfres_bak = wfres;

% Optional variable recovery for when manually tuning fudge factors below
% wfres = wfres_bak;

fprintf('Reformatting wfr\n')


%-- Split wavefront reformatting into sub-chunks to reduce memory usage
    % Otherwise runtime skyrockets since we become RAM limited

%NOTE: indexing in the following section gets a bit complicated due to the fact
%that we append the pupils to the end of the cube as well as the checkpoints.
%The index format was determined using the following testcode.
% Nitr = ( floor(wfSamps/SAVE_CHK) + ceil(mod(wfSamps,SAVE_CHK)/SAVE_CHK) )
% for ii = 1: Nitr
% cbstr = (ii-1)*SAVE_CHK+1;    % start index for this cube
% if (ii == Nitr)   % on last iteration, bring in full cube (including pupils)
% cbendin = size(wfres,3);  % Index used to pull rest of cube
% cbendsv = wfSamps;        % Index used to only save wfres (not pupils)
% else              % on every other iteration, just pull a normal cube piece
% cbendin = ii*500;         % Index used to pull current cube chunk
% cbendsv = cbendin;        % Index to save just the current cube chunk
% end
% fprintf('ii: %2d  |  cbstr: %5d  |  cbendin: %5d  |  cbendsv = %5d\n', ii, cbstr, cbendin, cbendsv)
% if cbendin-cbendsv        % deal with pupils on last iter
% fprintf('extracting last %d frames\n', cbendin-cbendsv)
% end
% end

% Rename full cube so I don't need to go through and fix variable names below
wfresFULL = wfres;

% Calculate how many total output cubes we'll need
Nitr = ( floor(wfSamps/SAVE_CHK) + ceil(mod(wfSamps,SAVE_CHK)/SAVE_CHK) );

% Import matlab library for saving to FITS with headers
import matlab.io.*

for i = 1:Nitr
    %- Setup indexing
    cbstr = (i-1)*SAVE_CHK+1;   %start index for this cube
    if (i == Nitr)  % on last iteration, bring in full cube (including pupils)
        cbendin = size(wfresFULL,3);    % Index used to pull rest of cube
        cbendsv = wfSamps;          % Index used to only save wfres (not pupils)
    else            % on every other iteration, just pull a normal cube piece
        cbendin = i*SAVE_CHK;       % Index used to pull current cube chunk
        cbendsv = cbendin;          % Index to save just the current cube chunk
    end

    %- Provide status update
    fprintf('Cube %3d of %3d  |  from %5d to %5d |  save til %5d\n', i, Nitr, cbstr, cbendin, cbendsv)
    
    %- Pull relevant cube chunk
    wfres = wfresFULL(:,:,cbstr:cbendin); 
    
    %- Upsample slightly to allow better centering
    upsp = 2;   % upsampling factor
    wfres = interpn(wfres,1:1/(2^upsp):size(wfres,1),(1:1/(2^upsp):size(wfres,2))',1:size(wfres,3),'nearest');
    % Recenter
    wfTMP = wfres(:,:,1);
    [rws,cls] = find(wfTMP);
    rw1 = min(rws);rw2 = max(rws);
    cl1 = min(cls); cl2 = max(cls);
    % Note: negative corrections below may require larger padding in WF Extraction
        % section above (padarray(wfres, [1 1 0]))
    cntfdg = [0 -2];   % centering fudge factor  [row col] ** must be even
    cnt = [rw1+cntfdg(1)/2, rw2+cntfdg(1)/2,   cl1+cntfdg(2)/2, cl2+cntfdg(2)/2];
    wfres = wfres(cnt(1):cnt(2),cnt(3):cnt(4),:);
    %- Pad matrix before up-sampling to prevent cropping from 'nearest' in interp2
    wfres = padarray(wfres, [1 1 0]);
    %- Upsample to PUPIL size   
    wfTMP = wfres(:,:,1);
    [rws,cls] = find(wfTMP);
    crc = 2*max(sqrt((rws-size(wfTMP,1)/2-1).^2 + (cls-size(wfTMP,2)/2-1).^2));
    sclfdg = 50;    % fudge factor for pupil resizing ** (-)=grow mask; (+)=shrink mask
    puprt = (pupCircDiam+sclfdg)/crc;   % Ratio of the two pupils
    rw = -(size(wfTMP,1)-1)*puprt/2:puprt:(size(wfTMP,1)-1)*puprt/2;
    cl = -(size(wfTMP,2)-1)*puprt/2:puprt:(size(wfTMP,2)-1)*puprt/2;
    rwi = rw(1):rw(end);
    cli = cl(1):cl(end);
    z = 1:size(wfres,3);
    zi = z;
    [rwt, clt, zt] = ndgrid(rw,cl,z); %create vector arrays
    [rwit, clit, zit] = ndgrid(rwi,cli,zi); %create vector arrays
    wfres = interpn(rwt,clt,zt,wfres,rwit,clit,zit, 'nearest');
    %- Pad to match Simulated pupil diameter
    wfres = pad_crop(wfres,size(PUPIL));
    
    
    %- Save cube chunk    
    %-- Save raw, reformatted, cropped WF
    % Create fits instance
    fitmap= fits.createFile([svfld sprintf('wfr_reform%06d.fits',i)]);
    fits.createImg(fitmap,'double',size(wfres));
    %header data
    fits.writeKey(fitmap, 'N       ', N);
    fits.writeKey(fitmap, 'apRad   ', apRad);
    fits.writeKey(fitmap, 'pupCrcDm', pupCircDiam);
    fits.writeKey(fitmap, 'lambda0 ', lambda0);
    fits.writeKey(fitmap, 'recNMds ', recNMds);
    fits.writeKey(fitmap, 'wfStrt  ', wfstrt);
    fits.writeKey(fitmap, 'wfSamps ', wfSamps);
    fits.writeKey(fitmap, 'wfStepSz', wfStepSz);
    fits.writeKey(fitmap, 'wfSCL   ', wfSCL);
    fits.writeKey(fitmap, 'PUP_rot ', PUP_rot);
    fits.writeKey(fitmap, 'upsp    ', upsp);
    fits.writeKey(fitmap, 'cntfdgRW', cntfdg(1));
    fits.writeKey(fitmap, 'cntfdgCL', cntfdg(2));
    fits.writeKey(fitmap, 'sclfdg  ', sclfdg);
    % save 
    fits.setCompressionType(fitmap,'NOCOMPRESS');
    fits.writeImg(fitmap,wfres);    
    % close instance
    fits.closeFile(fitmap);
    if isSaveFit && (cbendin-cbendsv)
        %-- On last iteration; Save the pupil masks
        pups = wfres(:,:,end-2+1:end);
        pups = cat(3,pups,PUPIL);
        fitswrite(pups, [svfld 'wfr_pupils.fits']);
    end
    
end

% Optional plot to see the two pupils superimposed
% wfTMP = wfres(:,:,1);
% wfTMP(wfTMP == 0) = nan;
% figure(); imagesc((wfTMP+0.2).*PUPIL);      % change 0.2 factor to improve visibility
% axis image;
% pupMatchTst = wfTMP.*PUPIL;
% % When tuning params, line below should yield 0 (this ensures none of the Keck
%   % PUPIL is beyond the PyWFS pupil
% sum(isnan(pupMatchTst(logical(PUPIL))))

%-- Reclaim un-needed variables to reduce memory usage
clear cl cli clit cls clt rw rwi rwit rws rwt z zi zit zt coords wfres wfresFULL pups

%% Main Loop (iterate through WFs and extract results)
%-- Preallocate realWF zernike coeffs FOR PLOTTING
    % NOTE: WF reconstruction is still done over 36 zernikes
coeffs = nan(wfSamps, length(nolls));
%-- Preallocate total RMS WFE matrix
rmsWF = nan(wfSamps, 2);
rmsWF(:, 1) = rmsWF_ch;
clear rmsWF_ch
%-- Preallocate other values
recCof = nan(recNMds,1);
recWF = zeros(size(PUPIL));

for i = 1:wfSamps
%% Load wfres chunk
if ~mod(i-1, SAVE_CHK)
    % read in new cube of WF data
    cbIND = floor(i/SAVE_CHK)+1 - ~mod(i,SAVE_CHK); % Index of current file
    wfres = fitsread([svfld sprintf('wfr_reform%06d.fits',cbIND)]);
end

%-- Extract frame from existing cube
frIND = mod(i-1,SAVE_CHK)+1;
wfTMP = wfres(:,:, frIND);
    
fprintf('wfSamp %06d of %06d (frame %3d in cube %3d)\n', i, wfSamps, frIND, cbIND);
    
%% Define pupil field

%-- DECOMPOSE WAVEFRONT INTO ZERNIKES;
recCof = recCof*nan;    % Clear the instance coeff vector
% Preallocate reconstructed WF matrix
recWF = recWF*0;        % Clear the instance reconstructed WF
for j = jvls
    % Extract current mode
    tmpzrn = unitzrn(:,:,j);
    % Left-divide to get coefficient
    recCof(j) = tmpzrn(logical(PUPIL))\wfTMP(logical(PUPIL));
    recWF = recWF + tmpzrn*recCof(j);
end
coeffs(i,:) = recCof(nolls);

if isPlotResCof
    %-- DISPLAY REQUESTED ZERNIKE COEFFICIENTS (nolls)
    if i == 1
        % Show requested coefficients
        zrnFig = figure(100);
        set(zrnFig, 'color', 'w', 'Visible', visflg);
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
        legH = legend(legs, 'Location', 'eastoutside');
        title(legH, 'Noll Index')
    else
        for j = 1:length(nolls)
            set(cofPlt(j), 'YData', coeffs(:,j)*zrnSCL);
        end
    end
end

%-- Renormalize one last time to ensure RMS values match.
rmsWF(i,2) = sqrt(mean(wfTMP(logical(PUPIL)).^2));
recWF = recWF*rmsWF(i,2)/sqrt(mean(recWF(logical(PUPIL)).^2));
%recWF = recWF*(recWF(logical(PUPIL))\wfTMP(logical(PUPIL)));

%-- DISPLAY PUPIL PROJECTIONS
if isPlotSimp && (i == 1)
    crcPup = makeCircularPupil(pupCircDiam/2, N);
    wfPLT = wfTMP; wfPLT(wfPLT==0) = nan;
    pupfig = figure(103); 
    set(pupfig, 'color', 'w', 'Visible', visflg);
    imagesc(xvals/pupCircDiam, yvals/pupCircDiam, wfPLT+0.05*PUPIL+0.1*crcPup);
    clear wfPLT crcPup
    axis image
    %axis([-1 1 -1 1])
    title('Pupil Projections')
end

% Convert WFE from waves to nm
wfTMP = recWF*lambda0*1e9.*PUPIL;      %.*PUPIL Optional to show with pupil
% Set out-of-pupil values to NaN for plotting
wfTMP(wfTMP==0) = nan;
% Store back into wfres to avoid creating 2 large matrices
wfres(:,:,frIND) = wfTMP;

%% Plot if desired
if isPlotSimp
    %-- Display WFR part of error only (top left plot in SPIE figure)
    if i == 1
        %-- First iteration: make the figure
        fig = figure(7);
        set(fig, 'units', 'inches', 'Position', [0 0 7 6], 'color', 'w', 'Visible', visflg);
        % Plot
        datWFR = imagesc(xvals/pupCircDiam,yvals/pupCircDiam,wfres(:,:,frIND));  
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
        set(datWFR, 'CData', wfres(:,:,frIND));           
    end

    drawnow;
end

%% Save Results
if isSaveFit && ~(mod(i,SAVE_CHK))
    %-- Save reconstructed WF cube with FITS_N frames 
    fitswrite(wfres, [svfld sprintf('wfr_recons%06d.fits',floor(i/SAVE_CHK))]);
    % OLD METHOD: when wfres was the full cube, not just a chunk
    %fitswrite(wfres(:,:,i-SAVE_CHK+1:i), [svfld sprintf('wfr_recons%06d.fits',floor(i/SAVE_CHK))]);
elseif isSaveFit && i==wfSamps
    %-- Save final cube in case it is not an integer multiple of FITS_N
    fitswrite(wfres, [svfld sprintf('wfr_recons%06d.fits',floor(i/SAVE_CHK)+1)]);
    % OLD METHOD: when wfres was the full cube, not just a chunk
    %fitswrite(wfres(:,:,i-mod(i,SAVE_CHK)+1:i), [svfld sprintf('wfr_recons%06d.fits',floor(i/SAVE_CHK)+1)]);
end

if isSaveRMS_COF && ~(mod(i,SAVE_CHK))
    %-- Save zernike coeffs and rmsWF with SAVE_CHK frames 
    fitswrite(rmsWF(i-SAVE_CHK+1:i,:), [svfld sprintf('rms_WF%06d.fits',floor(i/SAVE_CHK))]);
    fitswrite(coeffs(i-SAVE_CHK+1:i,:), [svfld sprintf('zern_coeff%06d.fits',floor(i/SAVE_CHK))]);
elseif isSaveRMS_COF && i==wfSamps
    %-- Save final fits in case it is not an integer multiple of SAVE_CHK
    fitswrite(rmsWF(i-mod(i,SAVE_CHK)+1:i,:), [svfld sprintf('rms_WF%06d.fits',floor(i/SAVE_CHK)+1)]);
    fitswrite(coeffs(i-mod(i,SAVE_CHK)+1:i,:), [svfld sprintf('zern_coeff%06d.fits',floor(i/SAVE_CHK)+1)]);
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

if isFigsVisEnd
    h = findobj('type', 'figure');
    set(h, 'Visible', 'on')
end