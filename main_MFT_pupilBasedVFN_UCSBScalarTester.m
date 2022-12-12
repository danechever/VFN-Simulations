%{
%}

clear; %close all; 
addpath(genpath(fullfile('VFNlib')));
addpath(genpath(fullfile('..','falco-matlab')))

%% Input parameters 

%-- Define Filenames
  % NOTE::: Notice the "%d" replacing what should be the input wavelengths
  % I assume all files will have the same format so this %d lets us
  % programatically load the files while also knowing the input wavelengths
  % Input matrices should be square, powers of 2 (256x256, 2048x2048, 4096x4096, etc.)
  % TLDR ==>> only provide one filename for each type of file and use %d
  %           to replace wavelength values and only provide square matrices
ampflnm = 'UCSB_SVFN_MaskDemo\\transmission_%2dum_SiSilica.txt';     % Amplitude (transmission) file
phsflnm = 'UCSB_SVFN_MaskDemo\\phase_%2dum_SiSilica.txt';     % Amplitude (transmission) file

%-- Define spatial sampling in pupil plane
  % (This should be your input grid size: 256x256, 1024x1024, 4096x4096, etc.)
  % NOTE::: powers of 2 make the code run significantly faster due to optimizations in the backend
N = 2^8;

%-- Define wavelength simulation info
% Provide list of avalable input file wavelengths
  % These should only include the wavelengths for which there are files
  % available to load
inlambdas = [2.2, 2.3, 2.4];     % [microns]
% Now define number of wavelenth samples to use for the sim:
  % if != to length of inlambdas, will linearly interpolate vortex data
  % NOTE::: strongly recommend ODD NUMBER and at least as many as inlambdas
  % Look at output print statements to make sure sampling makes sense...
numWavelengths = 3; % number of discrete wavelengths  to use

%% (OPTIONAL) Advanced features to mess with in simulation

%-- Option to use a standard/traditional vortex mask 
    % (In case you want to compare a test vortex against a theoretical one)
% Flag to choose if theoretical vortex should be used
isTheory = false;   % false: use provided vortex files, true: simulate a theoretically ideal vortex
charge = 1;         % if isTheory=true, this will be the charge of the ideal vortex simulated
isTheoryScalar = false; % true: use classical scalar model (lam/lam0) decay, false: assuming perfect achromatic

%-- Define wavefront error at the central wavelength
  % 1: Pist, 2: Tip, 3: Tilt, defoc, ob.astig, ver.astig, ver.coma, hor.coma,
  % 9: ver.tref, ob.tref, spher
nolls = 4:8;    % Noll indices to simulate
coeffs = [0.0, 0.0, 0.0, 0.0, 0.0];  % amplitude in [Waves RMS]

% Since using falco MFT propogator, can choose our final image sampling
lambda0Fnum_samp = 50; %[samples/ (lam0 F#)] in units of samples in the image plane
im_size = 10;    %[lam0/D] field of view in final image plane

%% Load the vortex specification files (If necessary)
if ~isTheory
    % Preallocate matrices to hold vortex data
    EPMP = nan(N, N, length(inlambdas));      % vortex phase matrix
    EPMA = nan(N, N, length(inlambdas));      % vortex amplitude matrix

    % Iterate through provided input wavelengths and load files for each one
    for ch = 1:length(inlambdas)
        lam = inlambdas(ch);    % Current wavelength 
        % Load transmission file as amplitude directly
        EPMA(:,:,ch) = readmatrix(sprintf(ampflnm,round(lam*10)));
        % Load phase file as real-valued element that will go in exponential
        EPMP(:,:,ch) = readmatrix(sprintf(phsflnm,round(lam*10)));
        
        % Correct amplitude to be at most 1 (some files showed >1 transmission which is unphysical)
        EPMA = min(EPMA, ones(size(EPMA)));
        
    end
end
%% Upsample vortex data to get all desired wavelengths (linear interp)
%-- Upsample if necessary, create lambdas vector either way
if numWavelengths == length(inlambdas)
	% No resampling needed since already have desired number of wavelengths
    lambdas = inlambdas;
    fprintf('Using input wavelengths only (ie. no upsampling performed)\n')
else
    % Create vector of wavelengths to use for sampling
    lambdas = linspace(inlambdas(1), inlambdas(end), numWavelengths);
    fprintf('-- Increaseing wavelength sampling --\n\tInput Wavelengths:')
    disp(inlambdas)
    fprintf('\tNew Wavelengths:')
    disp(lambdas)

    % Won't be able to, and don't need to, upsample when running with
    % theoretical/generated vortex
    if ~isTheory
        % Create grid of coord points for input sample
        [X,Y,Z] = meshgrid(1:N, 1:N, inlambdas);

        % Create grid of query points for resample
        [Xq,Yq,Zq] = meshgrid(1:N, 1:N, lambdas);

        % Interpolate to upsample
        EPMP = interp3(X,Y,Z, EPMP, Xq,Yq,Zq);
        EPMA = interp3(X,Y,Z, EPMA, Xq,Yq,Zq);
        fprintf('-- VORTEX DATA UPSAMPLED--\n')
    end
end 

%-- Compute central wavelength
central_band_index = ceil(numWavelengths/2);
% Convert lambdas vector from um to m as needed by code
lambdas = lambdas * 1e-6;
lambda0 = lambdas(central_band_index);

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
Epup = nan(N,N,numWavelengths);
for ch = 1:numWavelengths
    Epup(:,:,ch) = exp(1i*phz*lambda0/lambdas(ch)).*PUPIL;
    
    subplot(1,numWavelengths,ch);
    imagesc(xvalsPP,yvalsPP,angle(Epup(:,:,ch)));
    axis image; 
    axis([-1 1 -1 1]);
    title(['Phase at ',num2str(lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(parula(256));
   
end
drawnow;

%% Get PSF without vortex mask
[iPSF_BB, PSF] = getPSF_mft(Epup, lambdas, foc, coordsPP, coordsFP);

figure(3)
imagesc(xvalsFP,yvalsFP,iPSF_BB);
axis image; 
%axis([-3 3 -3 3]);
title('broadband PSF w/o vortex');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Make vortex mask
if isTheory
    % Generate a theoretically ideal vortex
    if isTheoryScalar
        charge = charge*lambda0./lambdas;   % classical (chromatic) scalar model
    else
        charge = charge*ones(1,numWavelengths); % achromatic model
    end
    EPMP = generateVortexMask( charge, coordsPP, [0 0] );
    EPMA = ones(size(EPMP));    % Assume perfect transmission
    fprintf('Generated a theoretically ideal vortex\n-- DID NOT USE INPUT VORTEX FILES --\n')
else
    % Exponentiate phase matrix so that it is truly phase
    EPMP = exp(1i*EPMP);
    fprintf('Used input vortex files\n')
end

% Combine vortex phase and amplitude into single mask
EPM = EPMA.*EPMP;


figure(4);
for ch = 1:numWavelengths
    % Show phases
    subplot(2,numWavelengths,ch);
    imagesc(xvalsPP,yvalsPP,angle(EPM(:,:,ch).*Epup(:,:,ch)));
    axis image; 
    axis([-1 1 -1 1]);
    title(['Phase at ',num2str(lambdas(ch)*1e9),'nm']);
    colormap(gca(), hsv(256));
    colorbar; 
   
    % Show amplitudes
    subplot(2,numWavelengths,ch+numWavelengths);
    imagesc(xvalsPP,yvalsPP,abs(EPM(:,:,ch).*Epup(:,:,ch)));
    axis image; 
    axis([-1 1 -1 1]);
    title(['Transmission at ',num2str(lambdas(ch)*1e9),'nm']);
    colormap(gca(), parula(256));
    caxis([0,1])
    colorbar; 
    
end
drawnow;

%% Get PSF with vortex mask
[iPSFv_BB, PSFv] = getPSF_mft(Epup.*EPM, lambdas, foc, coordsPP, coordsFP);


figure(5)
imagesc(xvalsFP,yvalsFP,iPSFv_BB);
axis image; 
axis([-3 3 -3 3]);
title('broadband PSF w/ vortex');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Generate fibermode at each lambda     
%-- Iterate through wavelengths generating modes
fibmode = nan(Nxi, Nxi, numWavelengths);
for ch = 1:numWavelengths
    % Generate the fiber mode for this wavelength with proper scaling
	fibmode(:,:,ch) = generateSMFmode_mft( fiber_props, lambdas(ch), coordsFP.dx, coordsFP);
    
end

%% Calculate Throughputs
%-- Get null depth (using overlap integral)
eta_onAx = nan(1,numWavelengths);
for ch = 1:numWavelengths
    eta_onAx(ch) = (abs(sum(sum(PSFv(:,:,ch).*fibmode(:,:,ch)))).^2)/totalPower0;
end

%-- Compute 2D coupling map
eta_maps = zeros(Nxi,Nxi,numWavelengths);
for ch = 1:numWavelengths
    % Compute the monochromatic coupling map (spatial units of lambda/D)
    eta_maps(:,:,ch) = generateCouplingMap( fibmode(:,:,ch), PSFv(:,:,ch), totalPower0, 5*lambda0Fnum_samp, coordsFP);

end

%-- Compute azimuthally-averaged performance
for ch = 1:numWavelengths
    % Compute average
    [tmpAvg, qvec] = radialAverage(eta_maps(:,:,ch), [Nxi/2+1 Nxi/2+1]);
    
    % Yes, I know it's not efficient to recreate matrix on every itr., I
    % don't care in this particular application...
    eta_pAvgs(:,ch) = tmpAvg;
end
% Since all eta_maps have the same spatial scale, they all have the same qvec
qvec = qvec/lambda0Fnum_samp; % [lam/D] (converted from pix to lam/D directly)

disp('Key Coupling Points:')
for ch = 1:numWavelengths
    % Find peak (value and separation at which it occurs) in radial average
    [peak, pind] = max(eta_pAvgs(:,ch),[],'all', 'linear');
    fprintf('lambda = %0.1f nm,    on-axis coup = %0.2e,    peak = %0.2f %% (at %0.2f lam0/D)\n',lambdas(ch)*1e9, eta_onAx(ch), peak*100, qvec(pind));
end

%% Display Plots of results
%-- Coupling Curves
figure(6);
for ch = 1:numWavelengths
    plot(qvec, eta_pAvgs(:,ch)*100, '--')
    if ch == 1
        hold on;
    end
end
xlabel('Separation [\lambda_0/D]')
ylabel('Coupling Efficiency [%]')
title('Azimuthally-Averaged Coupling')
xlim([0, 3]);
legend(string(lambdas*1e6)+' um')
hold off;

%-- Coup map linear scale
figure(7);
for ch = 1:numWavelengths
    subplot(1,numWavelengths,ch);
    imagesc(xvalsFP,yvalsFP,eta_maps(:,:,ch));
    axis image; 
    axis([-3 3 -3 3]);
    title(['\eta at ',num2str(lambdas(ch)*1e9),'nm']);
    xlabel('\lambda_0/D')
    ylabel('\lambda_0/D')
    colorbar; 
    colormap(gray(256));
end
drawnow

%-- Coup map log Scale
figure(8);
for ch = 1:numWavelengths
    subplot(1,numWavelengths,ch);
    imagesc(xvalsFP,yvalsFP,log10(eta_maps(:,:,ch)));
    axis image; 
    axis([-3 3 -3 3]);
    title(['log10(\eta) at ',num2str(lambdas(ch)*1e9),'nm']);
    xlabel('\lambda_0/D')
    ylabel('\lambda_0/D')
    colorbar; 
    colormap(gray(256));
end
drawnow

%% Display broadband coupling map if applicable
% This is what the photodiode should see when broadband light is applied

if numWavelengths > 1
    % Compute BB coupling map as average of all coupling maps
      % (This assumes the input spectrum is is flat --> ie. equal power at
      % all wavelengths)
    eta_map_BB = mean(eta_maps,3);
    
    % Compute average
        % NOTE::: This averaging isn't necessary since we can just take the
        % mean() of eta_pAvgs and it will give the same result but I do it
        % here for completeness
    [eta_pAvg_BB, ~] = radialAverage(eta_map_BB, [Nxi/2+1 Nxi/2+1]);
    
    % Plot Linear Scale
    figure(9);
    imagesc(xvalsFP, yvalsFP, eta_map_BB);
    axis image;
    axis([-3 3 -3 3]);
    title('Broadband \eta');
    xlabel('\lambda_0/D')
    ylabel('\lambda_0/D')
    colorbar;
    colormap(gray(256));
    
    % Plot Log Scale
    figure(10);
    imagesc(xvalsFP, yvalsFP, log10(eta_map_BB));
    axis image;
    axis([-3 3 -3 3]);
    title('Broadband log10(\eta)');
    xlabel('\lambda_0/D')
    ylabel('\lambda_0/D')
    colorbar;
    colormap(gray(256));
    
    % Plot broadband coupling efficiency curve
    figure(11);
    plot(qvec, eta_pAvg_BB*100, 'LineWidth', 2)
    xlabel('Separation [\lambda_0/D]')
    ylabel('Coupling Efficiency [%]')
    title('Azimuthally-Averaged Broadband Coupling')
    xlim([0, 3]);
    
    % Print key values
    disp('---')
    [peak_BB, pind_BB] = max(eta_pAvg_BB);
    fprintf('Broadband Performance:     on-axis coup = %0.2e,    peak = %0.2f %% (at %0.2f lam0/D)\n', mean(eta_onAx), peak_BB*100, qvec(pind_BB));
    
    fprintf('Broadband Integration-Time Reduction (eta_s/eta_p^2):  %0.1e\n', mean(eta_onAx)/(peak_BB^2));
else
    % When only a single wavelength is considered, set _BB to single point
    % so that the next section still works
    peak_BB = peak;
    pind_BB = pind;
end

%% Now that we know where peak is in broadband, plot null and peak vs. lambda at that location
%-- Null vs. Wavelength
figure(12);
semilogy(lambdas*1e6, eta_onAx, 'o-')
xlabel('Wavelength [um]')
ylabel('Null')
title('Null vs. Wavelength')
% Set reasonable min on ylim so we don't scale plotting to machine error...
tmpy = ylim;
ylim([max(tmpy(1),1e-8), 1])

%-- Peak vs. Wavelength
figure(13);
plot(lambdas*1e6, eta_pAvgs(pind_BB,:)*100, 'o-')
xlabel('Wavelength [um]')
ylabel('Coupling Efficiency [%]')
title([sprintf('Peak vs. Wavelength at %0.2f', qvec(pind_BB)), ' \lambda_0/D'])
ylim([0, 23])

