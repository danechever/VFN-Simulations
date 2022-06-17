function [noisy_image, PSF] = CRED2ImSim_getNoisyIm(cred2sim, Tint, nolls, coeffs, tilts, verbose)
%[noisy_image, PSF] = setUp_CRED2ImSim(cred2sim, tint, nolls, coeffs, tilts)
%   Setup the KPIC Phase II CRED2 Image simulator
%
%   Inputs:
%       cred2sim: struct with elements needed for the simulator
%       tint: [s] integration time for image 
%       nolls: noll indices of zernikes to add (corresponding to coeffs)
%       ---- Each of the following can be multi-D where the each row holds
%               the values to be used for each frame)
%       coeffs: [lambda0 waves rms] zernike coefficient for each noll
%       tilts: (x, y) [mas] X and Y tilts to apply 
%   Outputs:
%       noisy_image: Final, CRED2 image with noise included
%       PSF: complex-valued electric field at CRED2 image plane (no noise)


%% Figue out how many times to iterate (ie. how many frames were requested)
if size(nolls,1) ~= 1
    error('nolls should be 1 row only; use coeffs rows to vary WF or samples')
end

if (size(nolls,2) ~= size(coeffs,2))
    error('nolls and coeffs must be of equal size along second dimension')
end

if (size(coeffs,1) ~= size(tilts,1))
    % User either provided a single value for one of these with the hopes
    % that we repeat it or they messed up
    if size(coeffs,1) == 1
        % User probably wants us to repeat the zernikes on all frames
        coeffs = repmat(coeffs, [size(tilts,1),1]);
    elseif size(tilts,1) == 1
        % User probably wants us to repeat the tilts on all frames
        tilts = repmat(tilts, [size(nolls,1),1]);
    else
        error('If unequal shape nolls/coeffs and tilts are provided, one must be singleton in first dim')
    end
end

NFrames = size(tilts,1);

%% Iterate through making frames
% Preallocate final arrays
noisy_image = nan(NFrames, cred2sim.coordsFP.N, cred2sim.coordsFP.N);
PSF = complex(noisy_image);

% Iterate
for frame = 1:NFrames
    if verbose && (mod(frame, floor(NFrames/10))==0) || (frame == 1)
        fprintf('Making img %4d of %4d\n', frame, NFrames);
    end
    %% Define pupil field (WFE)
    
    % TODO:::: IMPLEMENT THE ABILITY TO PROVIDE A STATIC SPECKLE FIELD
        % This could obviously be done by the user by repeating some
        % zernikes with the same amplitude in all frames but that would
        % require unnecessary computation to re-create the static field.
        %
        % Alternatively, could have use create the static phz and provide
        % it as a field in cred2sim. Then if the field is present, it's
        % added otherwise if it's not, it's set to 0. This static field
        % could also be made for the user in setUp_CRED2ImSim
    
    % Generate varying WFE based on provided zernikes
    phz = generateZernike_fromList( nolls, coeffs(frame,:), cred2sim.PUPIL, cred2sim.apRad, cred2sim.coordsPP); 

%     figure(2);
    Epup = nan(cred2sim.coordsPP.N,cred2sim.coordsPP.N,cred2sim.numWavelengths);
    for ch = 1:cred2sim.numWavelengths
        Epup(:,:,ch) = exp(1i*phz*cred2sim.lambda0/cred2sim.lambdas(ch)).*cred2sim.PUPIL;

%         subplot(1,cred2sim.numWavelengths,ch);
%         imagesc(cred2sim.xvalsPP,cred2sim.yvalsPP,angle(Epup(:,:,ch)));
%         axis image; 
%         axis([-1 1 -1 1]);
%         title(['Phase at ',num2str(cred2sim.lambdas(ch)*1e9),'nm']);
%         colorbar; 
%         colormap(parula(256));

    end
%    drawnow;

    %% Add tilts
    %%%%%% TODO:::::::::::::: Check if this should be /N or /(2*apRad)
    % Get tilt phase at central wavelength
    phz = 2*pi*tilts(frame,1)/cred2sim.lambda0OverKeckD_mas*cred2sim.coordsPP.X/(2*cred2sim.apRad);%cred2sim.coordsPP.N;
    phz = phz + 2*pi*tilts(frame,2)/cred2sim.lambda0OverKeckD_mas*cred2sim.coordsPP.Y/(2*cred2sim.apRad);%cred2sim.coordsPP.N;

%     figure(10)
    for ch = 1:cred2sim.numWavelengths
        Epup(:,:,ch) = Epup(:,:,ch).*exp(1i*phz*cred2sim.lambda0/cred2sim.lambdas(ch));

%         subplot(1,cred2sim.numWavelengths,ch);
%         imagesc(cred2sim.xvalsPP,cred2sim.yvalsPP,angle(Epup(:,:,ch)));
%         axis image; 
%         axis([-1 1 -1 1]);
%         title(['Phase With Tilt at ',num2str(cred2sim.lambdas(ch)*1e9),'nm']);
%         colorbar; 
%         colormap(parula(256));
    end
%    drawnow;


    %% Get PSF with vortex mask
    [~, PSFv] = getPSF_mft(Epup.*cred2sim.EPM, cred2sim.lambdas, cred2sim.foc, cred2sim.coordsPP, cred2sim.coordsFP);

%     % Compute broadband image intensity (without peak-normalization)
%     iPSFv_BB = mean(abs(PSFv).^2,3);
% 
%     figure(5)
%     imagesc(cred2sim.xvalsFP,cred2sim.yvalsFP,iPSFv_BB);
%     axis image; 
%     %axis([-3 3 -3 3]);
%     title('broadband PSF w/ vortex');
%     cb = colorbar;
%     cb.Label.String = 'Intensity [ph/s]';
%     %caxis([-3 0])
%     colormap(parula(256));
%     drawnow;

    %% simulate the image with noise
    % 1. get photons per pixel in the pupil plane
        % photons/second/m2/time_sample
    % 2. multiply by system throughput/transmission
    % 3. take sqrt() and treat that as the magnitude of the electric field

    %-- 0) Square the PSF (E-field) to get intensity 
    iPSFv = abs(PSFv).^2;   % [photons/s/pixel]
    
    %-- 4) Account for integration time 
    iPSFv = Tint * iPSFv; % [ph/pix/t_sample]
    
    %-- 1) Given intensity in focal plane, add photon noise (in photons/s/pix)
    % draw from poisson distribution to get shot noise (on per-pixel basis)
    noisy_im = poissrnd(iPSFv);     % [photons/t_sample/pixel]

        % TODO::: CONFIRM THAT THIS IS INDEED THE NOISY IMAGE AND NOT TGHE
        % NOISE ITSELF

    %-- 2) Convert from photons/s/pix to electrons/s/pix using QE
    noisy_im = cred2sim.CRED2_QE * noisy_im;     % [e-/t_sample/pixel]

    %-- 3) Add Dark current (in electrons/t/pix)
    % Sample poisson distribution to get dark current per pixel
    dark_noise = poissrnd(cred2sim.CRED2_DarkCurrent*Tint,size(noisy_im)); % [e-/t_sample/pix]
    % Add noise
    noisy_im = dark_noise + noisy_im;   % [e-/t_sample/pix]

    %-- 5) Add Read noise (in e-/pix)
    % Sample gaussian distribution to get read noise per pixel
       % First 0 is to set noise to be mean-0
    read_noise = normrnd(0, cred2sim.CRED2_ReadNoiseRMS, size(noisy_im));
    noisy_im = read_noise + noisy_im;   % [e-/t_sample/pix]

    %-- 6) Convert to DN
    noisy_im = noisy_im / cred2sim.CRED2_ADU;   % [count/pix/t_sample]



%     % Display the noisy image!
%     figure();
%     imagesc(noisy_im);
%     axis image;
%     title('Noisy Image');
%     cb = colorbar;
%     cb.Label.String = 'Intensity [DN]';
%     colormap(parula(256));
%     drawnow;

    
    % Store results in output arrays
    noisy_image(frame, :,:) = noisy_im;
    PSF(frame, :,:) = PSFv;
end
end
