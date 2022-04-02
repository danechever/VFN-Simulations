clear; close all; 
%addpath('/home/vfndev/Documents/MATLAB/VFN-Simulations/VFNlib');
addpath(['..' filesep 'VFNlib']);
%% Input parameters 

% Define path where output files will be saved
svPth = '/media/Data_Drive/VFN/ScalarOnPrimaryData/cmap_optimizations';

% Define Sampling Info - initialized as struct for easy access
vars.N = 2^11;
vars.apRad = 128;
vars.apDia0 = 2 * vars.apRad;

%Define wavelength info
vars.lambda0 = 2.2e-6;
vars.fracBW = 0.1818;
vars.numWavelengths = 5;
vars.lambdas = getWavelengthVec(vars.lambda0, vars.fracBW, vars.numWavelengths);% array of wavelengths (meters)

%Define charge of the vortex mask at each wavelength
%charge = ones(1,numWavelengths); % achromatic
vars.charge = vars.lambda0./vars.lambdas;

%Define wavefront error at the central wavelength*
vars.nolls = 3;
vars.coeffs = 0*0.1;

%Give offsets for the vortex mask
vars.offsetX = 0;
vars.offsetY = 0;

vars.numRings = 3;
vars.wGap = 25.4/10916*vars.apDia0/2;

vars.relTints = zeros(1);

%% Fiber Properties
% Parameters for Thorlabs SM600
%Using SM2000 in this version of the code
%link: https://www.thorlabs.com/NewGroupPage9_PF.cfm?ObjectGroup_ID=949

fiber_props.core_rad = 5.5e-6; % Core radius [um]
fiber_props.n_core = 1.4436; %1.4571; core index (interpolated from linear fit to 3 points)
fiber_props.n_clad = 1.4381; %1.4558; cladding index (interpolated from linear fit to 3 points)
fiber_props.type = 'bessel';
fiber_props.Fnum = getMFD(fiber_props,vars.lambda0)/(vars.lambda0*1.4); % focal ratio of the beam at the fiber

%% Generate the coordinate system

vars.coords = generateCoordinates(vars.N); % Creates NxN arrays with coordinates 
xvals = vars.coords.xvals; % Helpful for plotting
yvals = vars.coords.yvals;

%% Create array with pupil function

vars.PUPIL = makeKeckPupil(2*vars.apRad, vars.N );
lambdaOverD = vars.N/vars.apRad/2; % lam/D in units of samples in the image plane
[vars.normI, vars.totalPower0] = getNormalization(vars.PUPIL);% Normalization factors
vars.lambdaOverD = vars.N/vars.apRad/2; % lam/D in units of samples in the image plane

 figure(1);
 imagesc(xvals/vars.apRad,yvals/vars.apRad,vars.PUPIL);
 axis image; 
 axis([-1 1 -1 1]);
 title('Pupil');
 colorbar; 
 colormap(parula(256));
 %grid on;
 drawnow;
 
 %% Generate Fiber Modes
 fibmodes = NaN(vars.N, vars.N, vars.numWavelengths);
 for ch = 1:vars.numWavelengths
     fibmodes(:,:,ch) = generateSMFmode(fiber_props, vars.lambdas(ch), fiber_props.Fnum, vars.lambdaOverD, vars.coords);
 end
 
 %% Create matrix with tilt and pistons for the optimization seed
 %xtilt, ytilt, piston
 
 initial = zeros(36,3);

 loc = 1;
 for ringNum = 1:3 %adds absolute value of tilt to each segment in the position they are added.
     for seg = (loc):(loc+ringNum*6-1)
         initial(seg,1:2) = (1*2*pi)/(ringNum*6);
     end
     loc = (loc+ringNum*6);
 end
 
hexFlatDiam = (vars.apDia0-3*2*vars.wGap)/(2*3+1);
vars.hexSep = hexFlatDiam + vars.wGap;

vars.hexAmpConst = NaN(vars.N, vars.N, 36);
vars.hexPhzConst = NaN(vars.N, vars.N, 36);

segs = ones(1,36);
count = 1;
 for ringNum = 1:3
     crow = ringNum * vars.hexSep;
     ccol = 0;
     t = atan2(crow,ccol);
     t = round(t,3);
     
     if(segs(count) == 1)
         initial(count,1) = initial(count,1) * -sin(t);
         initial(count,2) = initial(count,2) * cos(t);
         initial(count,3) = t/(2*pi);
         
         [vars.hexAmpConst(:,:,count), vars.hexPhzConst(:,:,count)] = generateHexConstants(crow, ccol, vars.numRings, vars.apDia0, vars.wGap, zeros(vars.N));
         
     end

     count = count + 1;
     for face = 1:6
         step_dir = pi/6*(2*face+5);
         steprow = vars.hexSep*sin(step_dir);
         stepcol = vars.hexSep*cos(step_dir);
         stepnum = 1;
         
         while(stepnum <= ringNum && count <= 36)
             crow = crow + steprow;
             ccol = ccol + stepcol;
             
             t = atan2(crow,ccol);
             t = round(t,3);
             
             if(face==6 && stepnum ==ringNum)
%                  disp('finished ring');
             else
                 if(segs(count) == 1)
                    initial(count,1) = initial(count,1) * -sin(t);
                    initial(count,2) = initial(count,2) * cos(t);
                    initial(count,3) = t/(2*pi);
                    
                    [vars.hexAmpConst(:,:,count), vars.hexPhzConst(:,:,count)] = generateHexConstants(crow, ccol, vars.numRings, vars.apDia0, vars.wGap, zeros(vars.N));

                 end
             count = count + 1;
             end
             stepnum = stepnum + 1;
         end
     end
 end
 
%  initial(:,1) = 0;
%  initial(:,2) = 0;
 
%addpath('/home/vfndev/Documents/MATLAB/falco-matlab/lib/utils');
addpath(['..' filesep '..' filesep 'falco-matlab' filesep 'lib' filesep 'utils']);
 
 %-- Decrease matrix size in pupil plane to reduce runtime
 vars.PUP_CRP_SZ = round(2.1*vars.apRad);
 vars.hexAmpConst = pad_crop(vars.hexAmpConst,vars.PUP_CRP_SZ);
 vars.hexPhzConst = pad_crop(vars.hexPhzConst,vars.PUP_CRP_SZ);
 
%  disp(count); 
%  disp(initial);

%% Define the initial phase for the optimization seed
phz_initial = angle(makeKeckPupilInputs( vars, initial ));

figure(3); %Phase of the optimization seed
for ch = 1:vars.numWavelengths
    Epup_initial(:,:,ch) = exp(1i*phz_initial*vars.lambda0/vars.lambdas(ch)).*vars.PUPIL;
    
    subplot(1,vars.numWavelengths,ch);
    imagesc(xvals/vars.apRad,yvals/vars.apRad,angle(Epup_initial(:,:,ch)));
    axis image; 
    axis([-1 1 -1 1]);
    title(['Phase at ',num2str(vars.lambdas(ch)*1e9),'nm. Optimization Seed.']);
    colorbar; 
    colormap(hsv(256));
end
drawnow;

%% Relative integration time function handle
nItrCap = 1000;
vars.relTintArray = nan(nItrCap*1.2,1);
vars.n_calls = 1;
ssMin = @(inputs) makeKeckPupilLeastSquaresRelTintOpt( vars, inputs, fiber_props, fibmodes );

%% Optimization

options = optimoptions('fmincon','MaxFunctionEvaluations', nItrCap);
[localMinX,fvalY] = fmincon(ssMin, initial,[],[],[],[],[],[],[], options);
disp(localMinX);
disp(fvalY);

%disp(vars.relTintArray);

%% Display Optimum Phase Pattern

optimum = localMinX;
save([svPth filesep 'tt_opt.mat'], 'optimum');
optPhz = angle(makeKeckPupilInputs( vars, optimum ));

figure(5); %Displays the optimized phase on the Keck pupil at the resolution it was calculated at
for ch = 1:vars.numWavelengths
    Epup_opt(:,:,ch) = exp(1i*optPhz*vars.lambda0/vars.lambdas(ch)).*vars.PUPIL;
    
    subplot(1,vars.numWavelengths,ch);
    imagesc(xvals/vars.apRad,yvals/vars.apRad,angle(Epup_opt(:,:,ch)));
    axis image;
    axis([-1 1 -1 1]);
    title(['Phase at ',num2str(vars.lambdas(ch)*1e9),'nm. Optimized.']);
    colorbar; 
    colormap(hsv(256));
end

%% Generates Coupling Map for Optimised Solution

iPSF_BB = getPSF(Epup_opt,vars.lambda0,vars.lambdas,vars.normI,vars.coords);

figure(9)
imagesc(xvals/vars.lambdaOverD,yvals/vars.lambdaOverD,iPSF_BB);
axis image;
axis xy;
axis([-2 2 -2 2]);
title('broadband PSF');
colorbar;
colormap(parula(256));
drawnow;

% Generate coupling maps for each wavelength

% Parameters for Thorlabs SM600
%Using SM2000 in this version of the code
    % link: https://www.thorlabs.com/NewGroupPage9_PF.cfm?ObjectGroup_ID=949
fiber_props.core_rad = 5.5e-6;% Core radius [um]
fiber_props.n_core = 1.4436; %1.4571; core index (interpolated from linear fit to 3 points)
fiber_props.n_clad = 1.4381; %1.4558; cladding index (interpolated from linear fit to 3 points)
fiber_props.type = 'bessel';
Fnum = getMFD(fiber_props,vars.lambda0)/(vars.lambda0*1.4); % focal ratio of the beam at the fiber

eta_maps = generateCouplingMap_polychromatic( Epup_opt, fiber_props, vars.lambda0, Fnum, vars.lambdas, vars.totalPower0, vars.lambdaOverD, 3*vars.lambdaOverD, vars.coords);

%-- Find Relative Integration Time
eta_map_mean = mean(eta_maps,3);
crp = 2*0.5*vars.lambdaOverD;
c_map = eta_map_mean(end/2+1-floor(crp/2):end/2+1+floor(crp/2),end/2+1-floor(crp/2):end/2+1+floor(crp/2)); %the centroid
eta_s = min(c_map,[],'all'); %minimum value in the centroid
[min_ind(1),min_ind(2)] = find(eta_s==c_map); %indices of minimum value in the centroid
min_ind = min_ind + (vars.N/2-floor(crp/2)); %adjust min values to reflect position in map, not cmap
        
eta_sX = min_ind(2);
eta_sY = min_ind(1);

[radProf2, rvec] = VFN_An_radAverage(eta_map_mean,[eta_sX, eta_sY]);
[radMx, ind] = max(radProf2);
rad = rvec(ind);

%-- Calculate relative integration time
relTints = eta_s./(radMx.^2);

figure(10);
for ch = 1:vars.numWavelengths
    subplot(1,vars.numWavelengths,ch);
    imagesc(xvals/vars.lambdaOverD, yvals/vars.lambdaOverD, log10(eta_maps(:,:,ch)));
    hold on;
    viscircles((min_ind/vars.lambdaOverD - (vars.N/2+1)/vars.lambdaOverD), rad/vars.lambdaOverD);
    plot((min_ind(1)/vars.lambdaOverD - (vars.N/2+1)/vars.lambdaOverD), (min_ind(2)/vars.lambdaOverD - (vars.N/2+1)/vars.lambdaOverD), 'r+', 'LineWidth', 1);
    hold off;
    axis image; 
    axis([-2 2 -2 2]);
    caxis([-3 -0.5])
    title(['\eta at ', num2str(vars.lambdas(ch)*1e9), 'nm']);
    colorbar;
    colormap(gray(256));
end

figure(11);
imagesc(xvals/vars.lambdaOverD, yvals/vars.lambdaOverD, log10(eta_map_mean));
hold on;
viscircles((min_ind/vars.lambdaOverD - (vars.N/2+1)/vars.lambdaOverD), rad/vars.lambdaOverD);
plot((min_ind(1)/vars.lambdaOverD - (vars.N/2+1)/vars.lambdaOverD), (min_ind(2)/vars.lambdaOverD - (vars.N/2+1)/vars.lambdaOverD), 'r+', 'LineWidth', 1);
hold off;
axis image; 
axis([-2 2 -2 2]);
caxis([-3 -0.5])
title('\eta map averaged');
colorbar;
colormap(gray(256));







