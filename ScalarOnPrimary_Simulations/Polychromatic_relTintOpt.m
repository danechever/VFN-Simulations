clear; close all; 
addpath(['..' filesep 'VFNlib']);

%% Input parameters 

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
                 disp('finished ring');
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
 
 initial(:,1) = 0;
 initial(:,2) = 0;
 
 addpath(['..' filesep '..' filesep 'falco-matlab' filesep 'lib' filesep 'utils']);
 %-- Decrease matrix size in pupil plane to reduce runtime
 vars.PUP_CRP_SZ = round(2.1*vars.apRad);
 vars.hexAmpConst = pad_crop(vars.hexAmpConst,vars.PUP_CRP_SZ);
 vars.hexPhzConst = pad_crop(vars.hexPhzConst,vars.PUP_CRP_SZ);
 
 disp(count); 
 disp(initial);

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

%% relative integration time Min Function
nItrCap = 1000;
vars.relTintArray = nan(nItrCap*1.2,1);
vars.n_calls = 1;
ssMin = @(inputs) makeKeckPupilLeastSquaresRelTintOpt( vars, inputs, fiber_props, fibmodes );

%% Optimization

% global relTintArray;
reltintArray = zeros(1);

% history = struct('x',[],'fval',[]);
% history(1) = [];
% save('history.mat');
options = optimoptions('fmincon','MaxFunctionEvaluations', nItrCap);
[localMinX,fvalY,history] = fmincon(ssMin, initial,[],[],[],[],[],[],[], options);
disp(localMinX);
disp(fvalY);

disp(relTintArray);
% load('history.mat');
% disp([history.x]);

%% Display Optimum Phase Pattern

optimum = localMinX;
save('optVarsRelTints2.mat', 'optimum');
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

%% Initialize new Pupil with greater resolution

% hres.N = 2^11; % Size of computational grid (NxN samples) 
% hres.apRad = 128; % Aperture radius in samples
% hres.apDia0 = 2 * hres.apRad;
% 
% %Define wavelength info
% hres.lambda0 = 2.2e-6;
% hres.fracBW = 0.1818;
% hres.numWavelengths = 5;
% hres.lambdas = getWavelengthVec(hres.lambda0,hres.fracBW,hres.numWavelengths);% array of wavelengths (meters)
% 
% %Define charge of the vortex mask at each wavelength
% %charge = ones(1,numWavelengths); % achromatic
% hres.charge = hres.lambda0./hres.lambdas;
% 
% %Define wavefront error at the central wavelength
% hres.nolls = 4:8;
% hres.coeffs = 0*[0.01,-0.0099,-0.0095,-0.0008,0.0033];
% 
% %Give offsets for the vortex mask
% hres.offsetX = 0;
% hres.offsetY = 0;
% 
% hres.numRings = 3;
% hres.wGap = 25.4/10916*hres.apDia0/2;
% 
% %Coordinate system
% coords = generateCoordinates(hres.N);% Creates NxN arrays with coordinates 
% hres_xvals = coords.xvals;% Helpful for plotting
% hres_yvals = coords.yvals;
% 
% hres_PUPIL = makeKeckPupil(2*hres.apRad, hres.N );
% [normI, totalPower0] = getNormalization(hres_PUPIL);% Normalization factors
% hres.lambdaOverD = hres.N/hres.apRad/2; % lam/D in units of samples in the image plane
% 
% figure(6);
% imagesc(hres_xvals/hres.apRad, hres_yvals/hres.apRad, hres_PUPIL);
% axis image; 
% axis([-1 1 -1 1]);
% title('Pupil');
% colorbar; 
% colormap(parula(256));
% grid on;
% drawnow;

%% Display optimum phase map at higher resolution
% optPhzHighRes = angle(makeKeckPupilInputs( hres, optimum ));
% 
% figure(8);
% for ch = 1:hres.numWavelengths
%     Epup4(:,:,ch) = exp(1i*optPhzHighRes*hres.lambda0/hres.lambdas(ch)).*hres_PUPIL;
%     
%     subplot(1,hres.numWavelengths,ch);
%     imagesc(hres_xvals/hres.apRad,hres_yvals/hres.apRad,angle(Epup4(:,:,ch)));
%     axis image;
%     axis([-1 1 -1 1]);
%     title(['Phase at ',num2str(hres.lambdas(ch)*1e9),'nm. Optimized.']);
%     colorbar; 
%     colormap(hsv(256));
% end

%% Generates Coupling Map for Optimised Solution

iPSF_BB = getPSF(Epup_opt,vars.lambda0,vars.lambdas,vars.normI,vars.coords);

figure(9)
imagesc(xvals/vars.lambdaOverD,yvals/vars.lambdaOverD,iPSF_BB);
axis image;
axis xy;
axis([-2 2 -2 2]);
title('broadband PSF');
colorbar;%caxis([-3 0])
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

figure(10);
for ch = 1:vars.numWavelengths
    subplot(1,vars.numWavelengths,ch);
    imagesc(xvals/vars.lambdaOverD, yvals/vars.lambdaOverD, log10(eta_maps(:,:,ch)));
    axis image; 
    axis([-2 2 -2 2]);
    caxis([-3 -0.5])
    title(['\eta at ', num2str(vars.lambdas(ch)*1e9), 'nm']);
    colorbar;
    colormap(gray(256));
end

%% Show eta_s and eta_p

%relTints = makeKeckPupilLeastSquaresETAOPT( vars, optimum, fiber_props, fibmodes);







