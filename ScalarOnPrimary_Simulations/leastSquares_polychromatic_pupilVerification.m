clear; close all;
addpath(['..' filesep 'VFNlib']);

%% Input parameters 

% Define initial sampling info (smaller so minimizer takes less time)
vars.N = 2^11;
vars.apRad = 128;
vars.apDia0 = 2 * vars.apRad;

%Define wavelength info
vars.lambda0 = 2.2e-6;
vars.fracBW = 0.1818;
vars.numWavelengths = 1;
vars.lambdas = vars.lambda0;

%Define charge of the vortex mask at each wavelength
%charge = ones(1,numWavelengths); % achromatic
vars.charge = vars.lambda0./vars.lambdas;

%Define wavefront error at the central wavelength
vars.nolls = 3;
vars.coeffs = 0*0.1;

%Give offsets for the vortex mask
vars.offsetX = 0;
vars.offsetY = 0;

vars.numRings = 3;
vars.wGap = 25.4/10916*vars.apDia0/2;

%% Generate the coordinate system

coords = generateCoordinates(vars.N);% Creates NxN arrays with coordinates 
xvals = coords.xvals;% Helpful for plotting
yvals = coords.yvals;

%% Create array with pupil function

PUPIL = makeKeckPupil(2*vars.apRad, vars.N );
%[normI, totalPower0] = getNormalization(PUPIL);% Normalization factors
lambdaOverD = vars.N/vars.apRad/2; % lam/D in units of samples in the image plane

 figure(1);
 imagesc(xvals/vars.apRad,yvals/vars.apRad,PUPIL);
 axis image; 
 axis([-1 1 -1 1]);
 title('Pupil');
 colorbar; 
 colormap(parula(256));
 drawnow;
 
 %% Create matrix with tilt and pistons
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
     disp(count);
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
 disp(count);
 
 disp(initial);
 
 addpath(['..' filesep '..' filesep 'falco-matlab' filesep 'lib' filesep 'utils']);
 vars.PUP_CRP_SZ = round(2.1*vars.apRad);
 vars.hexAmpConst = pad_crop(vars.hexAmpConst,vars.PUP_CRP_SZ);
 vars.hexPhzConst = pad_crop(vars.hexPhzConst,vars.PUP_CRP_SZ);
 
%  ps = [0.2500,    0.4167,   -0.4167,   -0.2500,   -0.0834,    0.0834,    0.2500,    0.3333,    0.4167,    0.5001,   -0.4167,   -0.3333,   -0.2500,   -0.1666,   -0.0834, 0,    0.0834,    0.1666,    0.2500,    0.3030,    0.3637,    0.4167,    0.4697,   -0.4697,   -0.4167,   -0.3637,   -0.3030,   -0.2500,   -0.1969,   -0.1364, -0.0834,   -0.0302,    0.0302,    0.0834,    0.1364,    0.1969];
%  
% initial(:,3) = ps;
%  initial(35:end,:) = 2*pi*rand(3,3);

%% Define ideal pupil field
%Make more figures for as many algorithms are being tested

phz_wfe = generateZernike_fromList( vars.nolls, vars.coeffs, PUPIL, vars.apRad, coords); 

figure(2); %Field of the ideal pupil
for ch = 1:vars.numWavelengths
    Epup_wfe(:,:,ch) = exp(1i*phz_wfe*vars.lambda0/vars.lambdas(ch)).*PUPIL;
    
    subplot(1,vars.numWavelengths,ch);
    imagesc(xvals/vars.apRad,yvals/vars.apRad,angle(Epup_wfe(:,:,ch)));
    axis image; 
    axis([-1 1 -1 1]);
    title(['Phase on ideal vortex at ',num2str(vars.lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(hsv(256));
   
end
drawnow;


phz_initial = angle(makeKeckPupilInputs( vars, initial ));

figure(3); %Phase of the optimization seed
for ch = 1:vars.numWavelengths
    Epup_initial(:,:,ch) = exp(1i*phz_initial*vars.lambda0/vars.lambdas(ch)).*PUPIL;
    
    subplot(1,vars.numWavelengths,ch);
    imagesc(xvals/vars.apRad,yvals/vars.apRad,angle(Epup_initial(:,:,ch)));
    axis image; 
    axis([-1 1 -1 1]);
    title(['Phase at ',num2str(vars.lambdas(ch)*1e9),'nm. Optimization Seed.']);
    colorbar; 
    colormap(hsv(256));
end
drawnow;

%% Make ideal vortex mask and Display Phase Map

EPM = generateVortexMask( vars.charge, coords, [vars.offsetX vars.offsetY] );
Epup_ideal = EPM.*Epup_wfe;

phz = EPM.*abs(Epup_wfe);

central_band_index = ceil(vars.numWavelengths/2);

figure(4) %The ideal phase superimposed on the Keck Pupil
imagesc(xvals/vars.apRad,yvals/vars.apRad,angle(Epup_ideal(:,:,central_band_index)));
axis image; 
axis([-1 1 -1 1]);
title('Vortex Superimposed on Pupil');
colormap(hsv(256));
colorbar; 
drawnow;

%% Compute Initial Variable Residual and Display on Phase map
res2 = (angle(Epup_initial) - angle(Epup_ideal));
sqsum2 = sum(res2.^2,'all');
figure(3);
txt2 = ['SQSUM Seed - Ideal: ' num2str(sqsum2)];
text(-0.5,0,txt2);


%% Phase Min Function

ssMin = @(inputs) makeKeckPupilLeastSquares( vars, inputs, Epup_ideal );

%% Optimization

% options = optimoptions('fmincon','Display','iter');
[localMinX,fvalY] = fmincon(ssMin, initial,[],[],[],[],[],[],[]);
disp(localMinX);
disp(fvalY);

%% Display Optimum Phase Pattern

optimum = localMinX;
save('optAnal.mat', 'optimum');
optPhz = angle(makeKeckPupilInputs( vars, optimum ));

figure(5); %Displays the optimized phase on the Keck pupil at the resolution it was calculated at
for ch = 1:vars.numWavelengths
    Epup_opt(:,:,ch) = exp(1i*optPhz*vars.lambda0/vars.lambdas(ch)).*PUPIL;
    
    subplot(1,vars.numWavelengths,ch);
    imagesc(xvals/vars.apRad,yvals/vars.apRad,angle(Epup_opt(:,:,ch)));
    axis image;
    axis([-1 1 -1 1]);
    title(['Phase at ',num2str(vars.lambdas(ch)*1e9),'nm. Optimized.']);
    colorbar; 
    colormap(hsv(256));
end

%% Compute Optimised Residual and Display on Phase Map
res1 = (angle(Epup_opt) - angle(Epup_ideal));
sqsum1 = sum(res1.^2,'all');
figure(5);
txt1 = ['SQSUM, Optimum - Ideal: ' num2str(sqsum1)];
text(-0.5,0,txt1);

%% Compute squared sum difference between optimum and seed
res3 = (angle(Epup_opt) - angle(Epup_initial));
sqsum3 = sum(res3.^2,'all');
figure(3);
txt3 = ['SQSUM, Optimum - Seed : ' num2str(sqsum3)];
text(-0.5,-0.2,txt3);

%% Generates Coupling Map for Optimised Solution
[normI, totalPower0] = getNormalization(PUPIL);% Normalization factors
vars.lambdaOverD = vars.N/vars.apRad/2; % lam/D in units of samples in the image plane
iPSF_BB = getPSF(Epup4,vars.lambda0,vars.lambdas,normI,coords);

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

eta_maps = generateCouplingMap_polychromatic( Epup4, fiber_props, vars.lambda0, Fnum, vars.lambdas, totalPower0, vars.lambdaOverD, 3*vars.lambdaOverD, coords);

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

% %% Input parameters 
% 
% % Define smapling info
% N = 2^11; % Size of computational grid (NxN samples) 
% apRad = 128; % Aperture radius in samples 
% 
% % Define wavelength info
% lambda0 = 2.2e-6; %central wavelength
% fracBW = 0.1818; %\Delta\lambda/\lambda
% numWavelengths = 1;% number of discrete wavelengths 
% lambdas = lambda0; %getWavelengthVec(lambda0,fracBW,numWavelengths);% array of wavelengths (meters)
% 
% % Define charge of the vortex mask at each wavelength
% %charge = ones(1,numWavelengths); % achromatic
% charge = 1 *(lambda0./lambdas); % simple scalar model
% 
% % Define wavefront error at the central wavelength
% nolls = 4:8;
% coeffs = 0*[0.01,-0.0099,-0.0095,-0.0008,0.0033];
% 
% % Give offsets for the vortex mask
% offsetX = 0;%0.0952*apRad;
% offsetY = 0;%0.0524*apRad; 
% 
% %% Generate the coordinate system
% 
% coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
% xvals = coords.xvals;% Helpful for plotting
% yvals = coords.yvals;
% 
% %% Create array with pupil function
% 
% PUPIL = makeKeckPupil(2*apRad, N );
% [normI, totalPower0] = getNormalization(PUPIL);% Normalization factors
% lambdaOverD = N/apRad/2; % lam/D in units of samples in the image plane
% 
%  figure(1)
%  imagesc(xvals/apRad,yvals/apRad,PUPIL);
%  axis image; 
%  axis([-1 1 -1 1]);
%  title('Pupil');
%  colorbar; 
%  colormap(parula(256));
%  grid on;
%  drawnow;
% 
% %% Define pupil field
% 
% %phz = generateZernike_fromList( nolls, coeffs, PUPIL, apRad, coords);
% for ch = 1:numWavelengths
%     phz(:,:,ch) = angle(makeKeckPupilPhz(2*apRad,N,charge(ch)));
% end
% %phz = angle(makeKeckPupilPhase(2*apRad,N,chargeC));
% % phz2 = angle(makeKeckPupilField(2*apRad,N));
% 
% figure(2);
% for ch = 1:numWavelengths
%     Epup(:,:,ch) = exp(1i*phz(:,:,ch)*lambda0/lambdas(ch)).*PUPIL;
%     
%     subplot(1,numWavelengths,ch);
%     imagesc(xvals/apRad,yvals/apRad,angle(Epup(:,:,ch)));
%     axis image; 
%     axis xy;
%     axis([-1 1 -1 1]);
%     title(['Phase at ',num2str(lambdas(ch)*1e9),'nm']);
%     colorbar; 
%     colormap(hsv(256));
%    
% end
% drawnow;




