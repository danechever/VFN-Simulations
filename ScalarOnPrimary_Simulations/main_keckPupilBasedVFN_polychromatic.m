clear; close all;
addpath(['..' filesep 'VFNlib']);

%% Input parameters 

% Define smapling info
inputs.N = 2^11; % Size of computational grid (NxN samples) 
inputs.apRad = 128; % Aperture radius in samples 
inputs.apDia0 = 2 * inputs.apRad;

% Define wavelength info
inputs.lambda0 = 2.2e-6; %central wavelength
inputs.fracBW = 0.1818; %\Delta\lambda/\lambda
inputs.numWavelengths = 5;% number of discrete wavelengths 
inputs.lambdas = getWavelengthVec(inputs.lambda0,inputs.fracBW,inputs.numWavelengths);% array of wavelengths (meters)

% Define charge of the vortex mask at each wavelength
%charge = ones(1,numWavelengths); % achromatic
inputs.charge = 1*(inputs.lambda0./inputs.lambdas); % simple scalar model

% Define wavefront error at the central wavelength
inputs.nolls = 4:8;
inputs.coeffs = 0*[0.01,-0.0099,-0.0095,-0.0008,0.0033];

% Give offsets for the vortex mask
inputs.offsetX = 0;%0.0952*apRad;
inputs.offsetY = 0;%0.0524*apRad; 

inputs.numRings = 3;
inputs.wGap = 25.4/10916*inputs.apDia0/2;

%% Generate the coordinate system

coords = generateCoordinates(inputs.N);% Creates NxN arrays with coordinates 
xvals = coords.xvals;% Helpful for plotting
yvals = coords.yvals;

%% Create matrix with tilt and pistons for the analytical solution, generate constant coefficients for phase and amplitude for each hexagonal segment.
hexFlatDiam = (inputs.apDia0-3*2*inputs.wGap)/(2*3+1);
inputs.hexSep = hexFlatDiam + inputs.wGap;

inputs.hexAmpConst = NaN(inputs.N, inputs.N, 36);
inputs.hexPhzConst = NaN(inputs.N, inputs.N, 36);

initial = NaN(36,3);
 loc = 1;
 for ringNum = 1:3 %adds absolute value of tilt to each segment in the position they are added.
     for seg = (loc):(loc+ringNum*6-1)
         initial(seg,1:2) = (1*2*pi)/(ringNum*6);
     end
     loc = (loc+ringNum*6);
 end

segs = ones(1,36);
count = 1;
 for ringNum = 1:3
     crow = ringNum * inputs.hexSep;
     ccol = 0;
     t = atan2(crow,ccol);
     t = round(t,3);
     
     if(segs(count) == 1)
           %Uncomment below to use analytical solution
         initial(count,1) = initial(count,1) * -sin(t);
         initial(count,2) = initial(count,2) * cos(t);
         initial(count,3) = t/(2*pi);
         
         [inputs.hexAmpConst(:,:,count), inputs.hexPhzConst(:,:,count)] = generateHexConstants(crow, ccol, inputs.numRings, inputs.apDia0, inputs.wGap, zeros(inputs.N));
         
     end

     count = count + 1;
     for face = 1:6
         step_dir = pi/6*(2*face+5);
         steprow = inputs.hexSep*sin(step_dir);
         stepcol = inputs.hexSep*cos(step_dir);
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
                     %Uncomment below to use analytical solution
                    initial(count,1) = initial(count,1) * -sin(t);
                    initial(count,2) = initial(count,2) * cos(t);
                    initial(count,3) = t/(2*pi);
                    
                    [inputs.hexAmpConst(:,:,count), inputs.hexPhzConst(:,:,count)] = generateHexConstants(crow, ccol, inputs.numRings, inputs.apDia0, inputs.wGap, zeros(inputs.N));

                 end
             count = count + 1;
             end
             stepnum = stepnum + 1;
         end
     end
 end

%% Create array with pupil function

PUPIL = makeKeckPupil(2*inputs.apRad, inputs.N );
[normI, totalPower0] = getNormalization(PUPIL);% Normalization factors
inputs.lambdaOverD = inputs.N/inputs.apRad/2; % lam/D in units of samples in the image plane

 figure(1)
 imagesc(xvals/inputs.apRad,yvals/inputs.apRad,PUPIL);
 axis image; 
 axis([-1 1 -1 1]);
 title('Pupil');
 colorbar; 
 colormap(parula(256));
 grid on;
 drawnow;

 addpath(['..' filesep '..' filesep 'falco-matlab' filesep 'lib' filesep 'utils']);
  %-- Decrease matrix size in pupil plane to reduce runtime
 inputs.PUP_CRP_SZ = round(2.1*inputs.apRad);
 inputs.hexAmpConst = pad_crop(inputs.hexAmpConst,inputs.PUP_CRP_SZ);
 inputs.hexPhzConst = pad_crop(inputs.hexPhzConst,inputs.PUP_CRP_SZ);
%% Define pupil field

%phz = generateZernike_fromList( inputs.nolls, inputs.coeffs, PUPIL, inputs.apRad, coords);
%optSeed = load('exampleFile.mat','optimum');
%disp(optSeed);

phz = angle(makeKeckPupilInputs( inputs, initial));
%phz(:,:,ch) = angle(makeKeckPupilPhz(inputs.apDia0, inputs.N, inputs.charge));
%phz = angle(makeKeckPupilPhase(2*apRad,N,chargeC));
%phz2 = angle(makeKeckPupilField(2*apRad,N));

figure(2);
for ch = 1:inputs.numWavelengths
    Epup(:,:,ch) = exp(1i*phz*inputs.lambda0/inputs.lambdas(ch)).*PUPIL;
    
    subplot(1,inputs.numWavelengths,ch);
    imagesc(xvals/inputs.apRad,yvals/inputs.apRad,angle(Epup(:,:,ch)));
    axis image; 
    axis xy;
    axis([-1 1 -1 1]);
    title(['Phase at ',num2str(inputs.lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(hsv(256));
   
end
drawnow;

% figure(10);
% for ch = 1:numWavelengths
%     Epup(:,:,ch) = exp(1i*phz2*lambda0/lambdas(ch)).*PUPIL;
%     
%     subplot(1,numWavelengths,ch);
%     imagesc(xvals/apRad,yvals/apRad,angle(Epup(:,:,ch)));
%     axis image; 
%     axis([-1 1 -1 1]);
%     title(['Phase at ',num2str(lambdas(ch)*1e9),'nm']);
%     colorbar; 
%     colormap(hsv(256));
% end
% drawnow;

% Get PSF without vortex mask

%% Get broadband PSF
iPSF_BB = getPSF(Epup,inputs.lambda0,inputs.lambdas,normI,coords);

figure(3)
imagesc(xvals/inputs.lambdaOverD,yvals/inputs.lambdaOverD,iPSF_BB);
axis image; 
axis([-2 2 -2 2]);
title('broadband PSF w/o vortex');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Make vortex mask 

% EPM = generateVortexMask( inputs.charge, coords, [inputs.offsetX inputs.offsetY] );
% 
% central_band_index = ceil(inputs.numWavelengths/2);
% 
% figure(4)
% imagesc(xvals/inputs.apRad,yvals/inputs.apRad,angle(EPM(:,:,central_band_index).*Epup(:,:,central_band_index)));
% axis image; 
% axis([-1 1 -1 1]);
% title('Pupil phase at \lambda_0');
% colormap(hsv(256));
% colorbar; 
% drawnow;

%% Get PSF with vortex mask

% iPSFv_BB = getPSF(Epup,lambda0,lambdas,normI,coords);
% 
% figure(5)
% imagesc(xvals/lambdaOverD,yvals/lambdaOverD,iPSFv_BB);
% axis image; 
% axis([-2 2 -2 2]);
% title('broadband PSF w/ vortex');
% colorbar;%caxis([-3 0])
% colormap(parula(256));
% drawnow;

%% Generate coupling maps for each wavelength

% Parameters for Thorlabs SM600
%Using SM2000 in this version of the code
    % link: https://www.thorlabs.com/NewGroupPage9_PF.cfm?ObjectGroup_ID=949
fiber_props.core_rad = 5.5e-6;% Core radius [um]
fiber_props.n_core = 1.4436; %1.4571; core index (interpolated from linear fit to 3 points)
fiber_props.n_clad = 1.4381; %1.4558; cladding index (interpolated from linear fit to 3 points)
fiber_props.type = 'bessel';
Fnum = getMFD(fiber_props,inputs.lambda0)/(inputs.lambda0*1.4); % focal ratio of the beam at the fiber

eta_maps = generateCouplingMap_polychromatic( Epup, fiber_props, inputs.lambda0, Fnum, inputs.lambdas, totalPower0, inputs.lambdaOverD, 3*inputs.lambdaOverD, coords);

figure(6);
for ch = 1:inputs.numWavelengths
    subplot(1,inputs.numWavelengths,ch);
    imagesc(xvals/inputs.lambdaOverD,yvals/inputs.lambdaOverD,log10(eta_maps(:,:,ch)));
    axis image; 
    axis([-2 2 -2 2]);
    caxis([-3 -0.5])
    title(['\eta at ',num2str(inputs.lambdas(ch)*1e9),'nm']);
    colorbar;
    colormap(gray(256));
end

% find the centroid of eta_maps(:,:,ch)
% find min in centroid
% find coordinates of this min in eta_maps
% plot the difference in these coordinates in x and y position from the
% origin, with lambda corresponding to ch on the bottom and offset on the
% left
% 
% Finds the radius of the centroid to crop the image to depending on the maximum eta in 1 slice
% of eta_maps(:,:,ch)
Xshift = zeros(inputs.numWavelengths,1);
Yshift = zeros(inputs.numWavelengths,1);
etas = zeros(inputs.numWavelengths, 1);
%etas_offset = zeros(numWavelengths, 1);
for ch = 1:inputs.numWavelengths 
    map = eta_maps(:,:,ch); %one slice of the eta_maps cube
    map_max = max(map,[],'all'); %the maximum value in cmap
    [max_ind(1),max_ind(2)] = find(map_max==map,1); %linear coordinates of max value
    max_rho = sqrt(((inputs.N/2+1)-max_ind(1))^2 + ((inputs.N/2+1)-max_ind(2))^2);
    
    crp = 2*max_rho; %The length of one side of the cube to crop the image to

    cmap = map(end/2+1-floor(crp/2):end/2+1+floor(crp/2),end/2+1-floor(crp/2):end/2+1+floor(crp/2)); %the centroid
    cmap_min = min(cmap,[],'all'); %minimum value in the centroid
    [min_ind(1),min_ind(2)] = find(cmap_min==cmap); %indices of minimum value in the centroid
    min_ind = min_ind + (inputs.N/2-floor(crp/2)); %adjust min values to reflect position in map, not cmap
   
    Xshift(ch) = inputs.N/2-min_ind(2); %x value is wavelength, y value is offset
    Yshift(ch) = inputs.N/2-min_ind(1);%^
    
    etas(ch) = cmap_min;
    %etas_offset(ch) = sqrt((Xshift(ch))^2 + (Yshift(ch)^2));
end

figure(7);
subplot(2,2,1);
plot(inputs.lambdas/inputs.lambda0,Xshift/inputs.lambdaOverD, '-o', 'Color', 'r');
title('Xshift');
subplot(2,2,2);
plot(inputs.lambdas/inputs.lambda0,Yshift/inputs.lambdaOverD, '-o', 'Color', 'b');
title('Yshift');

px = polyfit(inputs.lambdas/inputs.lambda0,Xshift'/inputs.lambdaOverD,1);
pxy = polyval(px,inputs.lambdas/inputs.lambda0);
subplot(2,2,3);
plot(inputs.lambdas/inputs.lambda0,pxy,'-o','Color','m')
title('X Offset Trend')
txt = ['p value: ' num2str(px)];
text(mean(inputs.lambdas/inputs.lambda0),mean(pxy),txt);

py = polyfit(inputs.lambdas/inputs.lambda0,Yshift'/inputs.lambdaOverD,1);
pyy = polyval(py,inputs.lambdas/inputs.lambda0);
subplot(2,2,4);
plot(inputs.lambdas/inputs.lambda0,pyy,'-o','Color','g');
title('Y Offset Trend')
txt = ['p value: ' num2str(py)];
text(mean(inputs.lambdas/inputs.lambda0),mean(px),txt);

figure(8);
subplot(1,1,1);
semilogy(inputs.lambdas/inputs.lambda0,etas,'-o','Color','r'); %lambdas/lambda0,,'-o','Color','r'
title('Null Value vs \lambda/\lambda0')
xlabel('\lambda/\lambda0')
ylabel('\eta')
grid on

figure(9);
plot(inputs.lambdas/inputs.lambda0,Xshift/inputs.lambdaOverD, '-o', 'Color', 'r');
hold on
plot(inputs.lambdas/inputs.lambda0,Yshift/inputs.lambdaOverD, '-o', 'Color', 'b');
legend({'Xshift', 'Yshift'}, 'Location', 'SouthEast');
title('Null Movement')
xlabel('\lambda/\lambda_0')
ylabel('Null Shift [\lambda/D]')
grid on

figure(12);
plot(inputs.lambdas/inputs.lambda0,pxy,'-o','Color','r');
hold on
plot(inputs.lambdas/inputs.lambda0,pyy,'-o','Color','b');
legend({'Xshift', 'Yshift'}, 'Location', 'SouthEast');
title('Null Movement')
xlabel('\lambda/\lambda_0')
ylabel('Null Shift [\lambda/D]')
txt = ['y trend p value: ' num2str(py) newline 'x trend p value: ' num2str(px)];
text(mean(inputs.lambdas/inputs.lambda0),mean(px),txt);
grid on





