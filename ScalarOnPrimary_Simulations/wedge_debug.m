clear; close all;
addpath(['..' filesep 'VFNlib']);
py.importlib.import_module('sellmeir1');
py.importlib.import_module('triple_prism');
load pyyo;
load oldetas;
sellmeir1coeffs;
%% Input parameters

inpar.pyyCopy = pyyCopy;
inpar.etasCopy = etasCopy;

% Define sampling info
inpar.N = 2^12; % Size of computational grid (NxN samples) 
inpar.apRad = 128; % Aperture radius in samples 
inpar.apDia0 = 2 * inpar.apRad;
inpar.keckD = 10.949; %Meters

% Define wavelength info
inpar.lambda0 = 2.2e-6; %central wavelength
inpar.fracBW = 0.1818; %\Delta\lambda/\lambda
inpar.numWavelengths = 5;% number of discrete wavelengths 
inpar.lambdas = getWavelengthVec(inpar.lambda0,inpar.fracBW,inpar.numWavelengths);% array of wavelengths (meters)
inpar.lam0OverD_rad = inpar.lambdas(ceil(inpar.numWavelengths / 2)) / inpar.keckD;

%Define charge of the vortex mask at each wavelength
%charge = ones(1,numWavelengths); % achromatic
inpar.charge = 1*(inpar.lambda0./inpar.lambdas); % simple scalar model

% Define wavefront error at the central wavelength
inpar.nolls = 4:8;
inpar.coeffs = 0*[0.01,-0.0099,-0.0095,-0.0008,0.0033];

% Give offsets for the vortex mask
inpar.offsetX = 0;%0.0952*apRad;
inpar.offsetY = 0;%0.0524*apRad; 

inpar.numRings = 3;
inpar.wGap = 25.4/10916*inpar.apDia0/2;

%Wedge Data
%inPar.lin_offset = py_1Copy;
wType = 'wedge';

inpar.p = 1.1001;%inPar.lin_offset(1);%1.1404;
inpar.wedge_mat = 'CaF2';
% wedge_angle = atan(890.16*(inpar.p*inpar.lambda0)/(inpar.keckD*(getRefractiveIndex(inpar.wedge_mat ,1e6*inpar.lambda0) - 1)));
wedge_angle = pi/(180*3600);

%ADC Wedge Data
wedge_ADC_1 = 'BaF2';
wedge_ADC_2 = 'CaF2';
wedge_ADC_3 = 'ZnSe';

inpar.phi1_ADC = 7.0516 * pi/180;
inpar.phi2_ADC = 3.8050 * pi/180;
inpar.phi3_ADC = 1.1465 * pi/180;

wvs = py.numpy.array(1e6*inpar.lambdas);
inpar.n0 = 0;
inpar.n1 = py.sellmeir1.sellmeir1(wvs, 273.15 + 3, baf2_args(1,:),baf2_args(2,:),baf2_args(3,:),baf2_args(4,:),baf2_args(5,:),baf2_args(6,:));
inpar.n2 = py.sellmeir1.sellmeir1(wvs, 273.15 + 3, caf2_args(1,:),caf2_args(2,:),caf2_args(3,:),caf2_args(4,:),caf2_args(5,:),caf2_args(6,:));
inpar.n3 = py.sellmeir1.sellmeir1(wvs, 273.15 + 3, znse_args(1,:),znse_args(2,:),znse_args(3,:),znse_args(4,:),znse_args(5,:),znse_args(6,:));

inpar.clocking = deg2rad(76.3636);

inpar.tilt1 = 0;

% Provide path where the DM Basis file is located
dmBasisPath = '/media/Data_Drive/VFN/ScalarOnPrimaryData/';

%% Generate the coordinate system

coords = generateCoordinates(inpar.N);% Creates NxN arrays with coordinates
inpar.Ycoords = coords.Y;
inpar.Xcoords = coords.X;
inpar.xvals = coords.xvals;% Helpful for plotting
inpar.yvals = coords.yvals;

%% Create array with pupil function

% inpar.PUPIL = makeKeckPupil(2*inpar.apRad, inpar.N );
inpar.PUPIL = makeCircularPupil(inpar.apRad,inpar.N);
[normI, totalPower0] = getNormalization(inpar.PUPIL); % Normalization factors
inpar.lambdaOverD_cgrid = inpar.N/inpar.apRad/2; % lam/D in units of samples in the image plane
inpar.lambdaOverD = inpar.N/inpar.apRad/2; % lam/D in units of samples in the image plane

% Sinf.PUP_CRP_SZ = round(2.1*Sinf.apRad);
% Sinf.hexAmpConst = pad_crop(Sinf.hexAmpConst,Sinf.PUP_CRP_SZ);
% Sinf.hexPhzConst = pad_crop(Sinf.hexPhzConst,Sinf.PUP_CRP_SZ);

figure(1)
imagesc(inpar.xvals/inpar.apRad,inpar.yvals/inpar.apRad,inpar.PUPIL);
axis image; 
axis([-1 1 -1 1]);
title('Pupil');
colorbar; 
colormap(parula(256));
grid on;
drawnow;

%% Define Pupil Field without Vortex Mask
% phz = generateVortexMaskKeckPrimary(inPar);
%phz = generateDMVortex(dmBasisPath);
phz = generateZernike_fromList( inpar.nolls, inpar.coeffs, inpar.PUPIL, inpar.apRad, coords);

figure(2);
for ch = 1:inpar.numWavelengths

    Epup(:,:,ch) = exp(1i*phz*inpar.lambda0/inpar.lambdas(ch)).*inpar.PUPIL;

    subplot(1,inpar.numWavelengths,ch);
    imagesc(inpar.xvals/inpar.apRad,inpar.yvals/inpar.apRad,angle(Epup(:,:,ch)));
    axis image; 
    axis xy;
    axis([-1 1 -1 1]);
    title(['Phase at ',num2str(inpar.lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(hsv(256));
end
drawnow;
%% Get broadband PSF 
%iPSF_BB = getPSF(Epups,inPar.lambda0,inPar.lambdas,normI,coords);
%iPSF_BB = mean(iPSF_BB, 3);

PSF = getPSF(Epup,inpar.lambda0,inpar.lambdas,normI,coords);

figure(3)
imagesc(inpar.xvals/inpar.lambdaOverD_cgrid,inpar.yvals/inpar.lambdaOverD_cgrid,PSF);
axis image; 
axis([-2 2 -2 2]);
title('Broadband PSF w/o vortex & wedge');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Apply Wedge Effects

figure(4);
inpar.lam0OverD_rad = inpar.lambda0/inpar.keckD;

for ch = 1:inpar.numWavelengths
    
    disp(getRefractiveIndex('CaF2',inpar.lambdas(ch)*1e6));
    
    wphz = 2*pi*(getWedgeTilt(inpar.wedge_mat, wedge_angle, 1e6*inpar.lambdas(ch))/inpar.lam0OverD_rad*inpar.lambda0/inpar.lambdas(ch))...
        *inpar.lambdaOverD*inpar.Ycoords/inpar.N; %- 2*pi*inpar.p*inpar.lambdaOverD*inpar.Ycoords/inpar.N;
    
    phzW = phz - wphz;
    
    %Epup(:,:,ch) = exp(1i*phz*inPar.lambda0/inPar.lambdas(ch)).*inPar.PUPIL;
    Epupw(:,:,ch) = exp(1i*wphz*inpar.lambda0/inpar.lambdas(ch)).*inpar.PUPIL;
    Epup_Wedge(:,:,ch) = exp(1i*phzW*inpar.lambda0/inpar.lambdas(ch)).*inpar.PUPIL;
    
    subplot(1,inpar.numWavelengths,ch);
    imagesc(inpar.xvals/inpar.apRad,inpar.yvals/inpar.apRad,angle(Epupw(:,:,ch)));
    axis image; 
    axis xy;
    axis([-1 1 -1 1]);
    title(['Phase at ',num2str(inpar.lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(hsv(256));
end



% wphz = generateWedgePlate(inpar,wedge_angle,inpar.lambdas(ceil(inpar.numWavelengths/2)),'wedge');

% figure(3);
% if(wType == 'ADC__')
%     inpar.lam0OverD_rad = inpar.lambda0/inpar.keckD;
%     wphz = generateWedgePlate(inpar,wedge_angle,inpar.lambdas(ceil(inpar.numWavelengths/2)),'ADC__');
%     
%     for ch =1:inpar.numWavelengths
%                  
%         phzW = phz - wphz(:,:,ch);%phz + wphz(ch);
%     
%         %Epup(:,:,ch) = exp(1i*phz*inPar.lambda0/inPar.lambdas(ch)).*inPar.PUPIL;
%         Epupw(:,:,ch) = exp(1i*wphz(:,:,ch)*inpar.lambda0/inpar.lambdas(ch)).*inpar.PUPIL;
%         Epup_Wedge(:,:,ch) = exp(1i*phzW*inpar.lambda0/inpar.lambdas(ch)).*inpar.PUPIL;
%     
%         subplot(1,inpar.numWavelengths,ch);
%         imagesc(inpar.xvals/inpar.apRad,inpar.yvals/inpar.apRad,angle(Epupw(:,:,ch)));
%         axis image; 
%         axis xy;
%         axis([-1 1 -1 1]);
%         title(['Phase at ',num2str(inpar.lambdas(ch)*1e9),'nm']);
%         colorbar; 
%         colormap(hsv(256));
%     end
% elseif(wType == 'wedge')
%     
%     lam0OverD_rad = inpar.lambda0/inpar.keckD;
%     disp(pyyCopy * lam0OverD_rad);
% 
%     for ch = 1:inpar.numWavelengths
%     
%         inpar.lam0OverD_rad = inpar.lambda0/inpar.keckD;
%     
%         %Requisite difference in phase computed here
%         wphz = generateWedgePlate(inpar,wedge_angle,inpar.lambdas(ch),'wedge');
%     
%         phzW = phz - wphz;
%     
%         %Epup(:,:,ch) = exp(1i*phz*inPar.lambda0/inPar.lambdas(ch)).*inPar.PUPIL;
%         Epupw(:,:,ch) = exp(1i*wphz*inpar.lambda0/inpar.lambdas(ch)).*inpar.PUPIL;
%         Epup_Wedge(:,:,ch) = exp(1i*phzW*inpar.lambda0/inpar.lambdas(ch)).*inpar.PUPIL;
%     
%         subplot(1,inpar.numWavelengths,ch);
%         imagesc(inpar.xvals/inpar.apRad,inpar.yvals/inpar.apRad,angle(Epupw(:,:,ch)));
%         axis image; 
%         axis xy;
%         axis([-1 1 -1 1]);
%         title(['Phase at ',num2str(inpar.lambdas(ch)*1e9),'nm']);
%         colorbar; 
%         colormap(hsv(256));
%     end
% end

%% Get wavelength specific PSF 
%iPSF_BB = getPSF(Epups,inPar.lambda0,inPar.lambdas,normI,coords);
%iPSF_BB = mean(iPSF_BB, 3);

PSFw = getallPSF(Epup_Wedge,inpar.lambda0,inpar.lambdas,normI,coords);

figure(5)
for ch = 1:inpar.numWavelengths
    
    subplot(1,inpar.numWavelengths,ch);
    imagesc(inpar.xvals/inpar.lambdaOverD_cgrid,inpar.yvals/inpar.lambdaOverD_cgrid,PSFw(:,:,ch));
    axis image;
    axis xy;
    axis([-2 2 -2 2]);
    title(['PSF at ' num2str(inpar.lambdas(ch)*1e9) 'nm w/ wedge']);
    colorbar;%caxis([-3 0])
    colormap(parula(256));
end
drawnow;

%% Generate coupling maps for each wavelength

%-- Parameters for Thorlabs SM600
%-- Using SM2000 in this version of the code
    %-- link: https://www.thorlabs.com/NewGroupPage9_PF.cfm?ObjectGroup_ID=949
fiber_props.core_rad = 5.5e-6;% Core radius [um]
fiber_props.n_core = 1.4436; %1.4571; core index (interpolated from linear fit to 3 points)
fiber_props.n_clad = 1.4381; %1.4558; cladding index (interpolated from linear fit to 3 points)
fiber_props.type = 'bessel';
Fnum = getMFD(fiber_props,inpar.lambda0)/(inpar.lambda0*1.4); % focal ratio of the beam at the fiber

eta_maps = generateCouplingMap_polychromatic(Epup_Wedge, fiber_props, inpar.lambda0, Fnum, inpar.lambdas, totalPower0, inpar.lambdaOverD, 3*inpar.lambdaOverD, coords);

figure();
for ch = 1:inpar.numWavelengths
    subplot(1,inpar.numWavelengths,ch);
    imagesc(inpar.xvals/inpar.lambdaOverD,inpar.yvals/inpar.lambdaOverD,log10(eta_maps(:,:,ch)));
    axis image; 
    axis([-2 2 -2 2]);
    caxis([-3 -0.5])
    title(['\eta at ',num2str(inpar.lambdas(ch)*1e9),'nm']);
    colorbar;
    colormap(gray(256));
end

%-- find the centroid of eta_maps(:,:,ch)
%-- find min in centroid
%-- find coordinates of this min in eta_maps
%-- plot the difference in these coordinates in x and y position from the
%   origin, with lambda corresponding to ch on the bottom and offset on the
%   left
% 
%-- Finds the radius of the centroid to crop the image to depending on the maximum eta in 1 slice
%   of eta_maps(:,:,ch)
Xshift = zeros(inpar.numWavelengths,1);
Yshift = zeros(inpar.numWavelengths,1);
etas = zeros(inpar.numWavelengths, 1);

for ch = 1:inpar.numWavelengths 
    map = eta_maps(:,:,ch); %one slice of the eta_maps cube
    map_max = max(map,[],'all'); %the maximum value in cmap
    [max_ind(1),max_ind(2)] = find(map_max==map,1); %linear coordinates of max value
    max_rho = sqrt(((inpar.N/2+1)-max_ind(1))^2 + ((inpar.N/2+1)-max_ind(2))^2);
    
    crp = 2*max_rho; %The length of one side of the cube to crop the image to

    cmap = map(end/2+1-floor(crp/2):end/2+1+floor(crp/2),end/2+1-floor(crp/2):end/2+1+floor(crp/2)); %the centroid
    cmap_min = min(cmap,[],'all'); %minimum value in the centroid
    [min_ind(1),min_ind(2)] = find(cmap_min==cmap); %indices of minimum value in the centroid
    min_ind = min_ind + (inpar.N/2-floor(crp/2)); %adjust min values to reflect position in map, not cmap
   
    Xshift(ch) = inpar.N/2+1-min_ind(2); %x value is wavelength, y value is offset
    Yshift(ch) = inpar.N/2+1-min_ind(1);%^
    
    etas(ch) = cmap_min;
end

%%
%-- Null shift plots for X, Y, and the trendlines that result
figure();
subplot(2,2,1);
plot(inpar.lambdas/inpar.lambda0,Xshift/inpar.lambdaOverD, '-o', 'Color', 'r');
title('Xshift');
subplot(2,2,2);
plot(inpar.lambdas/inpar.lambda0,Yshift/inpar.lambdaOverD, '-o', 'Color', 'b');
title('Yshift');

px = polyfit(inpar.lambdas/inpar.lambda0,Xshift'/inpar.lambdaOverD,1);
pxy = polyval(px,inpar.lambdas/inpar.lambda0);
subplot(2,2,3);
plot(inpar.lambdas/inpar.lambda0,pxy,'-o','Color','m')
title('X Offset Trend')
txt = ['p value: ' num2str(px)];
text(mean(inpar.lambdas/inpar.lambda0),mean(pxy),txt);

py_1 = polyfit(inpar.lambdas/inpar.lambda0,Yshift'/inpar.lambdaOverD,1);
pyy = polyval(py_1,inpar.lambdas/inpar.lambda0);
subplot(2,2,4);
plot(inpar.lambdas/inpar.lambda0,pyy,'-o','Color','g');
title('Y Offset Trend')
txt = ['p value: ' num2str(py_1)];
text(mean(inpar.lambdas/inpar.lambda0),mean(px),txt);

%%
%-- Null value vs wavelength offset from central wavelength
figure();
subplot(1,1,1);
semilogy(inpar.lambdas/inpar.lambda0,etas,'-o','Color','r'); %lambdas/lambda0,,'-o','Color','r'
title('Null Value vs \lambda/\lambda0')
xlabel('\lambda/\lambda0')
ylabel('\eta')
grid on

%%
%-- Actual positional offset of null for x and y overlayed
figure();
plot(inpar.lambdas/inpar.lambda0,Xshift/inpar.lambdaOverD, '-o', 'Color', 'r');
hold on
plot(inpar.lambdas/inpar.lambda0,Yshift/inpar.lambdaOverD, '-o', 'Color', 'b');
legend({'Xshift', 'Yshift'}, 'Location', 'SouthEast');
title('Null Movement')
xlabel('\lambda/\lambda_0')
ylabel('Null Shift [\lambda/D]')
grid on

%%
%-- Overlay of trends in x and y null positional offset
figure();
plot(inpar.lambdas/inpar.lambda0,pxy,'-o','Color','r');
hold on
plot(inpar.lambdas/inpar.lambda0,pyy,'-o','Color','b');
legend({'Xshift', 'Yshift'}, 'Location', 'SouthEast');
title('Null Movement')
xlabel('\lambda/\lambda_0')
ylabel('Null Shift [\lambda/D]')
txt = ['y trend p value: ' num2str(py_1) newline 'x trend p value: ' num2str(px)];
%text(mean(inPar.lambdas/inPar.lambda0),mean(px),txt);
grid on