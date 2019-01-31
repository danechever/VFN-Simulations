clear; %close all; 
addpath('../VFN_simulations/VFNlib/');

N = 2^12; % Size of computational grid (NxN samples) 
apRad = N/16; % Aperture radius in samples 
charge = 1; % Charge of the vortex mask 

coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
xvals = coords.xvals;% Helpful for plotting
yvals = coords.yvals;

showPlots = false;
useGPU = false; % make 'true' if you have a useful GPU card 
%% Create array with pupil function

EP = makeCircularPupil( apRad, N );

if(showPlots)
    figure(1)
    imagesc(xvals/apRad,yvals/apRad,EP);
    axis image; 
    axis([-1 1 -1 1]);
    title('Pupil');
    colorbar; 
    colormap(parula(256));
end

%% Get PSF without vortex mask and normalization factors 

% Keeps the origin at the center of the array (keep N even)
myfft2 = @(x) fftshift(fft2(fftshift(x)));
myifft2 = @(x) fftshift(ifft2(fftshift(x)));

PSF = myfft2(EP); % PSF with a flat wavefront 

normI = max(max(abs(PSF).^2));% Normalization for PSF plots
totalPower0 = sum(sum(abs(PSF).^2));% Normalization for coupling fractions
lambdaOverD = N/apRad/2; % lam/D in units of samples in the image plane

if(showPlots)
    figure(2)
    imagesc(xvals/lambdaOverD,yvals/lambdaOverD,abs(PSF).^2/normI);
    axis image; 
    axis([-3 3 -3 3]);
    title('PSF w/o vortex');
    colorbar; 
    colormap(parula(256));
end

%% Make vortex mask 

offsetX = 0;
offsetY = 0; 

EPM = generateVortexMask( charge, coords, [offsetX offsetY] );

if(showPlots)
    figure(3)
    imagesc(xvals/apRad,yvals/apRad,angle(EPM.*EP));
    axis image; 
    axis([-1 1 -1 1]);
    title('Pupil phase');
    colormap(hsv(256));
    colorbar; 
end

%% Add Zernike aberration 

% Fiber diameter
fiberDiam = 1.4; % units of lambda/D

% Noll indicies to analyze
noll_indices = 2:105;
% noll_indices = [2,3,7,8,16,17,29,30,46,47,67,68,92,93];

% List of Zernike coefficients 
coeffs = logspace(-4,-2,-2-(-4)+1);% coefficients are in units of waves rms 

% 2D array with fiber mode 
fibermode0 = generateFiberMode(fiberDiam*lambdaOverD,coords);

disp('j     power b');

if(useGPU)
    EP = gpuArray(EP);
end

count = 1;
for noll_index = noll_indices % Noll index
    eta_s_list = [];

    for coeff = coeffs % waves rms
        
        [Z,n,m] = generateZernike(noll_index,apRad,coords.RHO,coords.THETA);
        Z = Z/sqrt(mean(Z(logical(EP)).^2)); % Re-normalize (useful when pupil is not a circle)

        ABER = exp(1i*2*pi*coeff*Z);
        PSFv = myfft2(EP.*EPM.*ABER);

        eta_s = abs(sum(sum(fibermode0.*PSFv))).^2/totalPower0;
        eta_s_list = [eta_s_list,eta_s];
    end
    
    p = polyfit(log10(coeffs),log10(eta_s_list),1);
    b = 10^(p(2)/2);
   
    
    disp([num2str(noll_index),'     ',num2str(p(1),2),'     ',num2str(b,2)]);
    M(count,:) = [noll_index,n,m,p(1),b];
    count = count + 1;
end

%%
if(useGPU)
    M = gather(M);
end
T = array2table(M,'VariableNames',{'noll_index','n','m','power','b'});
writetable(T,'VFNzernikeSensitivities.csv');
