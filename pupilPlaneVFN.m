clear; close all; 

N = 2^12; % Size of computational grid (NxN) 
apRad = 250; % Aperture radius in samples 
charge = 1; % Charge of the vortex mask 

useGPU = false; % Only true if you know what you are doing 

%% Create array with pupil function

EP = makeCircularPupil( apRad, N );

lambdaOverD = N/apRad/2; % samples 

% Create coordinate system 
[X,Y] = meshgrid(-N/2:N/2-1);
[THETA,RHO] = cart2pol(X,Y);
xvals = X(1,:);yvals = Y(:,1);

figure(1)
imagesc(xvals/apRad,yvals/apRad,EP);
axis image; 
axis([-1 1 -1 1]);
title('Pupil');
colorbar; 
colormap(parula(256));
%% Get PSF without vortex mask and normalization factors 

myfft2 = @(x) fftshift(fft2(fftshift(x)));

PSF = myfft2(EP);
normI = max(max(abs(PSF).^2));
totalPower0 = sum(sum(abs(PSF).^2));

figure(2)
imagesc(xvals/lambdaOverD,yvals/lambdaOverD,abs(PSF).^2/normI);
axis image; 
axis([-3 3 -3 3]);
title('PSF w/o vortex');
colorbar; 
colormap(parula(256));
%% Make vortex mask 

EPM = exp(1i*charge*THETA);

figure(3)
imagesc(xvals/apRad,yvals/apRad,angle(EPM.*EP));
axis image; 
axis([-1 1 -1 1]);
title('Pupil phase');
colormap(hsv(256));
colorbar; 

%% Get PSF with vortex mask

PSFv = myfft2(EP.*EPM);

figure(4)
imagesc(xvals/lambdaOverD,yvals/lambdaOverD,abs(PSFv).^2/normI);
axis image; 
axis([-3 3 -3 3]);
title('log10(PSF) w/ vortex');
colorbar; 
colormap(parula(256));

%% Get throughput curve for different F numbers 

angSeps = (0:0.1:2)*lambdaOverD;
fiberDiams = [1 1.4 1.8];
figure(5); 
for fiberDiam = fiberDiams
    
    fibermode0 = sqrt(2/(pi*(fiberDiam*lambdaOverD/2)^2))* ...
        exp(-(RHO/(fiberDiam*lambdaOverD/2)).^2);
    
    coupling_planet = [];
    for angSep = angSeps
        tilt = exp(1i*2*pi*angSep*X/N);

        FPp = myfft2(EP.*EPM.*tilt);
        coupling_planet = [coupling_planet,abs(sum(sum(fibermode0.*FPp))).^2/totalPower0];
    end

    plot(angSeps/lambdaOverD,coupling_planet);hold on;
    xlabel('Angular separation (\lambda/D)');
    ylabel('\eta_p');
    drawnow;
end
hold off;
legLabels = {};
for index = 1:numel(fiberDiams)
    legLabels{index} = ['D_f = ',num2str(fiberDiams(index)),' \lambda F#'];
end
legend(legLabels);
