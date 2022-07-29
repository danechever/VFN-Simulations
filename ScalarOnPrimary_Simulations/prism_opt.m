clear; close all;
addpath(['..' filesep 'VFNlib']);
load pyyo;
load oldetas;

% py.importlib.import_module('sellmeir1');
% py.importlib.import_module('triple_prism');
% sellmeir1coeffs;

%% Input parameters
sppvfnConfig;
load dzCopy;

%% Wedge and ADC Parameters

inpar.magfactor = 890.16;

%% Search for Solution

dz = dzCopy*inpar.magfactor;
% dz = py.tuple(dz);
% dz = py.numpy.array(dz);

for i = 1:inpar.numWavelengths
   
    n1(i) = getRefractiveIndex('mgf2',inpar.lambdas(i));
    
end

n1 = py.numpy.array(n1);
samp = linspace(0,90,100);
samp2 = linspace(-89,89,178);
sampall = [];
tilt_vals = [];
tilt_all = [];
tiltsums_all = [];

for j = 1:178
jv = samp2(j);
for i = 1:100
    
    alpha = deg2rad(samp(i));
%     j = deg2rad(j);
    
    n0 = 1;

    %%THROUGH PRISM
    %replace 0s with dz
    
%     u0x = py.numpy.sin([0,0,0,0,0]);
%     u0y = py.numpy.zeros(py.numpy.shape([0,0,0,0,0]));
%     u0z = py.numpy.cos([0,0,0,0,0]);
    
    u0x = py.numpy.sin([jv,jv,jv,jv,jv]);
    u0y = py.numpy.zeros(py.numpy.shape([jv,jv,jv,jv,jv]));
    u0z = py.numpy.cos([jv,jv,jv,jv,jv]);

%     u0x = py.numpy.sin(dz);
%     u0y = py.numpy.zeros(py.numpy.shape(dz));
%     u0z = py.numpy.cos(dz);
    
    u0 = {u0x,u0y,u0z};
    
    u0 = py.tuple(u0);
    norm1 = py.numpy.array([0, 0, 1]);
    
    u1 = py.snell_3d.snell_3d(u0,norm1,n0,n1);
    
    %%LEAVING PRISM
    normf = py.numpy.array([py.numpy.sin(alpha)*py.numpy.cos(0), py.numpy.sin(alpha)*py.numpy.sin(0), py.numpy.cos(0)]);
    u1f = py.snell_3d.snell_3d(u1,normf,n1,n0);

    u1f = u1f.tolist();
    tilt = py.numpy.arctan2(u1f(1),u1f(3));
    tilt_out = double(tilt);
    
    tilt_out = tilt_out - tilt_out(ceil(inpar.numWavelengths/2));
    tilt_out = tilt_out/inpar.magfactor;
    
    tilt_vals(:,i) = tilt_out'; 
    
    tiltsums(i) = sum((dz'/inpar.magfactor - tilt_out').^2); %sum(tilt_out'.^2);
    
end
tilt_all = horzcat(tilt_all,tilt_vals);
tiltsums_all = horzcat(tiltsums_all,tiltsums);
sampall = horzcat(sampall,samp);
end

dzfoc = dz/inpar.magfactor;
% slopes = tilt_vals(2,:) - tilt_vals(1,:);

%% For just prism optimization
% figure()
% hold on
% title(['Gridsearch Min']);
% xlabel(['Clocking Angle 0-90 deg']);
% ylabel(['Sum of Squares']);
% plot(samp, tiltsums, 'Color', 'b');
% hold off
% 
% [M,I] = min(tiltsums);
% disp(M);
% disp(I);
% disp(samp(I));
% disp(tilt_vals(:,I));
% 
% figure()
% hold on
% title(['Polychromatic Tilt 2 - 2.4 um']);
% xlabel(['Clocking Angle 0 - 90 deg']);
% ylabel(['Tilt (radians)']);
% plot(samp, tilt_vals,'Color','r');
% 
% for i = 1:inpar.numWavelengths
%     plot(samp, tilt_vals(i,:));
% end


%% For prism and angle optimization
figure()
hold on
title(['Gridsearch Min']);
xlabel(['Clocking Angle 0-90 deg']);
ylabel(['Sum of Squares']);
plot(sampall, tiltsums_all, 'Color', 'b');
hold off

[M,I] = min(tiltsums_all);
disp(M);
disp(I);
disp(sampall(I));
disp(tilt_all(:,I));

figure()
hold on
title(['Polychromatic Tilt 2 - 2.4 um']);
xlabel(['Clocking Angle 0 - 90 deg']);
ylabel(['Tilt (radians)']);
plot(sampall, tilt_all,'Color','r');

for i = 1:inpar.numWavelengths
    plot(sampall, tilt_all(i,:));
end
