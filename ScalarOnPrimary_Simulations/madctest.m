clear;
close all;

n0 = [1, 1, 1, 1, 1];
n1 = [1.46535758, 1.46504135, 1.46472145, 1.4643964,  1.46406498];
n2 = [1.42443228, 1.42391618, 1.42338359, 1.42283329, 1.42226429];
n3 = [2.4456937,  2.44435182, 2.44317241, 2.44212632, 2.44119048];

dz = [-2.452e-9, -1.247e-9, 0, 1.305e-9, 2.683e-9];

lambdas = 1e-6*[2.0,2.1,2.2,2.3,2.4];

phi1 = 7.0516 * pi / 180;
phi2 = 3.8050 * pi / 180;
phi3 = 1.1465 * pi / 180;

clocking = -80 * pi/180;

for i = 1:5
    dz_out(i) = ADC_PRISM(dz(i), n0(i), n1(i), n2(i), n3(i), phi1, phi2, phi3, clocking);
end

OUT = dz_out - dz_out(3);

%% Gridscan optimization
samp = linspace(0,90,100);
tilt1 = linspace(0,90,100);
tilt2 = linspace(0,90,100);
tilt3 = linspace(0,90,100);
tilt4 = linspace(0,90,100);
tilt5 = linspace(0,90,100);

for i = 1:100
    
    clocking = deg2rad(samp(i));
%     disp(clocking);
    
    for j = 1:5
        tilt_out(j) = ADC_PRISM(dz(j), n0(j), n1(j), n2(j), n3(j), phi1, phi2, phi3, clocking);
    end
    
%     disp(tilt_out);
    
    tilt_out = tilt_out - tilt_out(3);
    %disp(tilt_out);
    
    tilt1(i) = tilt_out(1);
    tilt2(i) = tilt_out(2);
    tilt3(i) = tilt_out(3);
    tilt4(i) = tilt_out(4);
    tilt5(i) = tilt_out(5);
    
    tiltsums(i) = sum(tilt_out.^2);
end

figure()
hold on
title(['Gridsearch Min']);
xlabel(['Clocking Angle 0-90 deg']);
ylabel(['Sum of Squares']);
% plot(samp, tilt1,'Color', 'r');
% plot(samp, tilt2,'Color', 'y');
% plot(samp, tilt3,'Color', 'g');
% plot(samp, tilt4,'Color', 'c');
% plot(samp, tilt5,'Color', 'm');
plot(samp, tiltsums, 'Color', 'b');
hold off

[M,I] = min(tiltsums);
disp([tilt1(I),tilt2(I),tilt3(I),tilt4(I),tilt5(I)]);
disp(M);
disp(I);

figure()
hold on
title(['Polychromatic Tilt 2 - 2.4 um']);
xlabel(['Clocking Angle 0 - 90 deg']);
ylabel(['Tilt (radians)']);

plot(samp, tilt1,'Color', 'r');
plot(samp, tilt2,'Color', 'y');
plot(samp, tilt3,'Color', 'g');
plot(samp, tilt4,'Color', 'c');
plot(samp, tilt5,'Color', 'm');
% plot(samp, tiltsums, 'Color', 'b');
hold off

%% 
% ADC_PRISM(4, n0(1), n1(1), n2(1), n3(1), phi1, phi2, phi3, clocking);


