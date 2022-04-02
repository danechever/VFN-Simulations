clear;
close all;

load pyyo;
load oldetas;

lambda0 = 2.2e-6; %central wavelength
fracBW = 0.1818; %\Delta\lambda/\lambda
numWavelengths = 5;% number of discrete wavelengths 
lambdas = getWavelengthVec(lambda0,fracBW,numWavelengths);% array of wavelengths (meters)
keckD = 10.949; %Meters
lam0OverD = lambdas(ceil(numWavelengths / 2)) / keckD;
disp('Lam0OverD:');
disp(lam0OverD);


n0 = py.tuple([1, 1, 1, 1, 1]);
n1 = py.tuple([1.46535758, 1.46504135, 1.46472145, 1.4643964,  1.46406498]);
n2 = py.tuple([1.42443228, 1.42391618, 1.42338359, 1.42283329, 1.42226429]);
n3 = py.tuple([2.4456937,  2.44435182, 2.44317241, 2.44212632, 2.44119048]);

n0 = py.numpy.array(n0);
n1 = py.numpy.array(n1);
n2 = py.numpy.array(n2);
n3 = py.numpy.array(n3);

% disp([-2.452e-9, -1.247e-9, 0, 1.305e-9, 2.683e-9]);
disp('VFN Dispersion:');
disp(pyyCopy*lam0OverD);

% dz = [0 0 0 0 0];
dz_i = 1.0e-09 * [0.5475    0.1599         0    0.0493    0.2960];
disp('ADC input:');
disp(dz_i);


% dz = py.tuple([-2.452e-9, -1.247e-9, 0, 1.305e-9, 2.683e-9]*890.16);
dz = py.tuple(dz_i);
% dz = py.tuple(890.16*pyyCopy * lam0OverD);



dz = py.numpy.array(dz);

% '''
% n1 = 1.4641
% n2 = 1.4152
% n3 = 2.4437
% '''

phi1 = 7.0516 * py.numpy.pi / 180;
phi2 = 3.8050 * py.numpy.pi / 180;
phi3 = 1.1465 * py.numpy.pi / 180;

%clocking = -86.3938*pi/180;
% clocking = deg2rad(88.1818);
clocking = deg2rad(76.3636);
tilt1 = 0;

% #print(clocking)

dz_out = py.triple_prism.triple_prism(dz, n0, n1, n2, n3, phi1, phi2, phi3, clocking, tilt1);
dz_out = dz_out.tolist();
    
for i = 1:numWavelengths
    OUT(i) = dz_out{i};
end

% OUT

% OUT/890.16

OUT = (OUT - OUT(ceil(numWavelengths/2)))/890.16;

disp('ADC output Dispersion to Correct VFN:');
disp(OUT);
% OUT(2) - OUT(1)

%OUT(1)/OUT(2)





