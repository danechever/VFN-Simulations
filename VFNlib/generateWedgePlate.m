function OUT = generateWedgePlate(inputs, wedge, wavelength)

wedge_mat = inputs.wedge_mat;
lambda0 = inputs.lambda0;
lam0OverD_rad = inputs.lam0OverD_rad;
lambdaOverD = inputs.lambdaOverD;
Ycoords = inputs.Ycoords;
N = inputs.N;
p_val = inputs.p;

OUT = 2*pi*(getWedgeTilt(wedge_mat, wedge, 1e6*wavelength)...
        /lam0OverD_rad*lambda0/wavelength)*lambdaOverD...
        *Ycoords/N - 2*pi*p_val*lambdaOverD*Ycoords/N;