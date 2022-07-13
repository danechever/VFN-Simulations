function OUT = generateWedgePlate(inputs, wedge, wavelength, wedgetype)

wedge_mat = inputs.wedge_mat;
lambda0 = inputs.lambda0;
lam0OverD_meters = inputs.lam0OverD_meters;
lambdaOverD_pixels = inputs.lambdaOverD; %Use this for OLD FFT based equation
Ycoords = inputs.yvals; %inputs.yvalsPP; %yvals
N = inputs.N;
p_val = inputs.p;
clocking = inputs.clocking;
tilt_valsCopy = inputs.tilt_valsCopy;
I = 85;

if(wedgetype == 'wedge')
    %% BEAM deviation INPUT TO WEDGE (radians)
    
%    disp(getWedgeTilt(wedge_mat, wedge, 1e6*wavelength)/890.16);

%     disp(getWedgeTilt(wedge_mat, wedge, 1e6*wavelength)); %Snell's Law Shift

%    disp(2*pi*(getWedgeTilt(wedge_mat, wedge, 1e6*wavelength)/890.16/lam0OverD_rad*lambda0/wavelength)*lambdaOverD);
    
    %% OUTPUT Phase
%     OUT = 2*pi*(1/890.16*getWedgeTilt(wedge_mat, wedge, 1e6*wavelength)/(wavelength/inputs.keckD)*(wavelength/lambda0))...
%         *lambdaOverD*Ycoords/N - 2*pi*p_val*lambdaOverD*Ycoords/N;
%     %%fix this
%     
%     OUT = 2*pi*(1/890.16*getWedgeTilt(wedge_mat, wedge, 1e6*wavelength)/lam0OverD_rad*lambda0/wavelength)*lambdaOverD*Ycoords/N;% - 2*pi*p_val*lambdaOverD*Ycoords/N;
%**************************************************************************       
%DFT Implementation
%     OUT = 2*pi*(1/890.16*getWedgeTilt(wedge_mat, wedge, 1e6*wavelength)/lam0OverD_meters*(lambda0/wavelength))...
%         *lambdaOverD_pixels*Ycoords/N - 2*pi*p_val*lambdaOverD_pixels*Ycoords/N;
%**************************************************************************    
%     ____
% lambda0OverD = lambda0Fnum_meters/foc;
% tiltFP = [-0.5 0 0.5] * lambda0OverD;
% wphz(:,:,i) = 2*pi*tiltFP(i)/lambda0OverD*lambda0/lambdas(i)*coordsPP.Y/N;
% sphz(:,:,i) = phz - wphz(:,:,i);
% Epup_Shift(:,:,i) = exp(1i*sphz(:,:,i)).*PUPIL;
%     ____
%**************************************************************************    
%MFT Implementation
    OUT = 2*pi*(1/890.16*getWedgeTilt(wedge_mat,wedge,1e6*wavelength)/lam0OverD_meters*lambda0/wavelength)*Ycoords/N; % - 2*pi*p_val*lambdaOverD_pixels*Ycoords/N; % Use this for the DM
%**************************************************************************    
%     disp(inputs.pyyo);
%     fprintf(['Adjusted wedgeTilt: ']);
%     disp(getWedgeTilt(wedge_mat, wedge, 1e6*wavelength)/890.16 - 2*pi*p_val*lambdaOverD);
        
elseif(wedgetype == 'ADC__')
    
    %% Input Parameters
    phi1 = 7.0516 * py.numpy.pi / 180;
    phi2 = 3.8050 * py.numpy.pi / 180;
    phi3 = 1.1465 * py.numpy.pi / 180;
    
    n0 = py.numpy.array(inputs.n0);
    n1 = py.numpy.array(inputs.n1);
    n2 = py.numpy.array(inputs.n2);
    n3 = py.numpy.array(inputs.n3);

%     clocking = -88.1818*pi/180;
    tilt1 = 0;
    
%     wedge_mat = 'CaF2';
%     wedge = wedge_angle;
%     inputs.lambdas = inPar.lambdas;
    %% Beam deviation INPUT to ADC (radians)
%     for ch = 1:inputs.numWavelengths
%         dz(ch) = (-1/890.16)*(getWedgeTilt(wedge_mat, wedge, 1e6*inputs.lambdas(ch)) - getWedgeTilt(wedge_mat, wedge, 1e6*inputs.lambdas(3)));
%     end
    
%     dz = inputs.pyyCopy * inputs.lam0OverD_meters;
%     
%     fprintf(['VFN Dispersion on Sky: ']);
%     disp(dz);
%     
%     vfn_disp_mas = (3600000*180/pi)*(dz(end) - dz(1));
% 
%     fprintf(['VFN Dispersion on sky across band (mas): ']);
%     disp(vfn_disp_mas);
%     
%     fprintf(['VFN Dispersion Magnified in Pupil Plane: ']);
%     disp(dz*890.16);
%     
%     dz_i = 1.0e-09 * [0.5475    0.1599         0    0.0493    0.2960]; %[0.5477    0.1601         0    0.0492    0.2967]; %Why is this hardcoded???
% 
%     fprintf(['ADC input on Sky: ']);
%     disp(dz_i);
%     
%     dz1 = dz_i * 890.16;
%     
%     fprintf(['ADC input Magnified in Pupil Plane: ']);
%     disp(dz1);
%     dz_dev = dz;
% %     dz = 1*10^-5*[-0.5290 -0.2481 0 0.2214 0.4209]; %[0 0 0 0 0 ];
% 
%     % Convert to Python
%     dz_p = py.tuple([dz1]);
%     dz_p = py.numpy.array(dz1);
    
    % #print(clocking)
    
    %Here we use the output dispersion from ADC_OPT to feed the ADC
    %function.
    %This "reverses" the direction of the beam, approximating the best-case
    %null shift 
    
%     dz_i = 1e-8 * [0.2193 0.0407 0 0.0717 0.2358]; %"minimal" offset determined by ADC_OPT
    dz_i = tilt_valsCopy(:,I)';
    dz1 = 890.16 * dz_i; %Multiply by the system magnification
    
    %Convert to Python
    dz_p = py.tuple([dz1]);
    dz_p = py.numpy.array(dz1);

    %% ADC OUTPUT beam deviation (radians)
    dz_out = py.triple_prism.triple_prism(dz_p, n0, n1, n2, n3, phi1, phi2, phi3, clocking, tilt1);
    dz_out = dz_out.tolist();
    
    % Convert to MATLAB
    for i = 1:inputs.numWavelengths
        pOUT(i) = dz_out{i};
    end
    
    pOUT = pOUT/890.16;
    pOUT = (pOUT - pOUT(ceil(inputs.numWavelengths/2)));
    
%     fprintf(['ADC Output Dispersion Magnified in Pupil Plane: ']);
%     disp(pOUT);
    
%     fprintf(['ADC output Dispersion on Sky: ']);
% 
%     disp(pOUT/890.16);
    
%     t_dif = pOUT/890.16 + dz;
%     t_dif_mas = (3600000*180/pi)*(t_dif(end) - t_dif(1));
%     
%     fprintf(['Peak to Valley Dispersion (on sky): ' num2str(t_dif_mas) '\n'] );
%     
%     disp(pOUT);
%     disp(dz_dev + pOUT);
    %% Compute OUTPUT Phase
%     dOUT = [];

%     2*pi*(1/890.16*getWedgeTilt(wedge_mat,wedge,1e6*wavelength)/lam0OverD_meters*lambda0/wavelength)*Ycoords/N;

%**************************************************************************    
% FFT Implementation
%     for i = 1:35
%         fOUT(:,:,i) = 2*pi*(pOUT(i)/(890.16)/lam0OverD_meters*lambda0/inputs.lambdas(i))*lambdaOverD_pixels*Ycoords/N; %- 2*pi*p_val*lambdaOverD_pixels*Ycoords/N;
        %disp(2*pi*(pOUT(i)/890.16/lam0OverD_rad*lambda0/inputs.lambdas(i))*lambdaOverD);
%         fOUT(:,:,i) = 2*pi*(pOUT(i)/lam0OverD_meters*(lambda0/inputs.lambdas(i)))...
%         *lambdaOverD_pixels*Ycoords/N;% - 2*pi*p_val*lambdaOverD_pixels*Ycoords/N;
%     end
%**************************************************************************    
% MFT Implementation
%     disp(inputs.numWavelengths);
    for i = 1:inputs.numWavelengths
        fOUT(:,:,i) = 2*pi*(pOUT(i)/(890.16)/lam0OverD_meters*lambda0/inputs.lambdas(i))*Ycoords/N; % - 2*pi*p_val*lambdaOverD_pixels*Ycoords/N;
        %disp(2*pi*(pOUT(i)/890.16/lam0OverD_rad*lambda0/inputs.lambdas(i))*lambdaOverD);
    end
%**************************************************************************    
    
    OUT = fOUT;
    
end

end
    