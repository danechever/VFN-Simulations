function wphz = simplePRISM(inputs)

angles = linspace(-89,89,178);
mag = inputs.magfactor;
Y = inputs.coordsPP.Y;
dx = inputs.coordsPP.dx;
I = inputs.I; %inputs.gsmin;
alpha = inputs.clocking;%angles(inputs.clocking)
dispersion = inputs.tilt_valsCopy(:,I);
% clocking = inputs.clocking;

dz = dispersion' * mag;
    
for i = 1:inputs.numWavelengths
        
wphz(:,:,i) = 2*pi*(dz(i)/mag/inputs.lambdas(i))*Y*dx;
        
end
    

% %     j = deg2rad(j);
%     
%     n0 = 1;
% 
%     %%THROUGH PRISM
%     %replace 0s with dz
%     
% %     u0x = py.numpy.sin([0,0,0,0,0]);
% %     u0y = py.numpy.zeros(py.numpy.shape([0,0,0,0,0]));
% %     u0z = py.numpy.cos([0,0,0,0,0]);
%     
%     u0x = py.numpy.sin(dz);
%     u0y = py.numpy.zeros(py.numpy.shape(dz));
%     u0z = py.numpy.cos(dz);
% 
% %     u0x = py.numpy.sin(dz);
% %     u0y = py.numpy.zeros(py.numpy.shape(dz));
% %     u0z = py.numpy.cos(dz);
%     
%     u0 = {u0x,u0y,u0z};
%     
%     u0 = py.tuple(u0);
%     norm1 = py.numpy.array([0, 0, 1]);
%     
%     u1 = py.snell_3d.snell_3d(u0,norm1,n0,n1);
%     
%     %%LEAVING PRISM
%     normf = py.numpy.array([py.numpy.sin(alpha)*py.numpy.cos(0), py.numpy.sin(alpha)*py.numpy.sin(0), py.numpy.cos(0)]);
%     u1f = py.snell_3d.snell_3d(u1,normf,n1,n0);
% 
%     u1f = u1f.tolist();
%     tilt = py.numpy.arctan2(u1f(1),u1f(3));
%     tilt_out = double(tilt);
%     
%     tilt_out = tilt_out - tilt_out(ceil(inpar.numWavelengths/2));
%     tilt_out = tilt_out/mag;








end