function [tilt_valsCopy, clocking, I] = simple_PRISMOPT(inpar)

inpar.lam0OverD = inpar.lambdas(ceil(inpar.numWavelengths / 2)) / inpar.keckD;

dz = inpar.pyyCopy*inpar.lam0OverD;
inpar.magfactor = 890.16;

dz = dz*inpar.magfactor;

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


[M,I] = min(tiltsums_all);
disp(M);
disp(I);
disp(sampall(I));
disp(tilt_all(:,I));

tilt_valsCopy = tilt_all;

end
