function adcphz = simpleADC(inputs)

    mag = inputs.magfactor;
    dispersion = inputs.tilt_valsCopy;
    Y = inputs.coordsPP.Y;
    dx = inputs.coordsPP.dx;
    I = inputs.I; %inputs.gsmin;
    clocking = inputs.clocking;
    
    phi1 = 7.0516 * py.numpy.pi / 180;
    phi2 = 3.8050 * py.numpy.pi / 180;
    phi3 = 1.1465 * py.numpy.pi / 180;
    
    n0 = py.numpy.array(inputs.n0);
    n1 = py.numpy.array(inputs.n1);
    n2 = py.numpy.array(inputs.n2);
    n3 = py.numpy.array(inputs.n3);

    tilt1 = 0;
    
    dz_i = dispersion(:,I)';
    dz1 = mag * dz_i; %Multiply by the system magnification
    
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
    
    pOUT = pOUT/mag;
    pOUT = (pOUT - pOUT(ceil(inputs.numWavelengths/2)));
    
    for i = 1:inputs.numWavelengths
        
        adcphz(:,:,i) = 2*pi*(pOUT(i)/mag/inputs.lambdas(i))*Y*dx;
        
    end
end