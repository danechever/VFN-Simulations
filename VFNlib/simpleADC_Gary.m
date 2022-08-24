function alphay_rad = simpleADC_Gary(inputs)
    
    alphaIn_rad = [0,0,0,0,0];
    nAir = [1.000293,1.000293,1.000293,1.000293,1.000293];
    n1 = double(inputs.n1);
    n2 = double(inputs.n2);
    n3 = double(inputs.n3);
    nPrisms = [n1;n2;n3];
    
%     prismWedgeAngs_deg = [-7.0516,3.8050,1.1465];
    prismWedgeAngs_deg = [inputs.phi1_ADC, inputs.phi2_ADC, inputs.phi3_ADC];
    clockAng_deg = inputs.clocking;
    
    [alphax_rad,alphay_rad,~] = ADCrayTrace(alphaIn_rad,nAir,nPrisms,prismWedgeAngs_deg,clockAng_deg);
    
    alphax_rad = alphax_rad';
%     alphay_rad = alphay_rad';

end