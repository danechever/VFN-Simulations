function OUT = ADC_PRISM(dz, n0, n1, n2, n3, phi1, phi2, phi3, clocking)
    % angle going through first prism
    
    %u0 = dz;
    
    u0 = [sin(dz), 0, cos(dz)];
    
    norm1 = [0, 0, -1];
    u1 = snell_in3d(u0, norm1, n0, n1);
%     disp(u1);
    
    % angle going into 2nd prism
    phitot = phi1;
    norm2 = [(sin(phitot) * cos(clocking)), (sin(phitot) * sin(clocking)), cos(phitot)];
%     disp(norm2);
    u2 = snell_in3d(u1, norm2, n1, n2);
%     disp(u2);
    
    % angle going into 3rd prism
    phitot = -1*(-1*phi1 + phi2);    
    norm3 = [sin(phitot) * cos(clocking), (sin(phitot) * sin(clocking)), cos(phitot)];
    u3 = snell_in3d(u2, norm3, n2, n3);
    
    % angle going into the air between triplets
    phitot = -1*(-1*phi1 + phi2 + phi3);
    norm3p = [(sin(phitot) * cos(clocking)), (sin(phitot) * sin(clocking)), cos(phitot)];
    u3p = snell_in3d(u3, norm3p, n3, n0);
    
    % angle going into 4th prism
    phitot = (phi3 + phi2 - phi1);
    norm4 = [(sin(phitot) * cos(clocking)), (-1*sin(phitot) * sin(clocking)), cos(phitot)]; % clock other way now
    u4 = snell_in3d(u3p, norm4, n0, n3);
    
    % angle going into 5th prism
    phitot = (phi2 - phi1);
    norm5 = [(sin(phitot) * cos(clocking)), (-1*sin(phitot) * sin(clocking)), cos(phitot)]; % clock other way now
    u5 = snell_in3d(u4, norm5, n3, n2);
    
    % angle going through 6th prism
    phitot = (-phi1);
    norm6 = [(sin(phitot) * cos(clocking)), (-1*sin(phitot) * sin(clocking)), cos(phitot)]; % clock other way now
    u6 = snell_in3d(u5, norm6, n2, n1);
    
    % angle leaving 6th prism relative to surface normal (which is also direction of chief ray)
    phitot = 0;
    norm6p = [0, 0, 1];
    u6p = snell_in3d(u6, norm6p, n1, n0);
    %disp(u6p);
    
    OUT = atan2(u6p(1), u6p(3));
end
