function dz_OUT = getADCTilt(dz, n0, n1, n2, n3, phi1, phi2, phi3, clock_angle)

    % angle going through first prism
    
    %tilt_ideal = 2.7630e-7;
    
    clocking = [1,1,1,1,1,1]*clock_angle*pi/180;
    
    dz = dz*890.16;
    
    u0x = sin(dz);
    u0y = 0;
    u0z = cos(dz);
    u0 = [u0x, u0y, u0z];
    
    norm1 = [0, 0, -1];
    u1 = snell_in3d(u0, norm1, n0, n1);

    % angle going into 2nd prism
    phitot = (phi1);
    norm2 = [(sin(phitot) * cos(clocking(1))), (sin(phitot) * sin(clocking(1))), -cos(phitot)];
    u2 = snell_in3d(u1, norm2, n1, n2);
    
    % angle going into 3rd prism
    phitot = (phi1 - phi2);    
    norm3 = [(sin(phitot) * cos(clocking(2))), (sin(phitot) * sin(clocking(2))), -cos(phitot)];
    u3 = snell_in3d(u2, norm3, n2, n3);
    
    % angle going into the air between triplets
    phitot = (phi1 - phi2 -phi3);
    norm3p = [(sin(phitot) * cos(clocking(3))), (sin(phitot) * sin(clocking(3))), -cos(phitot)];
    u3p = snell_in3d(u3, norm3p, n3, n0);
    
    % angle going into 4th prism
    phitot = (phi3 + phi2 - phi1);
    norm4 = [(sin(phitot) * cos(clocking(4))), (-1*sin(phitot) * sin(clocking(4))), -cos(phitot)]; % clock other way now
    u4 = snell_in3d(u3p, norm4, n0, n3);
    
    % angle going into 5th prism
    phitot = (phi2 - phi1);
    norm5 = [(sin(phitot) * cos(clocking(5))), (-1*sin(phitot) * sin(clocking(5))), -cos(phitot)]; % clock other way now
    u5 = snell_in3d(u4, norm5, n3, n2);
    
    % angle going through 6th prism
    phitot = (-phi1);
    norm6 = [(sin(phitot) * cos(clocking(6))), (-1*sin(phitot) * sin(clocking(6))), -cos(phitot)]; % clock other way now
    u6 = snell_in3d(u5, norm6, n2, n1);
    
    % angle leaving 6th prism relative to surface normal (which is also direction of chief ray)
    phitot = 0;
    norm6p = [0, 0, -1];
    u6p = snell_in3d(u6, norm6p, n1, n0);
    
    dz_OUT = atan2(u6p(1), u6p(3))/890.16;
    
end
    