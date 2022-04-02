function angdif = getADCTiltOpt(dz, n0, n1, n2, n3, phi1, phi2, phi3, clock_angle)

    %tilt_ideal = [-2.235e-07, -2.223e-07, -2.210e-7, -2.197e-07, -2.184e-07];
    
%     tilt_ideal = [1.9895e-4, 1.9788e-4, 1.9677e-4, 1.9561e-4, 1.9438e-4];
%     
%     tilt_ideal = tilt_ideal - 1.9677e-4;
    
    clocking = [1,1,1,1,1,1]*clock_angle*pi/180;
    
    dz = dz.*890.16;
    
    dz_OUT = zeros(1,5);
    
    for ch = 1:5
        
        u0x = sin(dz(ch));
        u0y = 0;
        u0z = cos(dz(ch));
        %u0 = [u0x, u0y, u0z];
        u0t = py.tuple([u0x,u0y,u0z]);
        u0 = py.tuple({u0t,u0t});
    
        norm1t = py.tuple([0, 0, -1]);
        norm1 = py.tuple({norm1t,norm1t});
        
        u1 = py.snell_3d(u0, norm1, n0, n1(ch));

        % angle going into 2nd prism
        phitot = (phi1);
        norm2 = [(sin(phitot) * cos(clocking(1))), (sin(phitot) * sin(clocking(1))), -cos(phitot)];
        u2 = snell_in3d(u1, norm2, n1(ch), n2(ch));
    
        % angle going into 3rd prism
        phitot = (phi1 - phi2);    
        norm3 = [(sin(phitot) * cos(clocking(2))), (sin(phitot) * sin(clocking(2))), -cos(phitot)];
        u3 = snell_in3d(u2, norm3, n2(ch), n3(ch));
    
        % angle going into the air between triplets
        phitot = (phi1 - phi2 -phi3);
        norm3p = [(sin(phitot) * cos(clocking(3))), (sin(phitot) * sin(clocking(3))), -cos(phitot)];
        u3p = snell_in3d(u3, norm3p, n3(ch), n0);
    
        % angle going into 4th prism
        phitot = (phi3 + phi2 - phi1);
        norm4 = [(sin(phitot) * cos(clocking(4))), (-1*sin(phitot) * sin(clocking(4))), -cos(phitot)]; % clock other way now
        u4 = snell_in3d(u3p, norm4, n0, n3(ch));
    
        % angle going into 5th prism
        phitot = (phi2 - phi1);
        norm5 = [(sin(phitot) * cos(clocking(5))), (-1*sin(phitot) * sin(clocking(5))), -cos(phitot)]; % clock other way now
        u5 = snell_in3d(u4, norm5, n3(ch), n2(ch));
    
        % angle going through 6th prism
        phitot = (-phi1);
        norm6 = [(sin(phitot) * cos(clocking(6))), (-1*sin(phitot) * sin(clocking(6))), -cos(phitot)]; % clock other way now
        u6 = snell_in3d(u5, norm6, n2(ch), n1(ch));
    
        % angle leaving 6th prism relative to surface normal (which is also direction of chief ray)
        phitot = 0;
        norm6p = [0, 0, -1];
        u6p = snell_in3d(u6, norm6p, n1(ch), n0);
    
        dz_OUT(ch) = atan2(u6p(1), u6p(3))./890.16;
    
        %disp(dz_OUT)
    end
    dz_OUT = dz_OUT - dz_OUT(3);
    %dz_OUT = dz_OUT./890.16;
    
    %angdif = sqrt(sum(dz_OUT - tilt_ideal)^2);

    % angle going through first prism
    
    angdif = sqrt((0 - sum(dz_OUT)/5)^2);
    
end
