function dz_out = triple_prism(dz, n0, n1, n2, n3, phi1, phi2, phi3, clocking, tilt1)
%     Model a triple prism ADC, assuming it is aligned with the direction of DAR
%     
%     Args:
%         dz: incoming light displacement from normal of first surface/chief ray (radians)
%         n0: index of refraction of air
%         n1: index of refraction of first prism
%         n2: index of refraction of second prism
%         n3: index of refraction of third prism
%         phi1: wedge angle of first prism (radians)
%         phi2: wedge angle of second prism (radians)
%         phi2: wedge angle of third prism (radians)
%         clocking: a value between 0 and pi/2 to represent clocking angles. or an array of 6 clocking angles, all positive. (radians)
%         tilt1: tilt of the first prism (in radians)
%         
%     Returns:
%         dz_out: outgoing light displacement from normal of _first_ surface (radians)

    if ~exist('clocking','var')
        clocking = 0;
    elseif ~exist('tilt1','var')
        tilt1 = 0;
    end
    
    c_size = size(clocking);
    
    if c_size(1) && c_size(2) == 1
        clocking = ones(1,6) * clocking;
    end

    n0 = 1;
    
    dz_size = size(dz);
    
    u0x = sin(dz);
    u0y = zeros(dz_size(1),dz_size(2));
    u0z = cos(dz);
    u0 = [u0x, u0y, u0z];
    
%     % angle going through first prism (normal to chief ray when zenith angle=0)
    norm1 = [0, 0, 1];
    u1 = snell_3d(u0, norm1, n0, n1);
%     %theta1 = np.arcsin(n0 / n1 * np.sin(dz)) - phi1
%     %print("theta1", theta1 - np.median(theta1), theta1 + phi1)
    
    
%     % angle going into 2nd prism
    phitot = -1*(-1*phi1);
    norm2 = [(sin(phitot) * cos(clocking(1))), (sin(phitot) * sin(clocking(1))), cos(phitot)];
    u2 = snell_3d(u1, norm2, n1, n2);
%     %theta2 = np.arcsin(n1 / n2 * np.sin(theta1)) + phi2
%     %print("theta2", theta2 - np.median(theta2), "{0:.10f}".format(np.cos(theta2 + phi1 - phi2)[1]))
    
%     % angle going into 3rd prism
    phitot = -1*(-1*phi1 + phi2);
        
    norm3 = [sin(phitot) * cos(clocking(2)), (sin(phitot) * sin(clocking(1))), cos(phitot)];
    u3 = snell_3d(u2, norm3, n2, n3);
%     %theta3 = np.arcsin(n2 / n3 * np.sin(theta2)) + phi3
%     %print("theta3", theta3 - np.median(theta3), "{0:.10f}".format(np.cos(theta3 + phi1 - phi2 - phi3)[1]))
    
    % angle going back into air (relative to normal of the 3/air prism boundary)
    phitot = -1*(-1*phi1 + phi2 + phi3);
    norm3p = [(sin(phitot) * cos(clocking(3))), (sin(phitot) * sin(clocking(3))), cos(phitot)];
    
    %norm3p = (np.sin(phitot)*np.cos(clocking[2]), np.sin(phitot)*np.sin(clocking[2]), np.cos(phitot));
    u3p = snell_3d(u3, norm3p, n3, n0);
    %theta3p = np.arcsin(n3 / n0 * np.sin(theta3))
    %print(theta3p - np.median(theta3p), np.degrees(theta3p))
    
    %print(u3p)
    %print("debug", u1, u2, u3, u3p)
    
    % angle going into 4th prism
    phitot = (phi3 + phi2 - phi1);
    norm4 = [(sin(phitot) * cos(clocking(4))), (-1*sin(phitot) * sin(clocking(4))), cos(phitot)]; % clock other way now
    u4 = snell_3d(u3p, norm4, n0, n3);
    %theta4 = np.arcsin(n0 / n3 * np.sin(theta3pp)) + phi3
    %print("theta4", theta4 - np.median(theta4), np.degrees(theta4  +phi2 -phi1))
    
    % angle going into 5th prism
    phitot = (phi2 - phi1);
    norm5 = [(sin(phitot) * cos(clocking(5))), (-1*sin(phitot) * sin(clocking(5))), cos(phitot)]; % clock other way now
    u5 = snell_3d(u4, norm5, n3, n2);
    %theta5 = np.arcsin(n3 / n2 * np.sin(theta4)) + phi2
    %print("theta5", theta5 - np.median(theta5), np.degrees(theta5 - phi1))
    
    % angle going through 6th prism
    phitot = (-phi1);
    norm6 = [(sin(phitot) * cos(clocking(6))), (-1*sin(phitot) * sin(clocking(6))), cos(phitot)]; % clock other way now
    u6 = snell_3d(u5, norm6, n2, n1);
    %print("theta6", theta6 - np.median(theta6), np.degrees(theta6))
    
    % angle leaving 6th prism relative to surface normal (which is also direction of chief ray)
    phitot = 0;
    norm6p = [0, 0, 1];
    u6p = snell_3d(u6, norm6p, n1, n0);
    %theta6p = np.arcsin(n1 / n0 * np.sin(theta6))
    %print("theta6p", theta6p - np.median(theta6p), np.degrees(theta6p))
    %print("dz_out", theta6p)
    
    dz_out = atan2(u6p(1,:), u6p(3,:));
    %print(u6p)
    %print("debug", u4, u5, u6, u6p)
end

