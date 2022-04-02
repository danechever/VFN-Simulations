function OUT = optClock(dz, n0, n1, n2, n3, phi1, phi2, phi3, clock_angle,tilt1)
    
    tilt = py.triple_prism.triple_prism(dz,n0,n1,n2,n3,phi1,phi2,phi3,clock_angle,tilt1);
    
%     tilt = py.triple_prism.triple_prism(dz,n0,n1,n2,n3,phi1,phi2,phi3,-86.3938,tilt1);
    tilt = tilt.tolist();
    
    for i = 1:5
        OUT(i) = tilt{i};
    end
    
    OUT = OUT - OUT(3);
    
    %OUT = sum(OUT.^2 - sum(OUT)/5);
    %OUT = sqrt(sum(OUT)/5);
    %OUT = sum(abs(OUT));
    
    OUT = sum(OUT.^2);
end