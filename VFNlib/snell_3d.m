function upvec = snell_3d(uvec, nvec, n1, n2)

%     Snell's law in cartesian coordinates in 3D. We're going to assume uvec dot nvec is positive
%     
%     Args:
%         uvec: 3 rows specifying ux, uy, uz of incident beam (unit vector). Shape (3, N)
%         nvec: 3 rows specifying nx, ny, nz of surface normal (unit vector)
%         n1: index of refraction of first element. Shape N. 
%         n2: index of refraction of second element. Shape N. 
%         
%     Returns
%         upvec: 3 elements specifying up_x, up_y, up_z of outgoing beam (unit vector)

    eta = n1/n2;
    
    udotn = sum(uvec.*nvec', 1);
    
    upvec = (sqrt(1 - eta^2 + eta^2 * udotn.^2) - udotn * eta).*nvec';
    upvec = upvec + eta * uvec;
    
    upvec = upvec./sqrt(sum(upvec.^2,1));
end