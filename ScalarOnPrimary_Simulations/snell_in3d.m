function [OUT,term2,term1] = snell_in3d(uv, nv, n1, n2)
    
    eta = n1/n2;
    
    uv_normal = sqrt(sum(uv.^2));
    nv_normal = sqrt(sum(nv.^2));
    
    
    term2 = nv/nv_normal*sqrt(1 - eta^2 * cross(nv/nv_normal,uv/uv_normal)*cross(nv/nv_normal,uv/uv_normal)');
    
    term1 = eta * (cross(nv/nv_normal,cross(-1*nv/nv_normal, uv/uv_normal)));
    
    OUT = term1 - term2;
    
end
