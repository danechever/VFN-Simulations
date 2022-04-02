import numpy as np
def snell_3d(uvec, nvec, n1, n2):
    """
    Snell's law in cartesian coordinates in 3D. We're going to assume uvec dot nvec is positive
    
    Args:
        uvec (tuple): 3 elements specifying ux, uy, uz of incident beam (unit vector). Shape (3, N)
        nvec (tuple): 3 elemnets specifying nx, ny, nz of surface normal (unit vector)
        n1: index of refraction of first element. Shape N. 
        n2: index of refraction of second element. Shape N. 
        
    Returns
        upvec (tuple): 3 elements specifying up_x, up_y, up_z of outgoing beam (unit vector)
    """
    '''
    uvec = (0,0,1)
    nvec = (0,0,1)
    '''
  
    uvec = np.array(uvec)
    nvec = np.array(nvec)
    
    eta = n1/n2
    udotn = np.sum(uvec*nvec[:,None], axis=0)
    upvec = (np.sqrt(1 - eta**2 + eta**2 * udotn**2) - udotn * eta)[None, :] * nvec[:, None]
    upvec = upvec + eta*uvec
    
    # make sure output is unit vector?
    upvec =  upvec / np.linalg.norm(upvec, axis=0)[None, :]
    
    return upvec
