import numpy as np
from snell_3d import snell_3d

def singleprism(dz, n0, n1, n2, n3, phi1, phi2, phi3, clocking=0, tilt1=0):
    """
    Model a triple prism ADC, assuming it is aligned with the direction of DAR
    
    Args:
        dz: incoming light displacement from normal of first surface/chief ray (radians)
        n0: index of refraction of air
        n1: index of refraction of first prism
        n2: index of refraction of second prism
        n3: index of refraction of third prism
        phi1: wedge angle of first prism (radians)
        phi2: wedge angle of second prism (radians)
        phi2: wedge angle of third prism (radians)
        clocking: a value between 0 and pi/2 to represent clocking angles. or an array of 6 clocking angles, all positive. (radians)
        tilt1: tilt of the first prism (in radians)
        
    Returns:
        dz_out: outgoing light displacement from normal of _first_ surface (radians)
        
    """
    if not hasattr(clocking, "__len__"):
        # turn this into a array of 6 elements
        clocking = np.ones(6) * clocking
    
    
    n0 = 1
        
    u0x = np.sin(dz)
    u0y = np.zeros(dz.shape)
    u0z = np.cos(dz)
    u0 = (u0x, u0y, u0z)
    
    # angle going through first prism (normal to chief ray when zenith angle=0)
    norm1 = (0, 0, 1)
    u1 = snell_3d(u0, norm1, n0, n1)
    #theta1 = np.arcsin(n0 / n1 * np.sin(dz)) - phi1
    #print("theta1", theta1 - np.median(theta1), theta1 + phi1)
    
    
    # angle going into 2nd prism
    phitot = -(-phi1)
    norm2 = (np.sin(phitot)*np.cos(clocking[0]), np.sin(phitot)*np.sin(clocking[0]), np.cos(phitot))
    u2 = snell_3d(u1, norm2, n1, n2)
    #theta2 = np.arcsin(n1 / n2 * np.sin(theta1)) + phi2
    #print("theta2", theta2 - np.median(theta2), "{0:.10f}".format(np.cos(theta2 + phi1 - phi2)[1]))
    
    # angle going into 3rd prism
    phitot = -(-phi1 + phi2)
    norm3 = (np.sin(phitot)*np.cos(clocking[1]), np.sin(phitot)*np.sin(clocking[1]), np.cos(phitot))
    u3 = snell_3d(u2, norm3, n2, n3)
    #theta3 = np.arcsin(n2 / n3 * np.sin(theta2)) + phi3
    #print("theta3", theta3 - np.median(theta3), "{0:.10f}".format(np.cos(theta3 + phi1 - phi2 - phi3)[1]))
    
    # angle going back into air (relative to normal of the 3/air prism boundary)
    phitot = -(-phi1 + phi2 + phi3)
    norm3p = (np.sin(phitot)*np.cos(clocking[2]), np.sin(phitot)*np.sin(clocking[2]), np.cos(phitot))
    u3p = snell_3d(u3, norm3p, n3, n0)
    #theta3p = np.arcsin(n3 / n0 * np.sin(theta3))
    #print(theta3p - np.median(theta3p), np.degrees(theta3p))
    
    dz_out = np.arctan2(u3p[0], u3p[2])
    #print(u6p)
    #print("debug", u4, u5, u6, u6p)
    
    return dz_out
