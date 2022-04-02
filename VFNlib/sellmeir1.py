def sellmeir1(wvs, temp, K1, K2, K3, L1, L2, L3):
    """
    Compute index of refraction of a material given some coefficients
    
    Returns:
        n: index of refraction
    """
    import numpy as np;
    wvs2 = wvs**2
    
    K1_total = np.sum([K1[i] * temp**i for i in range(len(K1))])
    K2_total = np.sum([K2[i] * temp**i for i in range(len(K1))])
    K3_total = np.sum([K3[i] * temp**i for i in range(len(K1))])
    
    '''
    print(K1_total)
    print(K2_total)
    print(K3_total)
    '''
    L1_total = np.sum([L1[i] * temp**i for i in range(len(K1))])
    L2_total = np.sum([L2[i] * temp**i for i in range(len(K1))])
    L3_total = np.sum([L3[i] * temp**i for i in range(len(K1))])

    '''
    print(L1_total)
    print(L2_total)
    print(L3_total)
    '''
    #pr'''int(K1_total, K2_total, K3_total, L1_total, L2_total, L3_total)
    
    arg1 = K1_total*wvs2/(wvs2 - L1_total**2)
    arg2 = K2_total*wvs2/(wvs2 - L2_total**2)
    arg3 = K3_total*wvs2/(wvs2 - L3_total**2)

    '''
    print(arg1)
    print(arg2)
    print(arg3)
    '''
    
    n = np.sqrt(arg1 + arg2 + arg3 + 1)
    
    return n

'''
    ### a bunch of coefficients ###
    baf2_args = [3.91393303e-1, 7.61840241e-1, 4.13994379, 2.22978228e-4, 1.04999877e-2, 2.29871020e3]
    baf2_args = [[8.285359e-1, -8.986505e-04, -1.884197E-06, -1.332822E-10, 3.650068E-12],
             [3.315039E-01, 9.091254E-04, 1.656780E-06, 5.257707E-10, -3.904140E-12],
             [4.367314E+00, -1.161161E-02, 7.204123E-05, -4.302326E-08, -1.764139E-10],
             [ 8.362026E-02, 8.880306E-04, -1.277585E-05, 5.231437E-08, -7.312824E-11],
             [-1.148764E-01, 3.381142E-03, -1.897870E-05, 4.686248E-08, -4.348650E-11],
             [4.921549E+01, -6.672202E-02, 4.283633E-04, -3.280396E-07, -8.848551E-10]]

    caf2_args = [1.03032805, 4.03776479e-3, 4.09367374, 5.60918206e-3, 5.54020715e-2, 1.264110033e3]
    caf2_args = [[1.04834, -2.21666E-04, -6.73446E-06, 1.50138E-08, -2.77255E-11],
             [-3.32723E-03, 2.34683E-04, 6.55744E-06, -1.47028E-08, 2.75023E-11],
             [3.72693, 1.49844E-02, -1.47511E-04, 5.54293E-07, -7.17298E-10],
             [7.94375E-02, -2.20758E-04, 2.07862E-06, -9.60254E-09, 1.31401E-11 ],
             [0.258039, -2.12833E-03, 1.20393E-05, -3.06973E-08 , 2.79793E-11],
             [34.0169 , 6.26867E-02 , -6.14541E-04, 2.31517E-06,  -2.99638E-09]]

    znse_args = [4.29801490, 6.27765570e-1, 2.89556330, 3.6888190e-2, 1.43476258e-1, 2.20849196e3]
    znse_args = [[4.41367E+00, -1.13389E-03, 2.00829E-05, -8.77087E-08, 1.26557E-10],
             [4.47774E-01, 1.11709E-03, -1.80101E-05, 8.10837E-08, -1.18476E-10],
             [6.70952E+00, -8.18190E-02, 5.77330E-04, -1.89210E-06, 2.15956E-09],
             [1.98555E-01, -3.62359E-05, 7.20678E-07, -3.12380E-09, 4.51629E-12],
             [3.82382E-01, -1.56654E-04, 2.56481E-06, -1.07544E-08, 1.53230E-11],
             [7.33880E+01, -5.06215E-01, 3.06061E-03, -8.48293E-06, 6.53366E-09]]
'''


#wvs = np.array([1.3, 1.5, 2.4])
#temp = 273.15 + 3 #Kelvin

### test getting the index of refraction as a function of wavelength for each prism
#print(sellmeir1(wvs, temp, baf2_args[0], baf2_args[1], baf2_args[2], baf2_args[3], baf2_args[4], baf2_args[5]))
#print(sellmeir1(wvs, temp, caf2_args[0], caf2_args[1], caf2_args[2], caf2_args[3], caf2_args[4], caf2_args[5]))
#print(sellmeir1(wvs, temp, znse_args[0], znse_args[1], znse_args[2], znse_args[3], znse_args[4], znse_args[5]))