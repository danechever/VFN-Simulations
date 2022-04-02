
n1 = py.numpy.array((1.46535758, 1.46504135, 1.46472145, 1.4643964,  1.46406498))
n2 = py.numpy.array((1.42443228, 1.42391618, 1.42338359, 1.42283329, 1.42226429))
n3 = py.numpy.array((2.4456937,  2.44435182, 2.44317241, 2.44212632, 2.44119048))

dz = py.numpy.array((-2.452e-9, -1.247e-9, 0, 1.305e-9, 2.683e-9))

n0 = 1

'''
n1 = 1.4641
n2 = 1.4152
n3 = 2.4437
'''

phi1 = 7.0516 * py.numpy.pi / 180;
phi2 = 3.8050 * py.numpy.pi / 180;
phi3 = 1.1465 * py.numpy.pi / 180;

clocking = 45*numpy.pi/180;
tilt1 = 0;

#print(clocking)


dz_out = py.triple_prism.triple_prism(dz, n0, n1, n2, n3, phi1, phi2, phi3, clocking, tilt1);

disp(dz_out);





