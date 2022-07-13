function wphz = simpleWedge_max(inputs, wavelength)

    wphz = 2*pi*(1/inputs.magfactor*(getRefractiveIndex(inputs.wedge_mat, wavelength) - 1) * tan(inputs.wedge_angle)) * inputs.yvals;

end