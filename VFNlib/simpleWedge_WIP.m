function wphz = simpleWedge(wedge_angle, n_wedge, lambda, coords)


wphz = (2*pi/lambda) * (n_wedge - 1) * tan(wedge_angle) * coords.X*coords.dx;

end