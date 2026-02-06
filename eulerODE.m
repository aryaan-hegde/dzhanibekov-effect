function wdot = eulerODE(~, w, I)
    Ix = I(1); Iy = I(2); Iz = I(3);
    wx = w(1); wy = w(2); wz = w(3);

    wdot = zeros(3, 1);
    wdotx = ((Iy - Iz)/Ix) * wy * wz; wdoty = ((Iz - Ix)/Iy) * wz * wx; wdotz = ((Ix - Iy)/Iz) * wx * wy;
    wdot(1) = wdotx; wdot(2) = wdoty; wdot(3) = wdotz;
end