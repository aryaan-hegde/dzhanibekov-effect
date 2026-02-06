function Rdot_vec = rotationODE(t, Rvec, w_interp)
    R = reshape(Rvec, 3, 3);

    % Angular velocity at time t
    w = w_interp(t)';   % column vector [wx; wy; wz]

    % Skew-symmetric matrix
    Omega = [  0   -w(3)  w(2);
             w(3)   0   -w(1);
            -w(2)  w(1)   0 ];

    % Matrix ODE
    Rdot = R * Omega;

    % Return as vector
    Rdot_vec = Rdot(:);
end
