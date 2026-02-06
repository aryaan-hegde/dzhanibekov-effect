clear; clc;
inertia_conversion = 1.82899785e-5;
epsilon = 0.01;
I = inertia_conversion * [64.726; 119.334; 70.944];
w0 = [epsilon; epsilon; 0.5];
tspan = [0, 300];

opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t_w, w] = ode45(@(t, w) eulerODE(t, w, I), tspan, w0);

w_interp = @(t) interp1(t_w, w, t)';

R0 = eye(3);
R0_vec = R0(:);

[t_R, Rvec] = ode45(@(t,R) rotationODE(t, R, w_interp), tspan, R0_vec);

N = length(t_R);
R = zeros(3,3,N);
for k = 1:N
    R(:,:,k) = reshape(Rvec(k,:), 3, 3);
end

e3_body = [0; 0; 1];

handle_lab = zeros(3, N);
for k = 1:N
    handle_lab(:,k) = R(:,:,k) * e3_body;
end

L = 1.0;
W = 0.4;

shaft_body = [ 0   0;
               0   0;
              -L/2 L/2 ];

cross_body = [ -W  W;
                0   0;
                L/2   L/2 ];

figure;
axis equal;
axis([-1.2 1.2 -1.2 1.2 -1.2 1.2]);
grid on;
view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Dzhanibekov Effect — T-Handle Animation');

hold on;

shaft_lab = R(:,:,1) * shaft_body;
cross_lab = R(:,:,1) * cross_body;

h_shaft = plot3(shaft_lab(1,:), shaft_lab(2,:), shaft_lab(3,:), ...
                'k', 'LineWidth', 3);
h_cross = plot3(cross_lab(1,:), cross_lab(2,:), cross_lab(3,:), ...
                'r', 'LineWidth', 3);

dt_sim = t_R(2) - t_R(1);
slow_factor = 2000;

for k = 1:5:length(t_R)

    shaft_lab = R(:,:,k) * shaft_body;
    cross_lab = R(:,:,k) * cross_body;

    set(h_shaft, 'XData', shaft_lab(1,:), ...
                 'YData', shaft_lab(2,:), ...
                 'ZData', shaft_lab(3,:));

    set(h_cross, 'XData', cross_lab(1,:), ...
                 'YData', cross_lab(2,:), ...
                 'ZData', cross_lab(3,:));

    drawnow;
    pause(slow_factor * dt_sim);
end