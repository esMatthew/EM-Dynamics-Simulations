clear; clc; close all;

% Parameters
Q = 1e-6; % C
epsilon0 = 8.85e-12; % F/m
mu0 = (4*pi) * 10^(-7); % H/m
k = 1/(4*pi*epsilon0);
n = 8;

v0 = 0.1; % m/s
vel = [0, 0, v0];

% Time parameters
Tmax = 100; % s
dt = 0.01;

% Space grid
[X, Y, Z] = meshgrid(linspace(-10, 10, n), linspace(-10, 10, n), linspace(-10, 10, n));

figure; hold on;

hE = quiver3(X, Y, Z, zeros(size(X)), zeros(size(Y)), zeros(size(Z))); % Initialize the quiver view for E
hE.LineWidth = 2;

hB = quiver3(X, Y, Z, zeros(size(X)), zeros(size(Y)), zeros(size(Z))); % Initialize the quiver view for B
hB.LineWidth = 2;

[spx, spy, spz] = sphere(20);
radius = 1;

% Initialize the surface view for the sphere representing the particle
hQ = surf(radius*spx, radius*spy, radius*spz, 'FaceColor','r', 'EdgeColor','none', 'FaceAlpha',0.8);

axis equal; grid on; view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');
title("Movement of an electromagnetic field");

xlim([-10 10])
ylim([-10 10])
zlim([-10 10])

axis manual

fps = 60;
render_dt = 1/fps;
next_render_time = 0;

for t = 0:dt:Tmax
    % Update the particle's position
    pos = [0, 0, v0*t - 5];

    % Calculate the position components
    Rx = X - pos(1);
    Ry = Y - pos(2);
    Rz = Z - pos(3);

    R = sqrt(Rx.^2 + Ry.^2 + Rz.^2);

    R(R == 0) = eps; % Avoid division by 0

    % Calculate the electric field components
    E = (k*Q)./(R.^3);
    Ex = E .* Rx;
    Ey = E .* Ry;
    Ez = E .* Rz;

    % Calculate the magnetic field components
    B = (mu0/(4*pi)) .* (Q./(R.^3));
    Bx = -B .* v0 .* Ry;
    By = B .* v0 .* Rx;
    Bz = zeros(size(Z));

    % Update the quiver plots for electric and magnetic fields
    if t >= next_render_time
        next_render_time = render_dt + next_render_time;
        set(hQ, 'XData', radius*spx + pos(1), 'YData', radius*spy + pos(2), 'ZData', radius*spz + pos(3));
        set(hE, 'UData', Ex, 'VData', Ey, 'WData', Ez);
        set(hB, 'UData', Bx, 'VData', By, 'WData', Bz);
        drawnow limitrate
    end
end