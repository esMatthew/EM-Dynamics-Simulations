close all; clear; clc;

function [Bx, By, Bz] = getField(I, a, x, y, z, nTh)
    % I: Wire current [A]
    % a: Wire Radius [m]
    % x, y, z: Observaytion point coordinates
    % nTh: number of divisions to be calculated

    mu = (4*pi) * 1e-7; % Permeability of free space

    th = linspace(0, 2*pi, nTh); % Theta step

    % Source points in loop
    rx0 = a*cos(th);
    ry0 = a*sin(th);
    rz0 = zeros(size(th));

    % Distance between field and source
    Rx = x - rx0;
    Ry = y - ry0;
    Rz = z - rz0;
    Rmag = sqrt(Rx.^2 + Ry.^2 + Rz.^2);

    % Field differentials calculations
    dBx = ((mu*I)/(4*pi)) .* ((a.*z.*cos(th))./(Rmag).^3);
    dBy = ((mu*I)/(4*pi)) .* ((a.*z.*sin(th))./(Rmag).^3);
    dBz = ((mu*I)/(4*pi)) .* ((a.^2 - a.*(x.*cos(th) + y.*sin(th)))./(Rmag).^3);

    % Numerically integrate the field potentials
    Bx = trapz(th, dBx);
    By = trapz(th, dBy);
    Bz = trapz(th, dBz);
end

nGrid = 10;
nTh = 100;
a = 5;
I = 1;

% Meshgrid
[X, Y, Z] = meshgrid(linspace(-10, 10, nGrid), linspace(-10, 10, nGrid), linspace(-5, 5, nGrid));

Bx = zeros(size(X));
By = zeros(size(Y));
Bz = zeros(size(Z));

for i = 1:nGrid
    for j = 1:nGrid
        for k = 1:nGrid
           [Bx(i, j, k), By(i, j, k), Bz(i, j, k)] = getField(I, a, X(i, j, k), Y(i, j, k), Z(i, j, k), nTh);
        end
    end
end

% Seed points for streamlines
nSeeds = 20;
theta_s = linspace(0, 2*pi, nSeeds);
z_seed = 0;

Xs = a * cos(theta_s);
Ys = a * sin(theta_s);
Zs = z_seed * ones(size(theta_s));

figure; hold on;
view(3);

quiver3(X, Y, Z, Bx, By, Bz, 'AutoScale','on', LineWidth=1.2, AutoScaleFactor=1);
grid on;

% Wire
th_plot = linspace(0, 2*pi, 400);
plot3(a*cos(th_plot), a*sin(th_plot), zeros(size(th_plot)), 'r', 'LineWidth', 3);

%streamline(X, Y, Z, Bx, By, Bz, Xs, Ys, Zs);

xlabel('x');
ylabel('y');
zlabel('z');
title('Magnetic Field from Current Loop');
axis equal;
view(3);