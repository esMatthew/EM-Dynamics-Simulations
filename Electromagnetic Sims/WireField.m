clear; close all; clc;

% Parameters

k = 8.854187817e-12;
L = 5;
n = 10;
lambda = 1; % C/m

[X, Y, Z] = meshgrid(linspace(-5, 5, n), linspace(-5, 5, n), linspace(-L, L, n));

rho = sqrt(X.^2 + Y.^2);

u1 = Z + L/2;
u2 = Z - L/2;

E_rho = (k*lambda) .* ( u1./sqrt((u1).^2 + (rho).^2) - u2./sqrt((u2).^2 + (rho).^2) ) ./ rho;
Ez = ((lambda*k)/3) * ( 1 ./ sqrt(rho.^2 + u2.^2) - 1 ./ sqrt(rho.^2 + u1.^2) );

% Avoid division by 0 in Ez
Ez(rho == 0) = eps; % Set a small value where rho is zero
E_rho(rho == 0) = eps;

Ex = E_rho .* X./rho;
Ey = E_rho .* Y./rho;

% Avoid division by 0 in Ex and Ey
Ex(rho == 0) = 0; % Set Ex to 0 where rho is zero
Ey(rho == 0) = 0; % Set Ey to 0 where rho is zero

% Calculate the electric field magnitude
E_magnitude = sqrt(Ex.^2 + Ey.^2 + Ez.^2);

figure; hold on;
quiver3(X, Y, Z, Ex, Ey, Ez, "LineWidth", 2);
z_wire = linspace(-L/2, +L/2, 100);
plot3( zeros(size(z_wire)), zeros(size(z_wire)), z_wire, 'r', 'LineWidth', 5);

xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Electric Field Visualization with a rod');
grid on;
view(3)
axis equal;
view(3);
hold off;