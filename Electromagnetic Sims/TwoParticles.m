clear;
clc;
close all;

% Parameters
Q1 = -1e-6;
Q2 = 1e-6;
k = 8.854187817e-12;
n = 25;
L = 2;

pos1 = [-L/2, 0, 0];
pos2 = [L/2, 0, 0];

% Grid definition
[X, Y, Z] = meshgrid(linspace(-L, L, n),linspace(-L, L, n),linspace(-L, L, n));

% Mathematical calculations

rx1 = X-pos1(1);
ry1 = Y-pos1(2);
rz1 = Z-pos1(3);
R1 = sqrt(rx1.^2 +ry1.^2 + rz1.^2);

E1 = (Q1*k)./ R1.^3;
E1(R1 == 0) = 0; % Avoid division by 0

Ex1 = E1 .* rx1;
Ey1 = E1 .* ry1;
Ez1 = E1 .* rz1;
                                                                    
rx2 = X-pos2(1);
ry2 = Y-pos2(2);
rz2 = Z-pos2(3);
R2 = sqrt(rx2.^2 + ry2.^2 + rz2.^2);

E2 = (Q2*k) ./ R2.^3;
E2(R2 == 0) = 0; % Avoid division by 0

Ex2 = E2 .* rx2;
Ey2 = E2 .* ry2;
Ez2 = E2 .* rz2;

Ex = Ex1 + Ex2;
Ey = Ey1 + Ey2;
Ez = Ez1 + Ez2;


% Field lines
numLines = 8;
theta = linspace(0, 2*pi, numLines);
phi   = linspace(0, pi, numLines);

[startX, startY, startZ] = sphere(numLines);
startX = startX(:)*0.3 + pos2(1);
startY = startY(:)*0.3 + pos2(2);
startZ = startZ(:)*0.3 + pos2(3);

% Visualization
figure; hold on;
quiver3(X, Y, Z, Ex, Ey, Ez, LineWidth=1.2);
streamline(X, Y, Z, Ex, Ey, Ez, startX, startY, startZ);

% Visualize the particle as a sphere
[particleX, particleY, particleZ] = sphere(20);
radius = 0.1; % Radius of the particle
surf(radius * particleX - pos1(1), radius * particleY - pos1(2), radius * particleZ - pos1(3), 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
surf(radius * particleX - pos2(1), radius * particleY - pos2(2), radius * particleZ - pos2(3), 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Set plot limits and labels
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Electric Field Visualization with Particle');
grid on;
view(3);
axis equal;
view(3);
hold off;