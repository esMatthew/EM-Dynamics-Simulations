close all; clear; clc;

eps0 = 8.854e-12; % Vaccum permittivity
c = 299792458; % Speed of light [m/s]
q = 1e-6; % Charge [C]

T = 30;
dt = 0.1;

% Mesh and frame of reference creation
L = 10;
n = 200;

[x, y] = meshgrid(linspace(-L, L, n), linspace(-L, L, n));

R = sqrt(x.^2 + y.^2);

R(R < 0.2) = NaN;

sinth = abs(y) ./ R;

% Initial values for the particle at rest
v = 0;
beta = 0;
Emag = (q/(4*pi*eps0)) .* (1-beta)./((1-beta.^2.*sinth.^2).^(3/2) .* R.^2);
Ex = x .* Emag;
Ey = y .* Emag;

% Frame setup
figure;
axis equal;
axis([-L L -L L]);
xlabel('x');
ylabel('y');
title('Contours of relativistic electric field magnitude');
colormap(turbo);
colorbar;
hold on;

Eplot = log10(Emag);

[C, hcont] = contourf(x, y, Eplot, 30, 'LineColor', 'none');

hp = plot(0, 0, 'ro', 'MarkerFaceColor','r');

he = quiver(x, y, Ex, Ey);

% Animation
for t = 0:dt:T
    % Update the particle's position and velocity
    v = (t/T) * c * 0.99;
    beta = v / c;

    Emag = (q/(4*pi*eps0)) .* (1-beta)./((1-beta.^2.*sinth.^2).^(3/2) .* R.^2);

    Eplot = log10(Emag);
    % Update contours
    delete(hcont)
    [C, hcont] = contourf(x, y, Eplot, 30, 'LineColor','none');

    title(sprintf('Relativistic electric field magnitude (v = %.2f c)', beta));

    drawnow;
end