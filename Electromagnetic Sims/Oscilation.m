close all; clear; clc;

% Parameters
eps0 = 8.854e-12; % Vaccum permittivity
c = 1; % Speed of light [m/s]
q = 1e-6; % Charge [C]

omega = 1; % Oscillation angular speed
A = 1;  % Oscillation Amplitude

T = 50; % Max time [s]
dt = 0.01; % Time step [s]

% Grid definition
n = 20;
L = 10;
[x, y] = meshgrid(linspace(-L, L, n), linspace(-L, L, n));

% Distance to particle (initial)
R = sqrt(x.^2 + y.^2);
R(R < 0.2) = NaN;

tr = 0;

% Plotting
figure;
axis equal;
axis([-L L -L L]);
xlabel('x');
ylabel('y');
title('Electromagnetic wave radiation caused by oscillating charge');

hold on;

Eabs = 0;
hImg = imagesc(linspace(-L,L,n), linspace(-L,L,n), Eabs);
clim([0 5e3])   % adjust upper limit as needed

set(gca,'YDir','normal');
colormap(turbo);
colorbar;

hE = quiver(x, y, zeros(size(R)), zeros(size(R)));
hE.LineWidth = 2;
hE.AutoScale = 'on';
hE.AutoScaleFactor = 1;

hQ = plot(0, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

fps = 60;
render_dt = 1/fps;
next_render_time = 0;

for t = 0:dt:T
    if t > next_render_time
        next_render_time = render_dt + next_render_time;

        tr = t - (R./c);

        % Update distance calculation
        Rx = x;
        Ry = y - (A.*sin(omega.*tr));
        R = sqrt(Rx.^2 + Ry.^2);
        R(R < 0.2) = NaN;

        Emag = -((A*q*omega^2)/(4*pi*eps0*c^2)) .* (sin(omega.*tr)./R.^3);
        Ex = Emag .* x .* (Ry);
        Ey = Emag .* (Ry.^2-R.^2);
    
        Eabs = sqrt(Ex.^2 + Ey.^2);
    
        yq = A * sin(omega * t);
    
        %set(hQ, 'XData', 0, 'YData', yq);
        set(hE, 'UData', Ex, 'VData', Ey);
        set(hImg, 'CData', Eabs);
    
        title(sprintf('Electromagnetic wave radiation caused by oscillating charge t = %.2f s', t));
        drawnow;
    end
end
