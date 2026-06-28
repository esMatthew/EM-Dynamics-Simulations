clear; clc; close all;

% Parameters
g   = 9.81;
L   = 1.0;
m   = 1.0;
mu  = 0.5;

% Initial Conditions
theta0    =  pi/2; % initial angle         (rad)
thetadot0 =  0;     % initial angular vel.  (rad/s)

% Equation of motion
eom = @(t, x) [ x(2) ;
               -(mu/m)*x(2) - (g/L)*cos(x(1)) ];

% Vector Field Grid
[TH, DTH] = meshgrid(linspace(-pi, 2*pi, 26), linspace(-10, 10, 22));

U = DTH;
V = -(mu/m)*DTH - (g/L)*cos(TH);

mag = sqrt(U.^2 + V.^2) + 1e-6;

% Normalize so all arrows are the same visual length
U_n = U ./ mag;
V_n = V ./ mag;

% Figure
figure('Name', 'Phase portrait — single trajectory')
hold on

mag_flat = mag(:);
mag_min  = min(mag_flat);
mag_max  = max(mag_flat);

cmap = colormap(jet(256));  

for i = 1:numel(TH)
    % Map this arrow's magnitude to a color index
    c_idx = round(1 + 255 * (mag(i) - mag_min) / (mag_max - mag_min));
    col   = cmap(c_idx, :);

    quiver(TH(i), DTH(i), U_n(i)*0.22, V_n(i)*0.22, 0, ...
           'Color', col, 'LineWidth', 0.8, 'MaxHeadSize', 2)
end

cb = colorbar;
clim([mag_min, mag_max]);

% Solve for a single Trajectory
x0    = [theta0; thetadot0];
tspan = [0, 100];
opts  = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

[T, X] = ode45(eom, tspan, x0, opts);

theta    = X(:,1);
thetadot = X(:,2);

% Draw the trajectory
plot(theta, thetadot, 'w-',  'LineWidth', 2.0)

% Mark start and end
plot(theta(1),   thetadot(1),   'wo', 'MarkerSize', 9, ...
     'MarkerFaceColor', 'w', 'DisplayName', 'Start')
plot(theta(end), thetadot(end), 'w^', 'MarkerSize', 9, ...
     'MarkerFaceColor', 'k', 'DisplayName', 'End (settled)')

% Formatting
xlabel('$\theta$  (rad)',        'FontSize', 13, 'Interpreter','latex')
ylabel('$\dot{\theta}$  (rad/s)', 'FontSize', 13, 'Interpreter','latex')
title('Phase portrait', 'FontSize', 12)

xticks([-pi, -pi/2, 0, pi/2, pi, 3*pi/2, 2*pi])
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi', '2*\pi/2', '3*\pi'})
xlim([-pi 2*pi]);  ylim([-10 10])
grid on;  box on
hold off

% Time Response

figure(Name='Time response');
plot(T, theta, LineWidth=2);
grid on;
xlabel("Time (s)");
ylabel("\theta (rad)");
title("Time response");