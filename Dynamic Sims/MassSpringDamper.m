clear; close all; clc;

% Parámetros
g    = 9.81; % m/s^2
m    = 1.0;  % kg
beta = 0.5;  % kg/s
k    = 2.0;  % N/m
F    = 0;    % N

% Parámetros del impulso

t_kick  = 10.0;   % Tiempo en el que sucede
dt_kick = 0.05;  % Duración del impácto
F_kick  = 10;    % Fuerza del impacto

J = F_kick * dt_kick; % Impulso [Ns]

% Condiciones iniciales
x0 = 0; % Posición inicial
v0 = 0; % Velocidad inicial

% Ecuación del movimiento (Tras la patada F = 0)
eom = @(t, x) [x(2); 1/m * (-beta*x(2)-k*x(1)+F)];

% Pre-Impacto
opts  = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

X0 = [x0; v0];
[T1, X1] = ode45(eom, [0, t_kick], X0, opts);

state_pre_kick = X1(end, :);
state_post_kick = [state_pre_kick(1); state_pre_kick(2) + (J/m)]; % Contiene el estado de valores iniciales luego del impulso

% Post-Impacto
tspan = [t_kick, 100];

[T2, X2] = ode45(eom, tspan, state_post_kick, opts);

T = [T1; T2];
pos = [X1(:, 1); X2(:, 1)];
vel = [X1(:, 2); X2(:, 2)];

% Campo vectorial Espacio Fase
[X, DX] = meshgrid(linspace(-5, 5, 22), linspace(-5, 5, 22));

U = DX;
V = 1/m * ((-beta*DX)-(k*X)+F);

mag = sqrt(U.^2 + V.^2) + 1e-6;

U_n = U ./ mag;
V_n = V ./ mag;

% Configuración de la figura
figure('Name', 'Diagrama de fase')
hold on

mag_flat = mag(:);
mag_min  = min(mag_flat);
mag_max  = max(mag_flat);

cmap = colormap(jet(256));

for i = 1:numel(X)
    c_idx = round(1 + 255 * (mag(i) - mag_min) / (mag_max - mag_min));
    col   = cmap(c_idx, :);

    quiver(X(i), DX(i), U_n(i) * 0.22, V_n(i) * 0.22, 0, ...
           'Color', col, 'LineWidth', 0.8, 'MaxHeadSize', 2)
end

cb = colorbar;
clim([mag_min, mag_max]);

% Gráficas
plot(pos, vel, LineStyle="-",LineWidth=2,Color='white');
plot(pos(1),   vel(1),   'wo', 'MarkerSize', 9, 'MarkerFaceColor', 'w', 'DisplayName', 'Start')
plot(pos(end), vel(end), 'w^', 'MarkerSize', 9, 'MarkerFaceColor', 'k', 'DisplayName', 'End (settled)')

xlim([-5, 5]); ylim([-5, 5]);
xlabel('Posición (m)');
ylabel('Velocidad (m/s)');
title('Diagrama de Fase del Sistema');
grid on;
hold off;

figure('Name','Respuesta temporal');
plot(T, pos, LineWidth=2);
grid on;
xlabel("Tiempo (s)");
ylabel("Posición (m)");
title("Respuesta temporal");
xline(t_kick, "--r", "Momento del impulso", "LabelVerticalAlignment","bottom");