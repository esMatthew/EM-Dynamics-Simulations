close all; clear; clc;

f = 1;
lambda = 1;
omega = 2*pi*f;
A = 1;
k = 2*pi/lambda;
v = lambda * f;
phi = 0;

t_max = 5;
dt = 0.05;
t = 0:dt:t_max;

x_max = 5;
dx = 0.03;
x = 0:dx:x_max;

[X, T] = meshgrid(x, t);

Y = A .* cos(k.*X - omega.*T);
Z = A .* cos(k.*X - omega.*T);

temp = zeros(size(X));

figure;
h1 = plot3(x, Y(1,:), temp, 'b', 'LineWidth', 1.5);
hold on;
h2 = plot3(x, temp, Z(1,:), 'r', 'LineWidth', 1.5);

grid on;
axis equal;
view(3);
xlabel('x'); ylabel('y'); zlabel('z');
title('Traveling electromagnetic-like wave');

for i = 1:length(t)
    set(h1, ...
        'XData', x, ...
        'YData', Y(i,:), ...
        'ZData', zeros(size(x)));

    set(h2, ...
        'XData', x, ...
        'YData', zeros(size(x)), ...
        'ZData', Z(i,:));

    title(sprintf('t = %.2f s', t(i)));
    drawnow;
end