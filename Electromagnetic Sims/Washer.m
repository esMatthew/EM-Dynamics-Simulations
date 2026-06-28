clear; clc; close all;

function [E_rho, E_z] = field_annulus(a1, a2, sigma, rho, th, z, Nr, Nth)
    % a1, a2 = radios arandela
    % sigma = carga superficial arandela (C/m^2)
    % rho, th, z = coordenadas cilíndricas pto de observación

    if nargin < 7
        Nr = 200; % Divisiones para discretizar la integral.
        Nth = 360;
    end

    % Valores constantes
    eps0 = 8.854187817e-12;
    k = 1 / (4*pi*eps0);

    rho_vec = linspace(a1, a2, Nr);
    th_vec = linspace(0, 2*pi, Nth).';

    % Coordenadas cilíndricas del orígen del campo
    [RHO_p, TH_p] = meshgrid(rho_vec, th_vec);

    % Distancia del origen del campo al punto de observación
    R = sqrt( rho.^2 + RHO_p.^2 - 2.*rho.*RHO_p.*cos(th - TH_p) + z.^2);

    % Factores comunes
    dA = RHO_p; % Factor jacobiano rho*drho*dtheta
    common = sigma * dA ./ (R.^3); % Factor común de ambas integrales

    % Cálculo de los integrandos para E_rho y E_z
    integrand_rho = common .* (rho - RHO_p .* cos(th - TH_p)); % Integrando radial
    integrand_z = common .* z; % Integrando axial

    % Resultados de las integrales angulares (con respecto a theta)
    int_th_rho = trapz(th_vec, integrand_rho, 1);   % 1×Nr result
    int_th_z = trapz(th_vec, integrand_z, 1);

    assert(length(rho_vec)==length(int_th_rho), 'rho_vec length mismatch');

    % resultados de las integrales radiales (con respecto a rho)
    int_rho = trapz(rho_vec, int_th_rho, 2); % Integral radial para E_rho
    int_z = trapz(rho_vec, int_th_z, 2); % Integral radial para E_z

    E_rho = k * int_rho; % Campo eléctrico radial
    E_z = k * int_z; % Campo eléctrico axial
end

% Parametros
a1 = 0.2;      % Radio interno
a2 = 0.5;      % Radio externo
sigma = 1e-6;  % Densidad de carga superficial C/m^2

% Malla para la visualización del campo
Nrho = 20;    Nth = 20;    Nz = 20; % Divisiones de rho, theta, y z
rho_vals = linspace(0, 2*a2, Nrho);
th_vals  = linspace(0, 2*pi, Nth);
z_vals   = linspace(-a2, a2, Nz);

% Inicialización de los elementos
X=[]; Y=[]; Z=[];
Ex=[]; Ey=[]; Ez=[];

% Cálculo del campo en cada punto del espacio
for i = 1:Nrho
  rho = rho_vals(i);
  for j = 1:Nth
    th = th_vals(j);
    for k = 1:Nz
      z = z_vals(k);
      [E_rho, E_zval] = field_annulus(a1, a2, sigma, rho, th, z); % Cálculo del campo

      % Conversión a coordenadas cartesianas
      x = rho * cos(th);
      y = rho * sin(th);

      X(end+1) = x;
      Y(end+1) = y;
      Z(end+1) = z;

      % Conversión de componentes a coordenadas cartesianas
      Ex(end+1) = E_rho * cos(th);
      Ey(end+1) = E_rho * sin(th);
      Ez(end+1) = E_zval;
    end
  end
end

% Convertir matrices a vectores de comlumna
X = X(:); Y = Y(:); Z = Z(:);
Ex = Ex(:); Ey = Ey(:); Ez = Ez(:);


% Malla regular (reshape)
Nx = 40; Ny = 40; Nz = 40;
[xg, yg, zg] = meshgrid(linspace(-a2*2, a2*2, Nx), linspace(-a2*2, a2*2, Ny), linspace(-a2*2, a2*2, Nz));

Exg = griddata(X, Y, Z, Ex, xg, yg, zg);
Eyg = griddata(X, Y, Z, Ey, xg, yg, zg);
Ezg = griddata(X, Y, Z, Ez, xg, yg, zg);

Nr = 5;
Nth = 10;
offset = 0.02 * a2;

rho_vals = linspace(a1, a2, Nr);
theta_vals = linspace(0, 2*pi, Nth);

xs_all = []; ys_all = []; zs_all = [];

for r = rho_vals
    for th = theta_vals
        x = r*cos(th);
        y = r*sin(th);

        xs_all = [xs_all, x, x];
        ys_all = [ys_all, y, y];
        zs_all = [zs_all, offset, -offset];
    end
end

streams = stream3(xg, yg, zg, Exg, Eyg, Ezg, xs_all, ys_all, zs_all);

% Gráfica del campo
figure;
hold on;
view(3);

streamline(streams);

Ntheta = 50;  Nr = 50; % Divisiones para graficar la arandela
[TH, R] = meshgrid(linspace(0, 2*pi, Ntheta), linspace(a1, a2, Nr));
XX = R .* cos(TH);
YY = R .* sin(TH);
ZZ = zeros(size(XX));
hW = surf(XX, YY, ZZ, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.4); % Visualización de la arandela

hQ = quiver3(X, Y, Z, Ex, Ey, Ez); % Visualización del campo

axis equal;
xlabel('x'); ylabel('y'); zlabel('z');
title('Electric field and charged washer geometry');
view(3);
grid on;
hold off;