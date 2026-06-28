%% Gear-Link Animation: Problem 16-128
% Animates a gear rolling on a fixed rack, with link AB connecting
% a pin on the gear to a slider constrained to a horizontal rail.
% Live velocity and acceleration vectors are drawn at A and B.
%
% Geometry (all in inches):
%   Gear outer radius  R = 3 in  (rolls on rack, center O rises/falls? No —
%   the gear rolls on a FIXED rack so O stays at constant height y=3)
%   Inner pin radius   r = 2 in  (pin A is 2 in below O, i.e. at gear contact)
%   Link length        L = 8 in
%   Rail height        y_rail = 3 + 3 = 6 in  (O is at 3, top of gear + some)
%
% Note: In the problem snapshot, at the shown instant:
%   - A is directly below O at distance 2 in (bottom of inner hub)
%   - B is on the horizontal rail at height = O_y + 3*sin60 ... 
%     Actually from the figure B is on a rail at the top.
%     We compute B position kinematically each frame.

clear; clc; close all;

%% --- Parameters ---
R      = 3;     % gear outer radius (in)
r_pin  = 2;     % pin A radial distance from O (in)
L      = 8;     % link AB length (in)
omega  = 6;     % gear angular velocity (rad/s), CCW positive
alpha  = 12;    % gear angular acceleration (rad/s^2)

% Gear center O stays at fixed height (rolling on rack below)
O_y = R;        % = 3 in above rack

% Rail height for slider B (from figure, rail is above the gear)
% B is constrained to move horizontally. We set rail at y = O_y + 3 in
% (the pin at top of outer gear). From the figure the slider is at the
% same height as the top attachment. We'll compute y_B from geometry.
% The link goes from A (on the gear at radius r_pin from O) to B on rail.
% In the given instant, theta_pin = -90 deg (A straight down from O),
% link angle = 60 deg from vertical => B is up and to the right.
% We define rail height as the y-coordinate of B at the reference instant.

theta0_pin = -pi/2;   % initial pin angle (A straight below O, in gear frame)

% At reference instant, A position:
A0 = [0 + r_pin*cos(theta0_pin), O_y + r_pin*sin(theta0_pin)];
% A0 = [0, 3-2] = [0, 1]

% Link angle from vertical = 60 deg, so direction from A to B:
% B = A + L*[sin60, cos60] (going up-right at 60 from vertical)
link_angle0 = pi/2 - deg2rad(60);  % angle from +x axis
B0 = A0 + L*[cos(link_angle0), sin(link_angle0)];
y_rail = B0(2);   % fix rail height

fprintf('Reference B position: (%.2f, %.2f) in\n', B0(1), B0(2));
fprintf('Rail height: %.2f in\n', y_rail);

%% --- Time setup ---
% Simulate over enough time for ~1.5 revolutions of the gear
T_total = 2*pi/omega * 1.5;
dt = 0.01;
t_vec = 0:dt:T_total;
N = length(t_vec);

%% --- Precompute kinematics ---
% Gear rolls to the LEFT (omega CCW, rack below, so O moves left)
% O_x(t) = -R*omega*t  (rolling without slip: v_O = -omega*R in x)
% Actually omega CCW means top of gear moves LEFT, so O moves LEFT.
% Wait: for a gear on a fixed rack, rolling CCW means the gear translates
% to the LEFT (like a wheel spinning backward). Let's use:
%   O_x(t) = O_x0 - R*omega*t + 0.5*(-R*alpha)*t^2
% But alpha here is angular accel of the gear, not of O.
% v_O = -R*omega (leftward), a_O = -R*alpha
% We start O at x=0.

O_x = -R*omega*t_vec - 0.5*R*alpha*t_vec.^2;

% Pin angle in world frame: starts at theta0_pin, rotates by -omega*t
% (CCW omega means pin sweeps CCW, so angle increases)
theta_pin = theta0_pin + omega*t_vec + 0.5*alpha*t_vec.^2;

% Pin A world position
A_x = O_x + r_pin*cos(theta_pin);
A_y = O_y + r_pin*sin(theta_pin);  % O_y is constant

% Slider B: constrained to rail at y = y_rail, x found from link constraint
% |AB| = L  =>  (B_x - A_x)^2 + (y_rail - A_y)^2 = L^2
% B_x = A_x + sqrt(L^2 - (y_rail - A_y)^2)   [take + root: B to the right]
dy = y_rail - A_y;
dx_sq = L^2 - dy.^2;
% Clamp to avoid sqrt of negative (geometry limit)
dx_sq = max(dx_sq, 0);
B_x = A_x + sqrt(dx_sq);
B_y = y_rail * ones(size(t_vec));

%% --- Velocity kinematics (analytical) ---
% A velocity:
%   v_A = v_O + omega x r_{A/O}
%   v_O = (-R*omega - R*alpha*t) in x (scalar), 0 in y
%   omega_gear = (omega + alpha*t) CCW = (omega+alpha*t) k
%   r_{A/O} = r_pin*[cos(theta_pin), sin(theta_pin)]
omega_t = omega + alpha*t_vec;   % instantaneous gear angular speed

vO_x = -R*(omega + alpha*t_vec);
vA_x = vO_x + omega_t .* (-r_pin*sin(theta_pin));   % -omega*r_y
vA_y =        omega_t .* ( r_pin*cos(theta_pin));   %  omega*r_x

% Link AB angular velocity (omega_AB):
% v_B = v_A + omega_AB x r_{B/A}
% B is constrained: v_B = [vB_x, 0]
% r_{B/A} = [B_x-A_x, y_rail-A_y]
rBA_x = B_x - A_x;
rBA_y = B_y - A_y;

% omega_AB k x [rBA_x, rBA_y] = omega_AB*[-rBA_y, rBA_x]
% y-component: 0 = vA_y + omega_AB*rBA_x  =>  omega_AB = -vA_y/rBA_x
omega_AB = -vA_y ./ rBA_x;
vB_x = vA_x + omega_AB.*(-rBA_y);
vB_y = zeros(size(t_vec));  % constrained

%% --- Acceleration kinematics ---
% a_O = -R*alpha (constant in this problem, we use instantaneous values)
aO_x = -R*alpha * ones(size(t_vec));
alpha_t = alpha * ones(size(t_vec));  % gear angular acceleration (constant)

% a_A = a_O + alpha x r_{A/O} - omega^2 * r_{A/O}
%   alpha x r_{A/O} = alpha_t * [-r_{A/O}_y, r_{A/O}_x]
rAO_x = r_pin*cos(theta_pin);
rAO_y = r_pin*sin(theta_pin);

aA_x = aO_x + alpha_t.*(-rAO_y) - omega_t.^2.*rAO_x;
aA_y =        alpha_t.*( rAO_x) - omega_t.^2.*rAO_y;

% alpha_AB from acceleration constraint at B (a_B has only x component):
% a_B = a_A + alpha_AB x r_{B/A} - omega_AB^2 * r_{B/A}
% y: 0 = aA_y + alpha_AB*rBA_x - omega_AB^2*rBA_y
alpha_AB = (omega_AB.^2.*rBA_y - aA_y) ./ rBA_x;
aB_x = aA_x + alpha_AB.*(-rBA_y) - omega_AB.^2.*rBA_x;
aB_y = zeros(size(t_vec));

%% --- Draw gear teeth (polygon) ---
function pts = gear_shape(cx, cy, R_outer, R_inner, n_teeth, angle_offset)
    % Returns x,y points for a simple gear outline
    pts_x = []; pts_y = [];
    dtheta = 2*pi/n_teeth;
    tooth_w = dtheta*0.3;
    tooth_h = (R_outer - R_inner)*0.8;
    for i = 0:n_teeth-1
        th = angle_offset + i*dtheta;
        % root
        pts_x(end+1) = cx + R_inner*cos(th - tooth_w/2);
        pts_y(end+1) = cy + R_inner*sin(th - tooth_w/2);
        % tooth base left
        pts_x(end+1) = cx + (R_inner+tooth_h*0.1)*cos(th - tooth_w*0.4);
        pts_y(end+1) = cy + (R_inner+tooth_h*0.1)*sin(th - tooth_w*0.4);
        % tooth tip left
        pts_x(end+1) = cx + R_outer*cos(th - tooth_w*0.35);
        pts_y(end+1) = cy + R_outer*sin(th - tooth_w*0.35);
        % tooth tip right
        pts_x(end+1) = cx + R_outer*cos(th + tooth_w*0.35);
        pts_y(end+1) = cy + R_outer*sin(th + tooth_w*0.35);
        % tooth base right
        pts_x(end+1) = cx + (R_inner+tooth_h*0.1)*cos(th + tooth_w*0.4);
        pts_y(end+1) = cy + (R_inner+tooth_h*0.1)*sin(th + tooth_w*0.4);
        % root next
        pts_x(end+1) = cx + R_inner*cos(th + tooth_w/2);
        pts_y(end+1) = cy + R_inner*sin(th + tooth_w/2);
    end
    pts = [pts_x; pts_y];
end

%% --- Animation ---
fig = figure('Name','Gear-Link Animation 16-128','Color','k',...
    'Position',[100 100 1100 650]);

ax = axes('Parent',fig,'Color',[0.08 0.08 0.12],...
    'XColor','w','YColor','w','FontSize',10);
hold on; axis equal; grid on;
ax.GridColor = [0.3 0.3 0.3];
ax.GridAlpha = 0.4;

x_min = min(O_x) - R - 2;
x_max = max(O_x) + R + L + 4;
y_min = -1;
y_max = y_rail + 3;
axis([x_min x_max y_min y_max]);
xlabel('x (in)','Color','w'); ylabel('y (in)','Color','w');
title('Problem 16-128 — Gear + Link Animation','Color','w','FontSize',13);

% Scale for vectors
v_scale = 0.08;
a_scale = 0.004;

% --- Static rack ---
rack_x = linspace(x_min-5, x_max+5, 3);
rack_y = [0 0 0];
fill([x_min-5 x_max+5 x_max+5 x_min-5],[0 0 -1 -1],[0.35 0.22 0.12],'EdgeColor','none');
% Rack teeth
n_rack = ceil((x_max-x_min)/1.2)+10;
for i = 1:n_rack
    rx = x_min - 2 + (i-1)*1.0;
    fill([rx rx+0.5 rx+0.5 rx+0.7 rx+0.7 rx+1.0 rx+1.0 rx],...
         [0 0 0.25 0.25 0.25 0.25 0 0],[0.5 0.35 0.2],'EdgeColor',[0.6 0.5 0.4]);
end

% --- Static rail ---
rail_h = 0.25;
fill([x_min x_max x_max x_min],...
     [y_rail-rail_h/2 y_rail-rail_h/2 y_rail+rail_h/2 y_rail+rail_h/2],...
     [0.5 0.4 0.3],'EdgeColor',[0.7 0.6 0.5],'LineWidth',1);

% Animated objects (initialize)
n_teeth = 12;
gear_pts = gear_shape(0, O_y, R, R*0.72, n_teeth, 0);
h_gear  = fill(gear_pts(1,:), gear_pts(2,:), [0.25 0.45 0.7],...
    'EdgeColor',[0.5 0.75 1],'LineWidth',1.2);
h_hub   = fill(0,O_y,1,'FaceColor',[0.15 0.25 0.4],'EdgeColor',[0.5 0.75 1],'LineWidth',1.5);
theta_c = linspace(0,2*pi,40);
set(h_hub,'XData',0+0.4*cos(theta_c),'YData',O_y+0.4*sin(theta_c));

h_link  = plot([0 B0(1)],[1 B0(2)],'-','Color',[1 0.6 0.1],'LineWidth',3);
h_pinA  = plot(0,1,'o','MarkerSize',8,'MarkerFaceColor',[1 0.3 0.3],'MarkerEdgeColor','w','LineWidth',1.5);
h_pinB  = plot(B0(1),y_rail,'s','MarkerSize',10,'MarkerFaceColor',[1 0.3 0.3],'MarkerEdgeColor','w','LineWidth',1.5);

% Velocity arrows (quiver)
h_vA = quiver(0,1,vA_x(1)*v_scale,vA_y(1)*v_scale,0,'Color',[0.2 1 0.4],'LineWidth',2,'MaxHeadSize',0.5);
h_vB = quiver(B0(1),y_rail,vB_x(1)*v_scale,0,0,'Color',[0.2 1 0.4],'LineWidth',2,'MaxHeadSize',0.5);

% Acceleration arrows
h_aA = quiver(0,1,aA_x(1)*a_scale,aA_y(1)*a_scale,0,'Color',[1 0.3 0.8],'LineWidth',2,'MaxHeadSize',0.5);
h_aB = quiver(B0(1),y_rail,aB_x(1)*a_scale,0,0,'Color',[1 0.3 0.8],'LineWidth',2,'MaxHeadSize',0.5);

% Legend text
legend({'','','','','','Link AB','Pin A','Pin B',...
    'v_A','v_B','a_A','a_B'},...
    'TextColor','w','Color',[0.1 0.1 0.15],'EdgeColor',[0.4 0.4 0.4],...
    'Location','northeast','FontSize',9);

% Info box
h_txt = text(x_min+0.5, y_max-0.5, '', 'Color','w','FontSize',9,...
    'VerticalAlignment','top','FontName','Courier New',...
    'BackgroundColor',[0.05 0.05 0.1]);

%% --- Animate ---
for i = 1:1:N   % step by 1 for smoothest
    if ~ishandle(fig), break; end

    cx = O_x(i);
    th = theta_pin(i);

    % Update gear
    gp = gear_shape(cx, O_y, R, R*0.72, n_teeth, th);
    set(h_gear,'XData',gp(1,:),'YData',gp(2,:));
    set(h_hub,'XData',cx+0.4*cos(theta_c),'YData',O_y+0.4*sin(theta_c));

    % Update link
    set(h_link,'XData',[A_x(i) B_x(i)],'YData',[A_y(i) B_y(i)]);
    set(h_pinA,'XData',A_x(i),'YData',A_y(i));
    set(h_pinB,'XData',B_x(i),'YData',B_y(i));

    % Velocity vectors (scaled for visibility)
    set(h_vA,'XData',A_x(i),'YData',A_y(i),...
        'UData',vA_x(i)*v_scale,'VData',vA_y(i)*v_scale);
    set(h_vB,'XData',B_x(i),'YData',B_y(i),...
        'UData',vB_x(i)*v_scale,'VData',0);

    % Acceleration vectors
    set(h_aA,'XData',A_x(i),'YData',A_y(i),...
        'UData',aA_x(i)*a_scale,'VData',aA_y(i)*a_scale);
    set(h_aB,'XData',B_x(i),'YData',B_y(i),...
        'UData',aB_x(i)*a_scale,'VData',0);

    % Scroll view with gear
    ax.XLim = [cx - 8, cx + 14];

    % Info text
    info = sprintf(['t = %.2f s\n'...
                    'ω_gear = %.1f rad/s\n'...
                    '|v_A| = %.1f in/s\n'...
                    '|v_B| = %.1f in/s\n'...
                    'ω_AB = %.2f rad/s\n'...
                    '|a_A| = %.1f in/s²\n'...
                    '|a_B| = %.1f in/s²\n'...
                    'α_AB = %.1f rad/s²'],...
        t_vec(i), omega_t(i),...
        norm([vA_x(i) vA_y(i)]), abs(vB_x(i)),...
        omega_AB(i),...
        norm([aA_x(i) aA_y(i)]), abs(aB_x(i)),...
        alpha_AB(i));
    set(h_txt,'String',info,'Position',[cx-7.5, y_max-0.3]);

    drawnow limitrate;
    pause(0.04);
end

disp('Animation complete.');