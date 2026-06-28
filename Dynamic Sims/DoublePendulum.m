clear; close all; clc;

% Parameters

m1 = 1.0;   % [kg]
m2 = 1.0;   % [kg]
L1 = 1.0;   % [m]
L2 = 1.0;   % [m]
g = 9.81;   % [m/s^2]

% Vector holding the initial positions and velocities of both masses
initialConditions = [2.726117025152543;  % x1
                     0;     % v1
                     4.440052167226624;  % x2
                     0];    % v2

% This function takes in the values from the vector state x, ans uses them
% to calculate it's derivative, which will be the equation of motion to
% solve through ode45
function dxdt = derivatives(~, x, g, L1, L2, m1, m2)
    % x is the vector state that holds theta1, theta1_dot, theta2,
    % theta2_dot

    t1     = x(1);
    t1_dot = x(2);
    t2     = x(3);
    t2_dot = x(4);

    M = [(m1 + m2)*L1,      m2*L2*cos(t1-t2);
         L1*cos(t1-t2),     L2               ];

    F = [-m2*L2*(t2_dot)^2*sin(t1-t2)-(m1 + m2)*g*sin(t1);
         L1*(t1_dot)^2*sin(t1-t2)-g*sin(t2)               ];

    % Solve for the angular acceleration

    accel = M\F;

    theta1_ddot = accel(1);
    theta2_ddot = accel(2);

    % Find derivative for eom
    dxdt = [t1_dot;
             theta1_ddot;
             t2_dot;
             theta2_ddot];
end

tspan = [0, 30];
opts = odeset('RelTol',1e-6, 'AbsTol',1e-8);

% Solving for the eom
[T, X] = ode45(@(t, x) derivatives(t, x, g, L1, L2, m1, m2), tspan, initialConditions, opts);

fps = 120;
t_uniform = 0:(1/fps):T(end); % Create a perfectly spaced time grid

% Interpolate the chaotic solver outputs onto the uniform time grid
theta1 = interp1(T, X(:, 1), t_uniform);
theta2 = interp1(T, X(:, 3), t_uniform);

theta1_dot = interp1(T, X(:, 2), t_uniform);
theta2_dot = interp1(T, X(:, 4), t_uniform);

figure('Name', 'System Phase Portraits')

% Left Subplot: Upper Pendulum
subplot(1,2,1)
plot(theta1, theta1_dot, 'b-', LineWidth=2)
xlabel('\theta_1 (rad)'); ylabel('\omega_1 (rad/s)');
title('Upper Mass Phase Space')
grid on

% Right Subplot: Lower Pendulum
subplot(1,2,2)
plot(theta2, theta2_dot, 'r-', LineWidth=2)
xlabel('\theta_2 (rad)'); ylabel('\omega_2 (rad/s)');
title('Lower Mass Phase Space')
grid on

figure('Name','System time response')

subplot(2, 1, 1);
plot(t_uniform, theta1, 'b-', LineWidth=2);
ylabel("\theta_1 (rad)"); xlabel("Time (s)");
title("Upper mass time response");
grid on;

subplot(2, 1, 2);
plot(t_uniform, theta2, 'r-', LineWidth=2);
ylabel("\theta_2 (rad)"); xlabel("Time (s)");
title("Lower mass time response");
grid on;

% ----Animation----------------------

% Cartesian Coordinates
x1 = L1 * sin(theta1);
y1 = -L1 * cos(theta1);

x2 = x1 + L2 * sin(theta2);
y2 = y1 - L2 * cos(theta2);

% Animation Window
figure('Name', 'Double Pendulum Physical Space', 'Position', [100, 100, 700, 700])
ax = axes;
hold on; grid on; axis equal;

% Set axes limits slightly larger than the total combined reach of both arms
max_reach = (L1 + L2) * 1.1;
xlim([-max_reach, max_reach]); 
ylim([-max_reach, max_reach]);
xlabel('X Position (m)'); ylabel('Y Position (m)');
title('Chaotic Trajectory Animation');

% Initialize Graphics Objects (Static Frames)
% The trace line for Mass 2 (grows over time, doesn't clear)
trace = animatedline('Color', [0.5 0.5 0.5], 'LineWidth', 1.0, 'LineStyle', ':');

% The structural frame
% Link the coordinates: [Pivot_x, Mass1_x, Mass2_x]
pendulum_frame = plot([0, x1(1), x2(1)], [0, y1(1), y2(1)], ...
                      'w-', 'LineWidth', 2.0);

% The physical masses
bob1 = plot(x1(1), y1(1), 'bo', 'MarkerSize', 12, 'MarkerFaceColor', 'b');
bob2 = plot(x2(1), y2(1), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');

% Playback Animation Loop
frame_skip = 4; % Skip frames so the animation matches real-world time pacing
filename = 'DoublePendulum.gif';
delayTime = 0.05; % Time between frames in seconds (e.g., 20 frames/sec)
for k = 1:frame_skip:length(T)
    % Append the current position of mass 2 to the trace path
    addpoints(trace, x2(k), y2(k));
    
    % Update the skeletal lines connecting the joints
    set(pendulum_frame, 'XData', [0, x1(k), x2(k)], 'YData', [0, y1(k), y2(k)]);
    
    % Move the bob markers to the new joints
    set(bob1, 'XData', x1(k), 'YData', y1(k));
    set(bob2, 'XData', x2(k), 'YData', y2(k));
    drawnow;
    
    frame = getframe(gcf);          % Takes a snapshot of the current figure
    im = frame2im(frame);           % Converts the frame to an image matrix
    [imind, cm] = rgb2ind(im, 256); % Converts RGB image to an indexed image

    if k == 1
        % Create the file on the very first frame
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', delayTime);
    else
        % Append to the existing file on all subsequent frames
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
    end
end