clc; clear; close all;

%% Model
n_x = 8;
phs = PortHamiltonianSystem(@mass_matrix, @damping_matrix, @input_matrix, @potential_energy, n_x);
foot = 'left';
modes.left = @constraint_left;
modes.right = @constraint_right;
events.left = @left_foot_touchdown;
events.right = @right_foot_touchdown;

%% Simulate
x0 = [-0.2170;0.9762;0.2187;-0.3234;-16.9340;1.8667;0;0];
total_sim_time = 20;
current_time = 0;
step_count = 0;
max_steps = 4;
x_history = [];
time = [];
mode_indices = [1]; 
mode_label = {foot}; 

while(current_time < total_sim_time) && (step_count < max_steps)
    % Set constraints based on the current foot
    if strcmp(foot, 'left')
        phs.set_constraints(modes.left);
        options = odeset('Events', events.left,'RelTol',1e-5,'AbsTol',1e-5);
        dxdt = @(t,x) phs.forward_dynamics(x, [0; 0; 0; 0]);
    else
        phs.set_constraints(modes.right);
        options = odeset(options, 'Events', events.right,'RelTol',1e-5,'AbsTol',1e-5);
        dxdt = @(t,x) phs.forward_dynamics(x, [0; 0; 0; 0]);
    end
       
    % Simulate
    [t, x, te, xe, ie] = ode45(dxdt, [current_time total_sim_time], double(x0), options);
    x_history = [x_history; x];
    time = [time; t];
    
    % Check if an event occurred
    if ~isempty(te)
        fprintf('Swing foot impact at time = %0.2f\n', te(end));
        current_time = te(end);
        
        % Switch foot and apply impact mapping
        if strcmp(foot, 'left')
            foot = 'right';
            x0 = phs.impact_mapping(xe(end,:)', modes.right);
        else
            foot = 'left';
            x0 = phs.impact_mapping(xe(end,:)', modes.left);
        end
    else
        current_time = t(end);
    end

    mode_indices = [mode_indices; size(x_history, 1) + 1]; 
    mode_label{end + 1} = foot;
    
    step_count = step_count + 1;
end

%% Plot

%% Plot Configuration Variables with mode colors
fig1 = figure;
colors = {'#8cc5e3', '#1a80bb'}; % blue for left, red for right
for i = 1:4
    subplot(4, 1, i);
    hold on;
    for j = 1:length(mode_indices)-1
        idx_start = mode_indices(j);
        idx_end = mode_indices(j+1)-1;
        mode = mode_label{j};
        color = colors{strcmp(mode, 'left')+1};
        plot(time(idx_start:idx_end), x_history(idx_start:idx_end, i), 'LineWidth', 4, 'Color', color);
    end
    hold off;
    ylabel(['$q_' num2str(i) '$'], 'Interpreter', 'latex');
    grid on;
    if i == 1
        legend('Left contact', 'Right contact', 'Location', 'northwest');
        title('Configuration', 'Interpreter', 'latex');
    end
end
xlabel('Time (s)', 'Interpreter', 'latex');
format_plot(fig1);

%% Plot momentum variables with mode colors
fig2 = figure;
for i = 1:2
    subplot(2, 1, i);
    hold on;
    plot(time, x_history(:, i+4), 'k--', 'LineWidth', 2);
    for j = 1:length(mode_indices)-1
        idx_start = mode_indices(j);
        idx_end = mode_indices(j+1)-1;
        mode = mode_label{j};
        color = colors{strcmp(mode, 'left')+1};
        plot(time(idx_start:idx_end), x_history(idx_start:idx_end, i+4), 'LineWidth', 4, 'Color', color);
    end
    hold off;
    ylabel(['$p_' num2str(i) '$'],'FontSize',16, 'Interpreter', 'latex');
    grid on;
    if i == 1
        legend('Impact','Left', 'Right', 'Location', 'southeast','FontSize',10);
        title('Momentum', 'Interpreter', 'latex','FontSize',14)
    end
end
xlabel('Time (s)','FontSize',14, 'Interpreter', 'latex');

%% Plot Limit Cycle (Momentum vs. Configuration)
fig3 = figure; 
plot(x_history(:,3), x_history(:,5), 'LineWidth', 4, 'Color', colors{2});
xlabel('$q_3$', 'Interpreter', 'latex'); 
ylabel('$p_1$', 'Interpreter', 'latex');
title('Limit Cycle', 'Interpreter', 'latex');
grid on; 
format_plot(fig3);

fig4 = figure; 
plot(x_history(:,4), x_history(:,6), 'LineWidth', 4, 'Color', colors{2});
xlabel('$q_4$', 'Interpreter', 'latex'); 
ylabel('$p_2$', 'Interpreter', 'latex');
title('Limit Cycle', 'Interpreter', 'latex');
grid on; 
format_plot(fig4);

% Ensure the 'results' folder exists or create it
if ~exist('results', 'dir')
    mkdir('results');
end

% Save the figures
exportgraphics(fig1, 'results/configuration_variables.png', 'Resolution', 600);
exportgraphics(fig2, 'results/momentum_variables.png', 'Resolution', 600);
exportgraphics(fig3, 'results/limit_cycle_1.png', 'Resolution', 600);
exportgraphics(fig4, 'results/limit_cycle_2.png', 'Resolution', 600);

% %% Animate
% videoFile = 'results/robot_animation.mp4';
% v = VideoWriter(videoFile);
% v.FrameRate = 60; 
% v.Quality = 100; 
% open(v);
% 
% fig5 = figure;
% robot = importrobot('urdf/cbr.urdf');
% robot.DataFormat = 'column';
% nq = size(homeConfiguration(robot),1); 
% 
% axis([-1 10 -2 2]); 
% for i = 1:size(x_history,1)
%     show(robot, x_history(i, 1:nq)', 'PreservePlot', false);
%     drawnow;
%     frame = getframe(fig5); 
%     writeVideo(v, frame); 
% end
% 
% % Close the video file
% close(v);
% close all;

%% Helper
function format_plot(fig)
    set(findall(fig,'type','text'),'FontSize',16, 'Interpreter', 'latex');
    set(findall(fig,'type','axes'),'FontSize',16);
    set(findall(fig,'type','line'),'LineWidth',4);
end

%% Function for the mass matrix
function M = mass_matrix(q)
    mH = 10;
    m = 5;
    a = 0.5;
    b = 0.5;
    l = a + b;
    
    M = sym(zeros(length(q), length(q)));
    M(1,1) = mH + 2*m;
    M(1,2) = 0;
    M(1,3) = 0.5*m*l*cos(q(3));
    M(1,4) = 0.5*m*l*cos(q(4));
    
    M(2,1) = 0;
    M(2,2) = mH + 2*m;
    M(2,3) = 0.5*m*l*sin(q(3));
    M(2,4) = 0.5*m*l*sin(q(4));
    
    M(3,1) = 0.5*m*l*cos(q(3));
    M(3,2) = 0.5*m*l*sin(q(3));
    M(3,3) = 0.25*m*l^2;
    M(3,4) = 0;
    
    M(4,1) = 0.5*m*l*cos(q(4));
    M(4,2) = 0.5*m*l*sin(q(4));
    M(4,3) = 0;
    M(4,4) = 0.25*m*l^2;
end

%% Function for the damping matrix
function D = damping_matrix(q)
    D = zeros(length(q), length(q));
end

%% Function for the input matrix
function G = input_matrix(q)
    G = eye(length(q));
end

%% Function for the potential energy
function V = potential_energy(q)
    mH = 10;
    m = 5;
    a = 0.5;
    b = 0.5;
    g = 9.81;
    V = (2*m+mH)*g*q(2) - m*g*a*cos(q(3)) - m*g*b*cos(q(4));
end

%% Function for the constraints for the left foot
function G_cl = constraint_left(q)
    G_cl = [1, 0;
            0, 1;
            cos(q(3)), sin(q(3));
            0, 0];
end

%% Function for the constraints for the right foot
function G_cr = constraint_right(q)
    G_cr = [1, 0;
            0, 1;
            0, 0; 
            cos(q(4)), sin(q(4))];
end

%% Function for the left foot touchdown event 
function [value, isterminal, direction] = left_foot_touchdown(t,x)
    q1     = x(1);
    q2     = x(2);
    q3     = x(3);
    q4     = x(4);
    
    Hws = [cos(q4), 0, -sin(q4), q1 + sin(q4);
                  0, 1,        0,            0;
            sin(q4), 0,  cos(q4), q2 - cos(q4);
                  0, 0,        0,            1];

    psi= -deg2rad(3);

    Hgw = [cos(psi), 0,sin(psi),0; ...
              0,1,0,0; ...
              -sin(psi), 0, cos(psi),0;...
              0,0,0,1
              ];

    Hgs = Hgw*Hws;
    value = Hgs(3,4);
    direction = -1;
    Hwp = [cos(q3), 0, -sin(q3), q1 + sin(q3);
          0, 1,        0,            0;
          sin(q3), 0,  cos(q3), q2 - cos(q3);
          0, 0,        0,            1];
    Hgp = Hgw*Hwp;

    % Only if swing foot in front of planted foot, terminate
    if (Hgs(1,4) - Hgp(1,4)) > 0
        isterminal = 1;
    else
       isterminal = 0;
    end
end

%% Function for the right foot touchdown event 
function [value, isterminal, direction] = right_foot_touchdown(t,x)
    q1     = x(1);
    q2     = x(2);
    q3     = x(3);
    q4     = x(4);
    
    Hws = [cos(q3), 0, -sin(q3), q1 + sin(q3);
           0, 1,        0,            0;
           sin(q3), 0,  cos(q3), q2 - cos(q3);
           0, 0,        0,            1];

    psi = -deg2rad(3);

    Hgw = [cos(psi), 0,sin(psi),0; ...
              0,1,0,0; ...
              -sin(psi), 0, cos(psi),0;...
              0,0,0,1
              ];

    Hgs = Hgw*Hws;
    value = Hgs(3,4);
    direction = -1;
    Hwp = [cos(q4), 0, -sin(q4), q1 + sin(q4);
              0, 1,        0,            0;
        sin(q4), 0,  cos(q4), q2 - cos(q4);
              0, 0,        0,            1];
    Hgp = Hgw*Hwp;

    % Only if swing foot in front of planted foot, terminate
    if (Hgs(1,4) - Hgp(1,4)) > 0
        isterminal = 1;
    else
       isterminal = 0;
    end
end