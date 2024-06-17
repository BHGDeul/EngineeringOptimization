% Main Script

% Initialization
clf; clear;
format long;

% Define data parameters
data.v_w_adj = linspace(0.5, 25, 100);
data.A_proj_list = linspace(7, 25, 25);
data.F_out_list = linspace(10, 150, 20);
data.v_w_n = 10;
data.CL_out = 1.06;
data.A_proj = 16.65;
data.T_out_target = 10300;
data.T_out_max = 10300;
data.rho = 1.18;
data.lc = 250;
data.CD_out = 0.15;
data.CD_in = 0.10;
data.eff_in = 0.639;
data.eff_out = 0.652;
data.P_avg_e_req = 20000;
data.max_reel_speed = 25;
data.a_elev_out = 30 * pi / 180;
data.a_elev_in = 70 * pi / 180;
data.F_out = data.CL_out^3 / data.CD_out^2;
data.F_in = data.CD_in;
data.P_w = 0.5 * data.v_w_n^3 * data.rho;

% Define gamma range
lim = data.max_reel_speed / data.v_w_n;
resolution = 100;
gamma_in = linspace(0.01, lim, resolution);
gamma_out = linspace(0.01, 1, resolution);

% Initialize power array
power_array_e = zeros(resolution, resolution);

% Compute power array
for cj = 1:resolution
    for ci = 1:resolution
        power_array_e(ci, cj) = data.P_w * data.A_proj * ...
            (data.eff_out * data.F_out * (cos(data.a_elev_out) - gamma_out(cj))^2 - ...
            (data.F_in / data.eff_in) * (gamma_in(ci)^2 + 2 * cos(data.a_elev_in) * gamma_in(ci) + 1)) * ...
            ((gamma_out(cj) * gamma_in(ci)) / (gamma_out(cj) + gamma_in(ci)));
    end
end

% Objective function (negated for maximization)
objective = @(gamma) -data.P_w * data.A_proj * ...
    (data.eff_out * data.F_out * (cos(data.a_elev_out) - gamma(1))^2 - ...
    (data.F_in / data.eff_in) * (gamma(2)^2 + 2 * cos(data.a_elev_in) * gamma(2) + 1)) * ...
    ((gamma(1) * gamma(2)) / (gamma(1) + gamma(2)));

% Define initial point
initial_point = [0.5, 0.5];

% Plotting setup for contour plot
figure;
contourf(gamma_out, gamma_in, power_array_e, 200, 'LineColor', 'none');
colorbar;
xlabel('Gamma Out');
ylabel('Gamma In');
title('Electrical Power Distribution');
hold on;

% Normal run with original values
[optimum_value, func_evals, path] = run_optimization(initial_point, 1e-10, 1e-6, objective, [0.01, 1; 0.01, lim]);

% Plot path for normal run
plot(path(:,1), path(:,2), 'r-o', 'LineWidth', 2, 'MarkerSize', 5);
text(path(end, 1) + 0.05, path(end, 2) + 0.05, sprintf('Max\n(%.2f, %.2f)', path(end, 1), path(end, 2)), ...
    'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'white');
hold off;

% Ensure the figure is drawn
drawnow;

% Explore different step sizes
step_sizes = logspace(-20, -0.3, 100);
optimum_values_step_size = zeros(size(step_sizes));
num_evaluations_step_size = zeros(size(step_sizes));

for i = 1:length(step_sizes)
    step_size = step_sizes(i);
    [optimum_value, func_evals, ~] = run_optimization(initial_point, 1e-10, step_size, objective, [0.01, 1; 0.01, lim]);
    optimum_values_step_size(i) = optimum_value;
    num_evaluations_step_size(i) = func_evals;
end

% Plot the results for step size exploration
figure;
[ax1, h1, h2] = plotyy(step_sizes, optimum_values_step_size, step_sizes, num_evaluations_step_size, 'semilogx');
xlabel('Step Size');
ylabel(ax1(1), 'Optimum Objective Function Value');
ylabel(ax1(2), 'Number of Function Evaluations');
title('Optimum Objective Function Value and Function Evaluations for a Range of Step Sizes');
grid on;
h1.LineWidth = 2;
h2.LineWidth = 2;

% Ensure the figure is drawn
drawnow;

% Explore different perturbation sizes (hi)
hi_values = logspace(-20, -1, 100);  % Range of hi values
optimum_values_hi = zeros(size(hi_values));
num_evaluations_hi = zeros(size(hi_values));

for i = 1:length(hi_values)
    hi = hi_values(i);
    [optimum_value, func_evals, ~] = run_optimization(initial_point, hi, 1e-6, objective, [0.01, 1; 0.01, lim]);
    optimum_values_hi(i) = optimum_value;
    num_evaluations_hi(i) = func_evals;
end

% Plot the results for hi exploration
figure;
[ax2, h3, h4] = plotyy(hi_values, optimum_values_hi, hi_values, num_evaluations_hi, 'semilogx');
xlabel('Perturbation Size (hi)');
ylabel(ax2(1), 'Optimum Objective Function Value');
ylabel(ax2(2), 'Number of Function Evaluations');
title('Optimum Objective Function Value and Function Evaluations for a Range of Perturbation Sizes');
grid on;
h3.LineWidth = 2;
h4.LineWidth = 2;

% Ensure the figure is drawn
drawnow;

% Function to run the main optimization loop
function [optimum_value, func_evals, path] = run_optimization(initial_point, hi, step_size, objective, gamma_range)
    % Line search initialization
    again = 1;
    cycle = 0;
    f_old = Inf;
    xq = initial_point;
    options = optimset('tolx', 1.0e-8, 'MaxFunEvals', 50);
    total_func_evals = 0;
    path = xq;

    while again >= 1 && cycle < 50
        cycle = cycle + 1;
        % Forward finite difference gradients
        alpha = 0.0;
        sq = [0, 0];
        % Objective function at current point
        fx = objective(xq);
        % Perturbated objective function values
        fx1plush = objective([xq(1) + hi, xq(2)]);
        fx2plush = objective([xq(1), xq(2) + hi]);
        % Objective function derivatives
        dfdx1 = (fx1plush - fx) / hi;
        dfdx2 = (fx2plush - fx) / hi;
        % Gradient vector
        df = [dfdx1, dfdx2];
        % Steepest descent search direction
        sq = -df * step_size;  % Step size for line search

        % Line search using fminbnd
        [alphaq, fval, exitflag, output] = fminbnd(@(alpha) objective(max(min(xq + alpha * sq, gamma_range(:,2)'), gamma_range(:,1)')), 0, 10, options);
        total_func_evals = total_func_evals + output.funcCount;

        % Compute new design point
        xnew = max(min(xq + alphaq * sq, gamma_range(:,2)'), gamma_range(:,1)');
        
        % Plot marker in current design point
        plot(xq(1), xq(2), 'o');
        plot([xq(1), xnew(1)], [xq(2), xnew(2)], 'r-');

        % Update design point
        xq = xnew;
        path = [path; xq];

        % Continue optimization?
        if abs(f_old - fval) < 10
            again = 0;
        end
        f_old = fval;
    end

    optimum_value = -fval;  % Negate back for maximization
    func_evals = total_func_evals;
end
