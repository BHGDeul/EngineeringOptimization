% main script
% initialization
clf; clear;
format long;

% data parameters
v_w_adj = linspace(0.5, 25, 100);
A_proj_list = linspace(7, 25, 25);
F_out_list = linspace(10, 150, 20);
v_w_n = 10;
CL_out = 1.06;
A_proj = 16.65;
T_out_target = 10300;
T_out_max = 10300;
rho = 1.18;
lc = 250;
CD_out = 0.15;
CD_in = 0.10;
eff_in = 0.639;
eff_out = 0.652;
P_avg_e_req = 20000;
max_reel_speed = 25;
a_elev_out = 30 * pi / 180;
a_elev_in = 70 * pi / 180;
F_out = CL_out^3 / CD_out^2;
F_in = CD_in;
P_w = 0.5 * v_w_n^3 * rho;

number_of_steps = 100;
% gamma range
lim = max_reel_speed / v_w_n;
resolution = 100;
gamma_in = linspace(0.01, lim, resolution);
gamma_out = linspace(0.01, 1, resolution);

power_array_e = zeros(resolution, resolution);

% power array
for cj = 1:resolution
    for ci = 1:resolution
        power_array_e(ci, cj) = P_w * A_proj * ...
            (eff_out * F_out * (cos(a_elev_out) - gamma_out(cj))^2 - ...
            (F_in / eff_in) * (gamma_in(ci)^2 + 2 * cos(a_elev_in) * gamma_in(ci) + 1)) * ...
            ((gamma_out(cj) * gamma_in(ci)) / (gamma_out(cj) + gamma_in(ci)));
    end
end

% objective function
objective = @(gamma) -P_w * A_proj * ...
    (eff_out * F_out * (cos(a_elev_out) - gamma(1))^2 - ...
    (F_in / eff_in) * (gamma(2)^2 + 2 * cos(a_elev_in) * gamma(2) + 1)) * ...
    ((gamma(1) * gamma(2)) / (gamma(1) + gamma(2)));

% initial point
initial_point = [0.5, 0.5];

% contour plot
figure;
contourf(gamma_out, gamma_in, power_array_e, 200, 'LineColor', 'none');
colorbar;
xlabel('Gamma Out');
ylabel('Gamma In');
title('electrical power distribution');
hold on;

% run with original values
tic;
[optimum_value, func_evals, path] = run_optimization(initial_point, 1e-10, 1e-6, objective, [0.01, 1; 0.01, lim]);
runtime = toc;
fprintf('runtime of main optimization loop: %.4f seconds\n', runtime);

% path for normal run
plot(path(:,1), path(:,2), 'r-o', 'LineWidth', 2, 'MarkerSize', 5);
text(path(end, 1) + 0.05, path(end, 2) + 0.05, sprintf('max\n(%.2f, %.2f)', path(end, 1), path(end, 2)), ...
    'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'white');
hold off;
drawnow;

% different step sizes
step_sizes = logspace(-20, -0.3, number_of_steps);
optimum_values_step_size = zeros(size(step_sizes));
num_evaluations_step_size = zeros(size(step_sizes));
runtimes_step_size = zeros(size(step_sizes));

for i = 1:length(step_sizes)
    step_size = step_sizes(i);
    tic;
    [optimum_value, func_evals, ~] = run_optimization(initial_point, 1e-10, step_size, objective, [0.01, 1; 0.01, lim]);
    runtimes_step_size(i) = toc;
    optimum_values_step_size(i) = optimum_value;
    num_evaluations_step_size(i) = func_evals;
end

% results for step size exploration
figure;
[ax1, h1, h2] = plotyy(step_sizes, optimum_values_step_size, step_sizes, num_evaluations_step_size, 'semilogx');
xlabel('Step size');
ylabel(ax1(1), 'Optimum objective function value');
ylabel(ax1(2), 'Number of function evaluations');
title('Optimum value and function evaluations for step sizes');
grid on;
h1.LineWidth = 2;
h2.LineWidth = 2;

drawnow;

figure;
[ax3, h5, h6] = plotyy(step_sizes, optimum_values_step_size, step_sizes, runtimes_step_size, 'semilogx');
xlabel('Step size');
ylabel(ax3(1), 'Optimum objective function value');
ylabel(ax3(2), 'Runtime (seconds)');
title('Optimum value and runtime for step sizes');
grid on;
h5.LineWidth = 2;
h6.LineWidth = 2;

drawnow;

% different perturbation sizes (hi)
hi_values = logspace(-20, -1, number_of_steps);
optimum_values_hi = zeros(size(hi_values));
num_evaluations_hi = zeros(size(hi_values));
runtimes_hi = zeros(size(hi_values));

for i = 1:length(hi_values)
    hi = hi_values(i);
    tic;
    [optimum_value, func_evals, ~] = run_optimization(initial_point, hi, 1e-6, objective, [0.01, 1; 0.01, lim]);
    runtimes_hi(i) = toc;
    optimum_values_hi(i) = optimum_value;
    num_evaluations_hi(i) = func_evals;
end

% results for hi exploration
figure;
[ax2, h3, h4] = plotyy(hi_values, optimum_values_hi, hi_values, num_evaluations_hi, 'semilogx');
xlabel('Perturbation size (hi)');
ylabel(ax2(1), 'optimum objective function value');
ylabel(ax2(2), 'Number of function evaluations');
title('Optimum value and function evaluations for perturbation sizes');
grid on;
h3.LineWidth = 2;
h4.LineWidth = 2;

drawnow;

figure;
[ax4, h7, h8] = plotyy(hi_values, optimum_values_hi, hi_values, runtimes_hi, 'semilogx');
xlabel('Perturbation size (hi)');
ylabel(ax4(1), 'Optimum objective function value');
ylabel(ax4(2), 'Runtime (seconds)');
title('Optimum value and runtime for perturbation sizes');
grid on;
h7.LineWidth = 2;
h8.LineWidth = 2;

drawnow;

% function to run the main optimization loop
function [optimum_value, func_evals, path] = run_optimization(initial_point, hi, step_size, objective, gamma_range)
    again = 1;
    cycle = 0;
    f_old = Inf;
    xq = initial_point;
    options = optimset('tolx', 1.0e-8, 'MaxFunEvals', 50);
    total_func_evals = 0;
    path = xq;

    while again >= 1 && cycle < 50
        cycle = cycle + 1;
        alpha = 0.0;
        sq = [0, 0];
        fx = objective(xq);
        fx1plush = objective([xq(1) + hi, xq(2)]);
        fx2plush = objective([xq(1), xq(2) + hi]);
        dfdx1 = (fx1plush - fx) / hi;
        dfdx2 = (fx2plush - fx) / hi;
        df = [dfdx1, dfdx2];
        sq = -df * step_size;

        % line search using fminbnd
        [alphaq, fval, exitflag, output] = fminbnd(@(alpha) objective(max(min(xq + alpha * sq, gamma_range(:,2)'), gamma_range(:,1)')), 0, 10, options);
        total_func_evals = total_func_evals + output.funcCount;

        xnew = max(min(xq + alphaq * sq, gamma_range(:,2)'), gamma_range(:,1)');

        plot(xq(1), xq(2), 'o');
        plot([xq(1), xnew(1)], [xq(2), xnew(2)], 'r-');

        xq = xnew;
        path = [path; xq];

        if abs(f_old - fval) < 10
            again = 0;
        end
        f_old = fval;
    end

    optimum_value = -fval;
    func_evals = total_func_evals;
end
