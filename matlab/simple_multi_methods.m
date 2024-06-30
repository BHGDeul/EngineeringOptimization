% initialization
clf, hold off, clear
format long

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

% gamma range
lim = max_reel_speed / v_w_n;
resolution = 100;
gamma_in = linspace(0.01, lim, resolution);
gamma_out = linspace(0.01, 1, resolution);

% power array
power_array_e = zeros(resolution, resolution);

% compute power array
for cj = 1:resolution
    for ci = 1:resolution
        power_array_e(ci, cj) = P_w * A_proj * ...
            (eff_out * F_out * (cos(a_elev_out) - gamma_out(cj))^2 - ...
            (F_in / eff_in) * (gamma_in(ci)^2 + 2 * cos(a_elev_in) * gamma_in(ci) + 1)) * ...
            ((gamma_out(cj) * gamma_in(ci)) / (gamma_out(cj) + gamma_in(ci)));
    end
end

% objective function (negated for maximization)
objective = @(gamma) -P_w * A_proj * ...
    (eff_out * F_out * (cos(a_elev_out) - gamma(1))^2 - ...
    (F_in / eff_in) * (gamma(2)^2 + 2 * cos(a_elev_in) * gamma(2) + 1)) * ...
    ((gamma(1) * gamma(2)) / (gamma(1) + gamma(2)));

% initial design point
xq = [0.5, 0.5];

% bounds
lb = [0.01, 0.01];
ub = [1, lim];

% options for optimization algorithms
options = optimset('TolX', 1.0e-3, 'MaxFunEvals', 50, 'Display', 'iter');

% subplot setup
figure;

% optimization using fminunc
subplot(2,2,1);
contourf(gamma_out, gamma_in, power_array_e, 200, 'LineColor', 'none');
colorbar;
xlabel('Gamma Out');
ylabel('Gamma In');
title('fminunc');
hold on;
plot(xq(1), xq(2), 'ro');
[x, fval_fminunc] = fminunc(objective, xq, options);
xnew_fminunc = x;
plot([xq(1) xnew_fminunc(1)], [xq(2) xnew_fminunc(2)], 'ro-');
plot(xnew_fminunc(1), xnew_fminunc(2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
text(xnew_fminunc(1), xnew_fminunc(2), sprintf('max: %.2f', -fval_fminunc), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% optimization using fmincon
subplot(2,2,2);
contourf(gamma_out, gamma_in, power_array_e, 200, 'LineColor', 'none');
colorbar;
xlabel('Gamma Out');
ylabel('Gamma In');
title('fmincon');
hold on;
plot(xq(1), xq(2), 'ro');
[x, fval_fmincon] = fmincon(objective, xq, [], [], [], [], lb, ub, [], options);
xnew_fmincon = x;
plot([xq(1) xnew_fmincon(1)], [xq(2) xnew_fmincon(2)], 'ro-');
plot(xnew_fmincon(1), xnew_fmincon(2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
text(xnew_fmincon(1), xnew_fmincon(2), sprintf('max: %.2f', -fval_fmincon), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% optimization using fminsearch
subplot(2,2,3);
contourf(gamma_out, gamma_in, power_array_e, 200, 'LineColor', 'none');
colorbar;
xlabel('Gamma Out');
ylabel('Gamma In');
title('fminsearch');
hold on;
plot(xq(1), xq(2), 'ro');
[x, fval_fminsearch] = fminsearch(objective, xq, options);
xnew_fminsearch = x;
plot([xq(1) xnew_fminsearch(1)], [xq(2) xnew_fminsearch(2)], 'ro-');
plot(xnew_fminsearch(1), xnew_fminsearch(2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
text(xnew_fminsearch(1), xnew_fminsearch(2), sprintf('max: %.2f', -fval_fminsearch), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% contour plot
subplot(2,2,4);
contourf(gamma_out, gamma_in, power_array_e, 200, 'LineColor', 'none');
colorbar;
xlabel('Gamma Out');
ylabel('Gamma In');
title('line search');
hold on;
[optimum_value, func_evals, path] = run_optimization(xq, 1e-10, 1.5e-6, objective, [0.01, 1; 0.01, lim]);
plot(path(:,1), path(:,2), 'r-o', 'LineWidth', 1, 'MarkerSize', 1);
text(path(end, 1) + 0.05, path(end, 2) + 0.05, sprintf('max\n(%.2f, %.2f)\n %.0f W', path(end, 1), path(end, 2), optimum_value), ...
    'Color', 'red', 'FontSize', 8, 'FontWeight', 'bold', 'BackgroundColor', 'white');

% subplot setup
figure;
% optimization using ga
subplot(2,2,1);
contourf(gamma_out, gamma_in, power_array_e, 200, 'LineColor', 'none');
colorbar;
xlabel('Gamma Out');
ylabel('Gamma In');
title('ga');
hold on;
plot(xq(1), xq(2), 'ro');
nvars = 2;
opts_ga = optimoptions('ga', 'Display', 'iter');
[x, fval_ga] = ga(objective, nvars, [], [], [], [], lb, ub, [], opts_ga);
xnew_ga = x;
plot([xq(1) xnew_ga(1)], [xq(2) xnew_ga(2)], 'ro-');
plot(xnew_ga(1), xnew_ga(2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
text(xnew_ga(1), xnew_ga(2), sprintf('max: %.2f', -fval_ga), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% optimization using patternsearch
subplot(2,2,2);
contourf(gamma_out, gamma_in, power_array_e, 200, 'LineColor', 'none');
colorbar;
xlabel('Gamma Out');
ylabel('Gamma In');
title('patternsearch');
hold on;
plot(xq(1), xq(2), 'ro');
opts_patternsearch = optimoptions('patternsearch', 'Display', 'iter');
[x, fval_patternsearch] = patternsearch(objective, xq, [], [], [], [], lb, ub, [], opts_patternsearch);
xnew_patternsearch = x;
plot([xq(1) xnew_patternsearch(1)], [xq(2) xnew_patternsearch(2)], 'ro-');
plot(xnew_patternsearch(1), xnew_patternsearch(2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
text(xnew_patternsearch(1), xnew_patternsearch(2), sprintf('max: %.2f', -fval_patternsearch), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% optimization using simulannealbnd
subplot(2,2,3);
contourf(gamma_out, gamma_in, power_array_e, 200, 'LineColor', 'none');
colorbar;
xlabel('Gamma Out');
ylabel('Gamma In');
title('simulannealbnd');
hold on;
plot(xq(1), xq(2), 'ro');
opts_simulannealbnd = optimoptions('simulannealbnd', 'Display', 'iter');
[x, fval_simulannealbnd] = simulannealbnd(objective, xq, lb, ub, opts_simulannealbnd);
xnew_simulannealbnd = x;
plot([xq(1) xnew_simulannealbnd(1)], [xq(2) xnew_simulannealbnd(2)], 'ro-');
plot(xnew_simulannealbnd(1), xnew_simulannealbnd(2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
text(xnew_simulannealbnd(1), xnew_simulannealbnd(2), sprintf('max: %.2f', -fval_simulannealbnd), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% optimization using particleswarm
subplot(2,2,4);
contourf(gamma_out, gamma_in, power_array_e, 200, 'LineColor', 'none');
colorbar;
xlabel('Gamma Out');
ylabel('Gamma In');
title('particleswarm');
hold on;
plot(xq(1), xq(2), 'ro');
opts_particleswarm = optimoptions('particleswarm', 'Display', 'iter');
[x, fval_particleswarm] = particleswarm(objective, nvars, lb, ub, opts_particleswarm);
xnew_particleswarm = x;
plot([xq(1) xnew_particleswarm(1)], [xq(2) xnew_particleswarm(2)], 'ro-');
plot(xnew_particleswarm(1), xnew_particleswarm(2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
text(xnew_particleswarm(1), xnew_particleswarm(2), sprintf('max: %.2f', -fval_particleswarm), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

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
