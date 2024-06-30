% initialization
clf, clear
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

% bounds
lb = [0.01, 0.01];
ub = [1, lim];

% options for optimization algorithms
options = optimset('TolX', 1.0e-3, 'MaxFunEvals', 50, 'Display', 'off');
nvars = 2;

% generate random initial points
num_initial_points = 500;
initial_points = [rand(num_initial_points, 1), rand(num_initial_points, 1) * lim];

% arrays to store results
results_fminunc = zeros(num_initial_points, 1);
results_fmincon = zeros(num_initial_points, 1);
results_fminsearch = zeros(num_initial_points, 1);
results_ga = zeros(num_initial_points, 1);
results_patternsearch = zeros(num_initial_points, 1);
results_simulannealbnd = zeros(num_initial_points, 1);
results_particleswarm = zeros(num_initial_points, 1);
results_linesearch = zeros(num_initial_points, 1);

% optimization loop
for i = 1:num_initial_points
    xq = initial_points(i, :);

    % fminunc
    [x, fval_fminunc] = fminunc(objective, xq, options);
    results_fminunc(i) = -fval_fminunc;

    % fmincon
    [x, fval_fmincon] = fmincon(objective, xq, [], [], [], [], lb, ub, [], options);
    results_fmincon(i) = -fval_fmincon;

    % fminsearch
    [x, fval_fminsearch] = fminsearch(objective, xq, options);
    results_fminsearch(i) = -fval_fminsearch;

    % line search
    [optimum_value, ~] = run_line_search(xq, 1e-10, 1e-6, objective, [0.01, 1; 0.01, lim]);
    results_linesearch(i) = optimum_value;

    % ga
    opts_ga = optimoptions('ga', 'Display', 'off');
    [x, fval_ga] = ga(objective, nvars, [], [], [], [], lb, ub, [], opts_ga);
    results_ga(i) = -fval_ga;

    % patternsearch
    opts_patternsearch = optimoptions('patternsearch', 'Display', 'off');
    [x, fval_patternsearch] = patternsearch(objective, xq, [], [], [], [], lb, ub, [], opts_patternsearch);
    results_patternsearch(i) = -fval_patternsearch;

    % simulannealbnd
    opts_simulannealbnd = optimoptions('simulannealbnd', 'Display', 'off');
    [x, fval_simulannealbnd] = simulannealbnd(objective, xq, lb, ub, opts_simulannealbnd);
    results_simulannealbnd(i) = -fval_simulannealbnd;

    % particleswarm
    opts_particleswarm = optimoptions('particleswarm', 'Display', 'off');
    [x, fval_particleswarm] = particleswarm(objective, nvars, lb, ub, opts_particleswarm);
    results_particleswarm(i) = -fval_particleswarm;
end

% boxplot for different optimization algorithms, removing outliers
figure;
boxplot([results_fminunc, results_fmincon, results_fminsearch, results_linesearch, results_ga, ...
    results_patternsearch, results_simulannealbnd, results_particleswarm], ...
    {'fminunc', 'fmincon', 'fminsearch', 'linesearch', 'ga', 'patternsearch', 'simulannealbnd', 'particleswarm'}, 'Symbol', '');
xlabel('Optimization Algorithm');
ylabel('Power Found');
title('Comparison of Optimization Algorithms');
grid on;

% function to run line search optimization
function [optimum_value, func_evals] = run_line_search(initial_point, hi, step_size, objective, gamma_range)
    again = 1;
    cycle = 0;
    f_old = Inf;
    xq = initial_point;
    options = optimset('tolx', 1.0e-8, 'MaxFunEvals', 50);
    total_func_evals = 0;

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
        xq = xnew;

        if abs(f_old - fval) < 10
            again = 0;
        end
        f_old = fval;
    end

    optimum_value = -fval;
    func_evals = total_func_evals;
end
