% MATLAB code for line search and plotting based on given model and equations

% Initialization
clf, hold off, clear
format long

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

% Initial design point
xq = [0.5, 0.5];

% Line search initialization
again = 1;
cycle = 0;
f_old = Inf;

% Plotting setup
figure;
contourf(gamma_out, gamma_in, power_array_e, 200, 'LineColor', 'none');
colorbar;
xlabel('Gamma Out');
ylabel('Gamma In');
title('Electrical Power Distribution');
hold on;

while again >= 1
    cycle = cycle + 1;
    % Forward finite difference gradients
    hi = 1e-10;  % Reduce the perturbation size for gradient approximation
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
    % stepsize 5e-6 too larg
    sq = -df * 1e-6;  % Reduce the step size for line search

    % Line search using fminbnd
    options = optimset('tolx', 1.0e-8, 'MaxFunEvals', 50);
    [alphaq, fval] = fminbnd(@(alpha) objective(max(min(xq + alpha * sq, [1, lim]), [0.01, 0.01])), 0, 10, options);

    % Compute new design point
    xnew = max(min(xq + alphaq * sq, [1, lim]), [0.01, 0.01]);
    
    % Plot marker in current design point
    plot(xq(1), xq(2), 'o');
    plot([xq(1), xnew(1)], [xq(2), xnew(2)], 'r-');

    % Update design point
    xq = xnew;

    % Continue optimization?
    if abs(f_old - fval) < 10
        again = 0;
    end
    f_old = fval;
end

% Final optimal point marker
plot(xq(1), xq(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
text(xq(1) + 0.05, xq(2) + 0.05, sprintf('Max\n(%.2f, %.2f)', xq(1), xq(2)), ...
    'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'white');
hold off;
