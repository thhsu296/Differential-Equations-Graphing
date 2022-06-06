%% SIN_epsilon.m
% Plot trajectories of the system with positive perturbation parameter epsilon.
% Requires 'fig_SIN_ex1_sing.mat' and 'fig_SIN_ex2_sing.mat' generated from 'SIN_chi_lambda.m'.

%% Execute Example 1
ex_num = 1;
init_values = {{60, 2, 120}};
fig1 = make_epsilon_figure(ex_num, 1e-5, init_values);
saveas(fig1, sprintf('fig_SIN_ex%d_epsilon', ex_num), 'png');

%% Execute Example 2
ex_num = 2;
init_values = {{40, 1.3, 80}, {40, 2.5, 80}};
fig2 = make_epsilon_figure(ex_num, 1e-5, init_values);
saveas(fig2, sprintf('fig_SIN_ex%d_epsilon', ex_num), 'png');

function fig1 = make_epsilon_figure(ex_num, epsilon, init_values)
    %% Load singular trajectories
    load(sprintf('fig_SIN_ex%d_sing.mat', ex_num), 'Ssing', 'Ising', 'Nsing', 'f_slow', 'params');
    [D, p, alpha, beta, r, K, ~] = params{:};

    %% Set Equations
    a = D + alpha + r;
    g = @(S, N) beta * S / (K + S);

    deq0 = @(S, I, N) [
                    D * N + epsilon * f_slow(N) - g(S, N) * I - (D + p) * S
                    (g(S, N) - a) * I
                    epsilon * f_slow(N) - alpha * I
                    ];

    %% Set Numerics
    de_opt = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
    umax = 10000;

    %% Preparation
    deq = @(t, u) (max(abs(u)) < umax) * deq0(u(1), u(2), u(3));

    %% Set Figure
    fig1 = figure(ex_num);
    clf;
    plot3(NaN, NaN, NaN);
    hold on;
    grid on;
    pbaspect([2 2 1.5]);
    xlim([0 400]);
    ylim([0 8.2]);
    zlim([0, 400]);
    set(gca, 'BoxStyle', 'full');
    set(gca, 'FontSize', 14);
    view([-68 25]);
    set(gca, 'Ydir', 'reverse')
    xticks(0:100:400);
    set(gca, 'XTickLabel', {'0', '', '200', '', '400'});

    %% Axes
    text(250, 9.5, 0, '$S$', 'Interpreter', 'latex', 'FontSize', 20);
    text(-120, 3.05, 0, '$I$', 'Interpreter', 'latex', 'FontSize', 20);
    text(0, -1, 460, '$N$', 'Interpreter', 'latex', 'FontSize', 22);

    %% Draw Z
    N = [0 450];
    I = 0 * N;
    S = D / (D + p) * N;
    plot3(S, I, N, 'k:', 'LineWidth', 2);

    %% Draw trajectories
    colors = 'rb';

    for iter = 1:numel(init_values)
        [Si, Ii, Ni] = init_values{iter}{:};
        Tmax = 18 / epsilon;
        w0 = [Si, Ii, Ni];
        [~, w] = ode23(deq, [0 Tmax], w0, de_opt);
        S = w(:, 1);
        I = w(:, 2);
        N = w(:, 3);
        plot3(S(1), I(1), N(1), 'o', 'Color', colors(iter));
        plot3(S, I, N, '-.', 'Color', colors(iter), 'LineWidth', 2);
    end

    %% Draw Singular
    for jN = 1:numel(Nsing)
        S = Ssing{jN};
        I = Ising{jN};
        N = Nsing{jN};
        plot3(S, I, N, 'k-', 'LineWidth', 2);
    end

end
