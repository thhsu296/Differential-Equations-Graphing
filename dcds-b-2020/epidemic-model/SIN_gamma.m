%% SIN_gamma.m
% Plot and save heteroclinic orbits (gamma).

%% Set parameters
D = .2;
p = .01;
alpha = .048;
beta = 1;
r = .75;
K = .1;
N_star = 400;
params = {D, p, alpha, beta, r, K, N_star};

%% Plot heteroclic orbits
[fig1, ~, ~, ~] = get_gamma(7, params);
set(gca, 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none');
saveas(fig1, 'fig_SIN_gamma', 'png');

%% Save heteroclic orbits
[~, S_gamma, I_gamma, N_gamma] = get_gamma(500, params);
save('SIN_gamma.mat', 'S_gamma', 'I_gamma', 'N_gamma', 'params');

function [fig1, S_gamma, I_gamma, N_gamma] = get_gamma(num_init, params)
    %% Get the given number of heteroclinic orbits.
    % Input num_init: the number of initial points.

    %% Parse parameters
    [D, p, alpha, beta, r, K, N_star] = params{:};
    
    %% Set equations
    a = D + alpha + r;
    g = @(S, N) beta * S / (K + S);

    deqI0 = @(S, I, N) [
                    (D * N - g(S, N) * I - (D + p) * S) / (abs(g(S, N) - a) * I)
                    sign(g(S, N) - a)
                    -alpha / abs(g(S, N) - a)
                    ];

    deqN0 = @(S, I, N) [
                    (D * N - g(S, N) * I - (D + p) * S) / (alpha * I)
                    (g(S, N) - a) / alpha
                    -1
                    ];

    %% Set numerics
    de_opt = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
    I_cut = 1e-3;

    %% Set figure
    fig1 = figure(num_init);
    clf;
    plot3(NaN, NaN, NaN);
    hold on;
    grid on;
    pbaspect([2 2 1.2]);
    xlim([0 450]);
    ylim([0 10]);
    zlim([0, 450]);
    set(gca, 'BoxStyle', 'full');
    set(gca, 'FontSize', 14);
    view([-72 33]);
    set(gca, 'Ydir', 'reverse')

    %% Draw Axes
    plot3([0 400], [0 0], [0 0], '-k', 'LineWidth', 1);
    plot3([0 0], [0 10], [0 0], '-k', 'LineWidth', 1);
    plot3([0 0], [0 0], [0 400], '-k', 'LineWidth', 1);
    text(400, 0, 0, '$S$', 'Interpreter', 'latex', 'FontSize', 22);
    text(0, 10.05, 0, '$I$', 'Interpreter', 'latex', 'FontSize', 22);
    text(-40, 0, 460, '$N$', 'Interpreter', 'latex', 'FontSize', 22);

    %% Draw Z
    N = [0 400];
    I = 0 * N;
    S = D / (D + p) * N;
    h_Z = plot3(S, I, N, 'k:', 'LineWidth', 2, 'DisplayName', '${{\mathcal{Z}}_{0}}$');

    %% Set Functions and Points
    N0 = fzero(@(N) g(D / (D + p) * N, N) - a, N_star / 2);

    S_null = @(N) D / (D + p) * N;
    S0 = S_null(N0);
    S_star = S_null(N_star);

    %% Draw points
    plot3(S0, 0, N0, 'ro', 'Markersize', 6, 'MarkerFAceColor', 'r');
    plot3(S_star, 0, N_star, 'r*', 'Markersize', 6);
    text(S0, - .9, N0 + 25, '$N_0$', 'Interpreter', 'latex', 'FontSize', 20, 'Color', 'r');
    text(S_star, - .9, N_star + 25, '$N_{\mathrm{max}}$', 'Interpreter', 'latex', 'FontSize', 20, 'Color', 'r');

    %% Set rescaled DEs
    L_upper = @(I, N) N - N0 - 40 * I;
    L_lower = @(I, N) N - N0 + 2 * I;
    deqIupper = @(t, u) (L_upper(u(2), u(3)) > 0) * deqI0(u(1), u(2), u(3));
    deqN = @(t, u) (L_lower(u(2), u(3)) > 0) * deqN0(u(1), u(2), u(3));
    deqIlower = @(t, u) (u(2) > I_cut) * deqI0(u(1), u(2), u(3));

    %% Set initial values
    N_init_list = linspace(N0, N_star, num_init + 1);
    N_init_list(1) = [];

    numIspan = 400;
    Ispan = linspace(0, 10, numIspan);

    numNspan = 1000;
    Nspan = linspace(0, 400, numNspan);

    len_gamma = numNspan + 2 * numIspan;
    S_gamma = NaN(num_init, len_gamma);
    I_gamma = S_gamma;
    N_gamma = S_gamma;

    I_init = 1e-4;
    S_center = @(I, N) S_null(N) + ((alpha * D) / (D + p) - g(D / (D + p) * N, N)) / (D + p + g(D / (D + p) * N, N) - a) * I;

    %% Run simulation
    for jN = 1:num_init
        % Get the upper segment
        N_init = N_init_list(jN);
        S_init = S_center(I_init, N_init);
        w0 = [S_init, I_init, N_init];
        [~, wA] = ode23(deqIupper, Ispan, w0, de_opt);

        % Get the middle segment
        w0 = wA(end, :);
        [~, wB] = ode23(deqN, Nspan, w0, de_opt);

        % Get the lower segment
        w0 = wB(end, :);
        [~, wC] = ode23(deqIlower, Ispan, w0, de_opt);

        % Plot trajectory
        w = [wA; wB; wC];
        S = w(:, 1);
        I = w(:, 2);
        N = w(:, 3);
        h_gamma = plot3(S, I, N, 'b', 'LineWidth', 2, 'DisplayName', '$\gamma(N_1)$');
        drawnow();

        % Save the trajectory into solution matrices
        S_gamma(jN, :) = S.';
        I_gamma(jN, :) = I.';
        N_gamma(jN, :) = N.';
        fprintf('Run %d of %d: N_init = %.1f\n', jN, num_init, N_init);
    end

    %% Set legend
    figs = [h_Z, h_gamma];
    legend(figs, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeast', 'AutoUpdate', 'off');
end
