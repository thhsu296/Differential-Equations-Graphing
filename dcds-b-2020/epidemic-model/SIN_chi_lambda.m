%% SIN_chi_lambda.m
% Plot chi and lambda, and then save the graphs.
% Requires 'SIN_gamma.mat' generated from 'SIN_gamma.m'.

%% Execute Example 1
N_star = 400;
f_slow = @(N) N .* (1 - N / N_star);
ex_num = 1;
fig = make_chi_lambda_fig(f_slow, ex_num);
saveas(fig, sprintf('fig_SIN_ex%d_chi_lambda', ex_num), 'png');

%% Execute Example 2
c1 = 60;
c2 = .04;
c3 = 90;
N_star = 400;
q1 = @(N) c1 * exp(- (c2 * (N - c3)).^2);
f_slow = @(N) N .* (1 - N / N_star) - q1(N);
ex_num = 2;
fig = make_chi_lambda_fig(f_slow, ex_num);
saveas(fig, sprintf('fig_SIN_ex%d_chi_lambda', ex_num), 'png');

function fig1 = make_chi_lambda_fig(f_slow, ex_num)
    %% Load gamma
    load('SIN_gamma.mat', 'S_gamma', 'I_gamma', 'N_gamma', 'params');
    [D, p, alpha, beta, r, K, N_star] = params{:};

    %% Set variables
    a = D + alpha + r;
    g = @(S) beta * S ./ (K + S);
    DhS = @(S) beta * K ./ (K + S).^2;
    N0 = fzero(@(N) g(D / (D + p) * N) - a, 200);

    %% Extract data from gamma
    numNList = size(S_gamma, 1);
    Ssol0 = S_gamma(1:end - 1, 1:end - 1);
    NList = N_gamma(:, 1);
    NOmega = N_gamma(:, end);

    %% Compute SdI
    Sd1 = S_gamma(1:end - 1, 2:end) - S_gamma(1:end - 1, 1:end - 1);
    Id1 = I_gamma(1:end - 1, 2:end) - I_gamma(1:end - 1, 1:end - 1);
    Nd1 = N_gamma(1:end - 1, 2:end) - N_gamma(1:end - 1, 1:end - 1);
    Sd2 = S_gamma(2:end, 1:end - 1) - S_gamma(1:end - 1, 1:end - 1);
    Id2 = I_gamma(2:end, 1:end - 1) - I_gamma(1:end - 1, 1:end - 1);
    Nd2 = N_gamma(2:end, 1:end - 1) - N_gamma(1:end - 1, 1:end - 1);
    Delta = Id1 .* Nd2 - Id2 .* Nd1;
    SdI = (Sd1 .* Nd2 - Sd2 .* Nd1) ./ Delta;

    ind = (Delta == 0);
    SdI(ind) = 0;

    %% Set Figure
    fig1 = figure(ex_num);
    clf;
    hold on;
    grid on;
    set(gca, 'FontSize', 16);
    plot([N0 N_star], [0 0], 'k:', 'LineWidth', 2);
    ylim([-5 3]);
    xlim([N0, N_star]);

    %% Draw boundary
    xline(N_star, 'k-', 'LineWidth', 1);

    %% Compute chi
    chi = 0 * NList;

    for jN = 1:numel(chi)
        N = NOmega(jN):.01:NList(jN);
        G = (g(D / (D + p) * N) - a) ./ f_slow(N);
        chi(jN) = trapz(N, G);
    end

    %% Draw chi
    chi1 = 10^3 * chi;
    h_chi = plot(NList, chi1, 'r-', 'LineWidth', 2);
    xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 20);

    %% Find roots of chi
    N_cut = N0 + 10;
    ind = find(chi(1:end - 1) .* chi(2:end) < 0 & NList(1:end-1) > N_cut);
    numRoots = numel(ind);
    Ssing = cell(1, numRoots);
    Ising = Ssing;
    Nsing = Ssing;

    Nx = zeros(1, numel(ind));

    for jN = 1:numel(ind)
        ind1 = ind(jN);
        N1 = NList(ind1);
        plot(N1, 0, 'b*', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
        Ssing{jN} = S_gamma(ind1, :);
        Ising{jN} = I_gamma(ind1, :);
        Nsing{jN} = N_gamma(ind1, :);
        Nx(jN) = N1;
    end

    save(sprintf('fig_SIN_ex%d_sing.mat', ex_num), 'Ssing', 'Ising', 'Nsing', 'f_slow', 'params');

    %% Compute lambda
    intLambda = @(S) -DhS(S) / alpha;
    lambda1 = arrayfun(@(jN)trapz(N_gamma(jN, 1:end - 1), intLambda(Ssol0(jN, :)) .* SdI(jN, :)), 1:numNList - 1).';
    lambda2 = log(f_slow(N_gamma(1:end - 1, 1)) ./ f_slow(N_gamma(1:end - 1, end)));
    lambda = lambda1 + lambda2;

    fprintf('\nRoots of chi for Example %d:\n', ex_num);

    for N1 = Nx
        [~, ind] = min(abs(NList - N1));
        lam = lambda(ind);
        fprintf('N1 = %.2f, lambda = %.2f\n', N1, lam);
    end

    %% Draw boundary
    tmpL = lambda;
    tmpN = NList(1:end - 1);
    h_lambda = plot(tmpN, tmpL, 'm--', 'LineWidth', 2);
    xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 20);

    %% Make Legend
    figs = [h_chi, h_lambda];
    set(h_chi, 'DisplayName', '$\chi(N)\cdot 10^3$');
    set(h_lambda, 'DisplayName', '$\lambda(N)$');
    legend(figs, 'FontSize', 22, 'Location', 'northwest', 'Interpreter', 'latex', 'AutoUpdate', 'off');

end
