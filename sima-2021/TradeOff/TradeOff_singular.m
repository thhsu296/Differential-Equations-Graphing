%% TradeOff_singular
% Find the singular trajectory of the trade-off model
% and compute the eigenvalues of the transition function.
% Requires 'TransP.m' and 'DTransP.m'

%% Set Parameters
a =- .1;
b = 3;
c = 1;
d = 2.8;
k = 1;
r = 10;
qMin = 0;
qMax = 1;

%% Set Equations
F = @(x, q) x .* (q + r - k * x);
G = @(x, y, q) x .* y .* (a * q.^2 + b * q + c) ./ (1 + x);
D = @(y, q) d * y;
E = @(x, y, q) 1 - y .* (2 * a * q + b) ./ (1 + x);

DxF = @(x, q) -2 * k * x + q + r;
DyF = @(x, q) 0;
DxG = @(x, y, q) y * (a * q^2 + b * q + c) / (1 + x) - x * y * (a * q^2 + b * q + c) / (1 + x)^2;
DyG = @(x, y, q) x .* (a * q.^2 + b * q + c) ./ (1 + x);
DxD = @(y, q) 0;
DyD = @(y, q) d;
DxE = @(x, y, q) y .* (2 * a * q + b) ./ (1 + x).^2;
DyE = @(x, y, q) - (2 * a * q + b) ./ (1 + x);

deqM = @(x, y, q) [
                F(x, q) - G(x, y, q)
                G(x, y, q) - D(y, q)
                E(x, y, q)
                ];

deqMw = @(x, y, q) [
                (F(x, q) - G(x, y, q)) / E(x, y, q)
                (G(x, y, q) - D(y, q)) / E(x, y, q)
                1
                1 / E(x, y, q)
                ];

deqDM = @(x, y, dx, dy, q) [
                        F(x, q) - G(x, y, q)
                        G(x, y, q) - D(y, q)
                        (DxF(x, q) - DxG(x, y, q)) * dx + (DyF(x, q) - DyG(x, y, q)) * dy
                        (DxG(x, y, q) - DxD(y, q)) * dx + (DyG(x, y, q) - DyD(y, q)) * dy
                        DxE(x, y, q) * dx + DyE(x, y, q) * dy
                        ];

%% Set Numerics
de_opt = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
deltaZ = .1;

%% Preparation
deqM0 = @(t, u) deqM(u(1), u(2), qMin);
deqM0w = @(t, u) deqMw(u(1), u(2), qMin);
deqDM0 = @(t, u) deqDM(u(1), u(2), u(3), u(4), qMin);

deqM1 = @(t, u) deqM(u(1), u(2), qMax);
deqM1w = @(t, u) deqMw(u(1), u(2), qMax);
deqDM1 = @(t, u) deqDM(u(1), u(2), u(3), u(4), qMax);

%% Set Figure
fig1 = figure(1);
clf;
plot(NaN, NaN);
hold on;
grid on;

%%
X0 = [5.5707 11.0320];
Tmax = 10;
P0 = @(X0) TransP(X0, deqM0, deqM0w, 'min', Tmax);
P1 = @(X1) TransP(X1, deqM1, deqM1w, 'max', Tmax);
PP = @(X0) P1(P0(X0));
tmp = @(X0) norm(PP(X0) - X0);
XA = fminsearch(tmp, X0);
XB = P0(XA);

%% Draw Trajectory
[~, ttA, uuA] = P0(XA);
uu = uuA;
xx = uu(:, 1);
ww = uu(:, end);
h0 = plot(xx, ww, 'b', 'LineWidth', 1, 'DisplayName', '$\alpha=0$');

[~, ttB, uuB] = P1(XB);
uu = uuB;
xx = uu(:, 1);
ww = uu(:, end);
h1 = plot(xx, ww, 'r-.', 'LineWidth', 1, 'DisplayName', '$\alpha=1$');

%% Make Legend
legend([h0, h1], 'FontSize', 20, 'Interpreter', 'LaTex', 'Location', 'NorthWest');

%% Compute Stability
tA = ttA(end);
tB = ttB(end);
gB0 = E(XB(1), XB(2), qMin);
gA1 = E(XA(1), XA(2), qMax);
DP0 = DTransP(XA, tA, deqDM0, gB0);
DP1 = DTransP(XB, tB, deqDM1, gA1);
lambda = eig(DP1 * DP0);

format long
disp('Eigenvalues:')
disp(lambda);

%% Text Eigenvalues
s1 = sprintf("$\\lambda_1=%.16f$", lambda(1));
s2 = sprintf("$\\lambda_2=%.16f$", lambda(2));
text(5.2, -.58, s1, 'Interpreter', 'LaTeX', 'FontSize', 18);
text(5.2, -.75, s2, 'Interpreter', 'LaTeX', 'FontSize', 18);

%% Save Figure
filename = 'fig_TradeOff_singular2D';
saveas(fig1, filename, 'png');

% Save Solutions and Derivatives
save('TradeOff_singular.mat', 'ttA', 'uuA', 'ttB', 'uuB', 'DP0', 'DP1');
