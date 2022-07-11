function J = DTransP(X0, tA, deqDM0, g1)
    % Compute the derivative of the singular trajectory according to the limiting system.

    %% Set Numerics
    de_opt = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

    %% Initialize Jacobian
    n = numel(X0);
    I = eye(n);
    Tspan = [0 tA];
    X0 = X0(:);
    J = 0 * I;

    %% Compute Directional Derivatives
    for jj = 1:n
        v0 = I(:, jj);
        u0 = [X0; v0; 0];
        [~, uu] = ode15s(deqDM0, Tspan, u0, de_opt);
        uuend = uu(end, :);
        v1 = uuend(end - n:end - 1);
        ffend = deqDM0(0, uuend);
        f1 = ffend(1:n);
        dw = ffend(end);
        J(:, jj) = v1(:) - (dw / g1) * f1(:);
    end

end
