function [X1, ttA, uuA] = TransP(X0, deqM0, deqM0w, which_side, Tmax)
    % Compute the singular trajectory according to the limiting system.

    %% Set Numerics
    Tspan = [0 Tmax];
    de_opt = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

    %% Compute Trajectories
    u0 = [X0(:); 0];
    [tt, uu] = ode23(deqM0, Tspan, u0, de_opt);
    tmp = uu(:, end);

    if isequal(which_side, 'min')
        tmp1 = find(tmp < min(tmp) / 2);
    elseif isequal(which_side, 'max')
        tmp1 = find(tmp > max(tmp) / 2);
    end

    if isempty(tmp1)
        error('w does not change sign');
    end

    ind0 = tmp1(end);
    tt = tt(1:ind0);
    uu = uu(1:ind0, :);

    U0 = [uu(end, :) tt(end)];
    w0 = uu(end, end);
    wspan = [w0 0];
    [~, UU] = ode23(deqM0w, wspan, U0, de_opt);
    X1 = UU(end, 1:end - 2);

    ttA = [tt; UU(:, end)];
    uuA = [uu; UU(:, 1:end - 1)];

end
