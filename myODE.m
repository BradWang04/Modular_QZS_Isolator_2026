function dydt = myODE(t, y)
    % y(1) 对应 y(t)
    % y(2) 对应 v(t) = y'(t)
    dydt = [y(2); % y'(t) = v(t)
            f(t, y(1), y(2))]; % v'(t) = f(t, y(t), v(t))
end



