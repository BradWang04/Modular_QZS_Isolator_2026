function dydt = odefun(t, y, M, g, k, b, omega, z_a)
    % Extract variables
    z_b = y(1);
    z_b_dot = y(2);
    l = y(3);
    l_dot = y(4);

    % Compute derivatives
    z_b_ddot = (M * z_a * omega^2 * sin(omega * t) - M * g + (2 * M * l * l_dot^2) / (3 * (3 * b^2 - l^2)^(3/2)) + (2 * M * l_dot^2) / (3 * sqrt(3 * b^2 - l^2))) / M;
    l_ddot = (-2 * M * z_b_ddot * l / (3 * sqrt(3 * b^2 - l^2)) + (2 * M * z_b_dot * l_dot) / (3 * (3 * b^2 - l^2)^(3/2)) + (4 * M * l_dot^2) / (9 * (3 * b^2 - l^2)^(3/2)) + (2 * M * g) / (3 * sqrt(3 * b^2 - l^2)) + (3 * k * l) / 2) / (4 * M / (9 * (3 * b^2 - l^2)) - 2 * M * l / (3 * sqrt(3 * b^2 - l^2)));

    % Return derivatives
    dydt = [z_b_dot; z_b_ddot; l_dot; l_ddot];
end