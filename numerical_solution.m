function numerical_solution
    % Parameters
    M = 1; % Mass
    g = 9.81; % Gravity
    k = (300/1E3*g)/(80); % Spring stiffness
    b = 1; % Constant
    f = 1; % frequency
    omega = 2*pi*f; % Excitation frequency
    z_a = 1; % Excitation amplitude
    l0 = 10; % Initial value of l
    l_dot0 = 0; % Initial value of l_dot
    z_b0 = 0; % Initial value of z_b
    z_b_dot0 = 0; % Initial value of z_b_dot

    % Initial conditions
    y0 = [z_b0; z_b_dot0; l0; l_dot0];

    % Time span
    tspan = [0 2];

    % Solve the system of ODEs
    [t, y] = ode45(@(t, y) odefun(t, y, M, g, k, b, omega, z_a), tspan, y0);

    % Extract solutions
    z_b = y(:, 1);
    z_b_dot = y(:, 2);
    l = y(:, 3);
    l_dot = y(:, 4);

    % Plot results
    figure;
    subplot(2, 1, 1);
    plot(t, z_b);
    xlabel('Time (s)');
    ylabel('z_b');
    title('Displacement of z_b');

    subplot(2, 1, 2);
    plot(t, l);
    xlabel('Time (s)');
    ylabel('l');
    title('Spring Length l');
end

