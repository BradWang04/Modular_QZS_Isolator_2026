%% QZS Isolator: Comparison of Linear, Hardening, and Softening Systems
% This script demonstrates and compares three types of vibration isolators:
% 1. A standard linear isolator.
% 2. A hardening QZS isolator (like the original geometric model).
% 3. A physically-plausible softening QZS isolator (e.g., a magnetic-spring model).

clear all; clc; close all;

%% ====================== 1. Common System Parameters ======================
M = 0.548;          % Mass (kg), taken from your previous model for consistency
zeta = 0.03;        % Damping ratio, kept low to clearly show nonlinear effects
W = 2e-3;           % Base excitation amplitude (m)
freq_range = [0, 10]; % Focus on the low-frequency range where effects are prominent
n_freq = 500;       % High resolution for smooth curves

%% ====================== 2. Define System Properties for Each Case ======================

% We will create a "struct" to hold the parameters for each case
systems = struct();

% --- Case 1: Standard Linear Isolator ---
% Let's design it to have a natural frequency around 3 Hz
f_n_linear = 3.0;
systems(1).name = 'Linear System';
systems(1).k1 = M * (2*pi*f_n_linear)^2; % k1 > 0
systems(1).k3 = 0;                      % k3 = 0
systems(1).color = [0.5 0.5 0.5];       % Gray color
systems(1).linestyle = ':';

% --- Case 2: Hardening QZS Isolator (Represents your original model) ---
% Designed for a very low linear natural frequency (QZS condition)
f_n_qzs = 0.5; % Target low natural frequency
systems(2).name = 'Hardening QZS (k3 > 0)';
systems(2).k1 = M * (2*pi*f_n_qzs)^2;    % k1 is small and positive
systems(2).k3 = 2.7e5;                  % A significant POSITIVE nonlinear term (from your model)
systems(2).color = 'b';
systems(2).linestyle = '-';

% --- Case 3: Softening QZS Isolator (Represents the new magnetic model) ---
% Same QZS condition (low k1), but with a negative k3
systems(3).name = 'Softening QZS (k3 < 0)';
systems(3).k1 = M * (2*pi*f_n_qzs)^2;    % k1 is small and positive, same as hardening case
systems(3).k3 = -2.7e5;                 % A significant NEGATIVE nonlinear term
systems(3).color = 'r';
systems(3).linestyle = '--';

%% ====================== 3. Dynamic Analysis Loop ======================
figure('Position', [100, 100, 900, 520]);
hold on; box on;

% Loop through each system defined above
for s_idx = 1:length(systems)
    
    current_system = systems(s_idx);
    fprintf('\n--- Calculating for: %s ---\n', current_system.name);
    
    k1 = current_system.k1;
    k3 = current_system.k3;
    
    % Damping coefficient 'c' is based on the linear part of the stiffness
    c = 2 * zeta * sqrt(abs(k1) * M);
    
    % Frequency array
    f_values = linspace(1e-6, freq_range(2), n_freq);
    omega_values = 2*pi*f_values;
    Tr_values_dB = zeros(size(f_values));
    
    U_prev = W;
    options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12);
    
    for i = 1:length(f_values)
        omega = omega_values(i);
        
        % The same Harmonic Balance equation works for all cases!
        HB_eq = @(U) ((k1 - M*omega^2 + 0.75*k3*U.^2).^2 + (c*omega).^2).*U.^2 - (M*omega^2*W).^2;
        
        U_linear = abs(M*omega^2*W) / sqrt((k1 - M*omega^2)^2 + (c*omega)^2);
        U0 = max([U_prev, U_linear, 1e-9]);
        [U_sol, ~, exitflag] = fsolve(HB_eq, U0, options);
        if exitflag <= 0, U_sol = U_linear; end
        
        U = abs(U_sol);
        U_prev = U;
        k_eff = k1 + 0.75*k3*U^2;
        phi = atan2(c*omega, (k_eff - M*omega^2));
        Tr_disp = sqrt(1 + (U/W)^2 + 2*(U/W)*cos(phi));
        Tr_values_dB(i) = 20*log10(Tr_disp);
    end
    
    % Plot the result for the current system
    plot(f_values, Tr_values_dB, 'Color', current_system.color, ...
         'LineStyle', current_system.linestyle, 'LineWidth', 2.5, ...
         'DisplayName', current_system.name);
end

%% ====================== 4. Finalize Plot ======================
yline(0, '-k', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Transmissibility (dB)', 'FontSize', 14, 'FontWeight', 'bold');
title('Comparison of Linear, Hardening, and Softening Isolators', 'FontSize', 16, 'FontWeight', 'bold');
xlim(freq_range);
ylim([-40, 30]);
grid on; grid minor;
legend('Location', 'northeast', 'FontSize', 12);
set(gca, 'FontSize', 12);
hold off;