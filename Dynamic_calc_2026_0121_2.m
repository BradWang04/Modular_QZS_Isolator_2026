%% QZS Vibration Isolator - Transmissibility Analysis (Corrected Version + Required Plots)
% Adds:
% 1) 0-90 Hz transmissibility curve vs frequency
% 2) Steady-state response: acceleration vs time (at a chosen frequency)
% 3) Comparison between theoretical transmissibility and experimental data
%
% Note:
% - This script computes displacement transmissibility (in dB) based on HB solution for U and phase phi.
% - For steady-state acceleration time history, we reconstruct u(t)=U cos(wt-phi),
%   then z_top(t)=u(t)+W cos(wt), and a_top(t)=d2(z_top)/dt2.

clear all; clc; close all;

%% ====================== 1. System Parameters ======================
% Geometric parameters (SI units: m)
b = 80e-3;          % Base radius (m)
l0 = 80e-3;         % Natural length of springs (m), l0 = b for QZS
x0 = 5e-3;          % Initial displacement offset (m)
g = 9.81;           % Gravitational acceleration (m/s^2)

% Spring stiffness (N/m)
k = 300;            % Individual spring stiffness (N/m)

% Damping parameters
zeta = 0.0311;      % Damping ratio (from experiment)

% Excitation parameters
W = 2e-3;           % Base excitation amplitude (m)
freq_range = [0, 90];
n_freq = 200;

%% ====================== 2. Static Equilibrium Analysis ======================
% Initial angle (unloaded position)
phi_0 = asin((b - x0) / b);

% Angular range for static analysis
phi_array = linspace(phi_0, asin((b - 50e-3) / b), 200);

% Displacement from initial position
z_disp = 2*b*(sin(phi_0) - sin(phi_array));

% Vertical restoring force (from Eq.18 in manuscript)
% F_v = (9k*b*sin(phi)/2) * (cos(phi) - cos(phi_0)) / cos(phi)
F_Ver = (9*k*b/2) .* sin(phi_array) .* (cos(phi_array) - cos(phi_0)) ./ cos(phi_array);

% Determine equilibrium position under payload
% Here you selected the maximum restoring force point as the supported load.
[F_max, idx_eq] = max(F_Ver);
phi_eq = phi_array(idx_eq);
z_eq = z_disp(idx_eq);

% Equivalent mass
M = F_max / g;

fprintf('System Parameters:\n');
fprintf('  Equilibrium angle φ_eq = %.4f rad (%.2f°)\n', phi_eq, phi_eq*180/pi);
fprintf('  Supported mass M = %.4f kg\n', M);
fprintf('  Equilibrium position z_eq = %.4f mm\n', z_eq*1000);

%% ====================== 3. Linearized Stiffness Coefficients ======================
% From manuscript Eq.(24-25):
% k1 = (9k/2) * [1 - cos(phi_0)/cos(phi_eq)]
% k3 = (9k*cos(phi_0)) / (4*b^2*cos^3(phi_eq))
k1 = (9*k/2) * (1 - cos(phi_0)/cos(phi_eq));
k3 = (9*k*cos(phi_0)) / (4*b^2 * cos(phi_eq)^3);

% Natural frequency (use abs(k1) for safety if k1<0)
omega_n = sqrt(abs(k1)/M);
f_n = omega_n / (2*pi);

% Damping coefficient (linearized about k1)
c = 2*zeta*sqrt(abs(k1)*M);

fprintf('  Linear stiffness k1 = %.4f N/m\n', k1);
fprintf('  Nonlinear stiffness k3 = %.4e N/m^3\n', k3);
fprintf('  Natural frequency f_n = %.4f Hz\n', f_n);
fprintf('  Damping coefficient c = %.4f N·s/m\n', c);

%% ====================== 4. Transmissibility Calculation (HB) ======================
% Frequency array (avoid exactly 0 Hz in omega to prevent divide-by-zero in HB initial guess)
f_values = linspace(max(freq_range(1), 1e-6), freq_range(2), n_freq);
omega_values = 2*pi*f_values;

% Initialize arrays
Tr_values_dB = zeros(size(f_values));   % displacement transmissibility (dB)
Tr_force_dB  = zeros(size(f_values));   % force transmissibility (dB), optional
U_values = zeros(size(f_values));       % relative amplitude
phi_values = zeros(size(f_values));     % phase (relative response lag)

fprintf('\nCalculating transmissibility...\n');

U_prev = W; % continuation initial
options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12, 'TolX', 1e-12, 'MaxIterations', 200);

for i = 1:length(f_values)
    omega = omega_values(i);

    % HB amplitude equation:
    % [(k1 - M*ω² + 3/4*k3*U²)² + (cω)²]*U² = (Mω²W)²
    HB_eq = @(U) ((k1 - M*omega^2 + 0.75*k3*U.^2).^2 + (c*omega).^2).*U.^2 - (M*omega^2*W).^2;

    % Initial guess: use previous solution (continuation) + linear estimate as fallback
    U_linear = abs(M*omega^2*W) / sqrt((k1 - M*omega^2)^2 + (c*omega)^2);
    U0 = max([U_prev, U_linear, 1e-9]);

    [U_sol, ~, exitflag] = fsolve(HB_eq, U0, options);
    if exitflag <= 0 || ~isfinite(U_sol)
        % Try alternative guesses
        guess_set = [U_linear, W, 0.1*W, 10*W, U0];
        U_sol = NaN;
        for gg = guess_set
            [tmp, ~, ef] = fsolve(HB_eq, max(gg,1e-9), options);
            if ef > 0 && isfinite(tmp)
                U_sol = tmp;
                break;
            end
        end
        if isnan(U_sol)
            U_sol = U0; % fallback (won't crash)
        end
    end

    U = abs(U_sol);
    U_values(i) = U;
    U_prev = U;

    % Effective stiffness at amplitude U
    k_eff = k1 + 0.75*k3*U^2;

    % Phase: u(t)=U cos(ωt-φ)
    phi = atan2(c*omega, (k_eff - M*omega^2));
    phi_values(i) = phi;

    % Force transmissibility (optional)
    Tr_force = sqrt(k_eff^2 + (c*omega)^2) / sqrt((k_eff - M*omega^2)^2 + (c*omega)^2);
    Tr_force_dB(i) = 20*log10(Tr_force);

    % Displacement transmissibility for base displacement excitation:
    % Z_top = u + z_b, with z_b=W cos(ωt)
    % => |Z_top|/W = sqrt(1 + (U/W)^2 + 2(U/W)cosφ)
    Tr_disp = sqrt(1 + (U/W)^2 + 2*(U/W)*cos(phi));
    Tr_values_dB(i) = 20*log10(Tr_disp);
end

fprintf('Calculation complete.\n');

%% ====================== 5. Plot 1: 0-90 Hz Transmissibility curve vs frequency ======================
figure('Position', [100, 100, 900, 520]);
hold on; box on;

plot(f_values, Tr_values_dB, 'b-', 'LineWidth', 2, 'DisplayName', 'Theory (HB, disp. Tr)');

yline(0, '--k', 'LineWidth', 1, 'DisplayName', '0 dB');
% yline(-3, ':r', 'LineWidth', 1, 'DisplayName', '-3 dB');

% resonance marker
[Tr_max_dB, idx_res] = max(Tr_values_dB);
f_res = f_values(idx_res);
plot(f_res, Tr_max_dB, 'ro', 'MarkerSize', 9, 'MarkerFaceColor', 'r', ...
    'DisplayName', sprintf('Peak: %.1f Hz', f_res));

% isolation onset: first < 0 dB
idx_isolation = find(Tr_values_dB < 0, 1, 'first');
if ~isempty(idx_isolation)
    f_isolation = f_values(idx_isolation);
    xline(f_isolation, '--g', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Isolation onset: %.1f Hz', f_isolation));
end

xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Transmissibility (dB)', 'FontSize', 14, 'FontWeight', 'bold');
title('Transmissibility Curve of QZS Vibration Isolator (0–90 Hz)', 'FontSize', 16, 'FontWeight', 'bold');
xlim([0, 90]);
ylim([-60, 20]);
grid on; grid minor;
legend('Location', 'northeast', 'FontSize', 11);
set(gca, 'FontSize', 12);

%% =========  3 cases for Tr calculation   ==============
% ====================== ADD-ON: Transmissibility Curve Comparison for Different QZS Parameter Combinations ======================
% 1. Define three parameter sets (b: base radius, k: individual spring stiffness)
param_sets = {
    struct('b', 60e-3, 'k', 80,  'label', 'b=60 mm, k=80 N/m',  'color', '#E74C3C', 'linestyle', '-'),   % Red solid line
    struct('b', 80e-3, 'k', 100, 'label', 'b=80 mm, k=100 N/m', 'color', '#2ECC71', 'linestyle', '--'),  % Green dashed line
    struct('b', 90e-3, 'k', 120, 'label', 'b=90 mm, k=120 N/m', 'color', '#9B59B6', 'linestyle', ':'),   % Purple dotted line
};

% 2. Shared parameters (consistent with the original script)
x0 = 5e-3;          % Initial displacement offset (m)
zeta = 0.0311;      % Damping ratio (from experiment)
W = 2e-3;           % Base excitation amplitude (m)
g = 9.81;           % Gravitational acceleration (m/s²)
freq_range = [0, 90]; % Frequency range (Hz)
n_freq = 200;       % Number of frequency points
f_values = linspace(max(freq_range(1), 1e-6), freq_range(2), n_freq);
omega_values = 2*pi*f_values;

% 3. Create figure window
figure('Position', [100, 100, 900, 520]);
hold on; box on; grid on; grid minor;

% 4. Loop to calculate transmissibility and plot for each parameter set
for p_idx = 1:length(param_sets)
    % Extract current parameter set
    b = param_sets{p_idx}.b;
    k = param_sets{p_idx}.k;
    label = param_sets{p_idx}.label;
    color = param_sets{p_idx}.color;
    linestyle = param_sets{p_idx}.linestyle;

    % ---------------------- Static Equilibrium Analysis (reuse original logic) ----------------------
    phi_0 = asin((b - x0) / b);
    phi_array = linspace(phi_0, asin((b - 50e-3) / b), 200);
    z_disp = 2*b*(sin(phi_0) - sin(phi_array));
    F_Ver = (9*k*b/2) .* sin(phi_array) .* (cos(phi_array) - cos(phi_0)) ./ cos(phi_array);
    
    [F_max, idx_eq] = max(F_Ver);
    phi_eq = phi_array(idx_eq);
    M = F_max / g;  % Equivalent supported mass

    % ---------------------- Linearized Stiffness Calculation (reuse original logic) ----------------------
    k1 = (9*k/2) * (1 - cos(phi_0)/cos(phi_eq));
    k3 = (9*k*cos(phi_0)) / (4*b^2 * cos(phi_eq)^3);
    c = 2*zeta*sqrt(abs(k1)*M);  % Damping coefficient

    % ---------------------- Transmissibility Calculation (Harmonic Balance Method) ----------------------
    Tr_values_dB = zeros(size(f_values));
    U_prev = W;  % Initial guess for continuation method
    options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12, 'TolX', 1e-12);

    for i = 1:length(f_values)
        omega = omega_values(i);
        % HB amplitude equation (core formula)
        HB_eq = @(U) ((k1 - M*omega^2 + 0.75*k3*U.^2).^2 + (c*omega).^2).*U.^2 - (M*omega^2*W).^2;
        % Initial guess: combine continuation value and linear solution
        U_linear = abs(M*omega^2*W) / sqrt((k1 - M*omega^2)^2 + (c*omega)^2);
        U0 = max([U_prev, U_linear, 1e-9]);
        
        % Solve HB equation with fsolve
        [U_sol, ~, exitflag] = fsolve(HB_eq, U0, options);
        if exitflag <= 0 || ~isfinite(U_sol)
            % Try alternative guesses if primary solution fails
            guess_set = [U_linear, W, 0.1*W, 10*W, U0];
            U_sol = NaN;
            for gg = guess_set
                [tmp, ~, ef] = fsolve(HB_eq, max(gg,1e-9), options);
                if ef > 0 && isfinite(tmp)
                    U_sol = tmp;
                    break;
                end
            end
            if isnan(U_sol); U_sol = U0; end % Fallback to avoid crash
        end

        U = abs(U_sol);
        U_prev = U; % Update continuation guess
        k_eff = k1 + 0.75*k3*U^2; % Effective stiffness at current amplitude
        phi = atan2(c*omega, (k_eff - M*omega^2)); % Phase lag
        % Calculate displacement transmissibility (in dB)
        Tr_disp = sqrt(1 + (U/W)^2 + 2*(U/W)*cos(phi));
        Tr_values_dB(i) = 20*log10(Tr_disp);
    end

    % ---------------------- Plot transmissibility curve for current parameters ----------------------
    plot(f_values, Tr_values_dB, 'Color', color, 'LineStyle', linestyle, ...
         'LineWidth', 2, 'DisplayName', label);
    
    % Mark resonance peak for each curve
    [Tr_max, res_idx] = max(Tr_values_dB);
    f_res = f_values(res_idx);
    plot(f_res, Tr_max, 'o', 'Color', color, 'MarkerSize', 8, ...
         'MarkerFaceColor', color, 'DisplayName', sprintf('%s (Peak: %.1f Hz)', label(1:8), f_res));
end

% 5. Figure styling and annotations
yline(0, '--k', 'LineWidth', 1, 'DisplayName', '0 dB (Isolation Threshold)');  % 0 dB reference line
xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Transmissibility (dB)', 'FontSize', 14, 'FontWeight', 'bold');
title('Transmissibility Curves of QZS Vibration Isolator with Different Parameters', 'FontSize', 16, 'FontWeight', 'bold');
xlim([0, 90]);
ylim([-60, 20]);
legend('Location', 'northeast', 'FontSize', 10);
set(gca, 'FontSize', 12);
hold off;







%% ====================== 6. Plot 2: Steady-state response (acceleration vs time) ======================
% Choose a frequency for time response:
% Option A: resonance peak frequency from transmissibility curve
f_time = f_res;
% Option B: specify manually
% f_time = 10;

omega = 2*pi*f_time;

% Use computed U and phi at closest frequency point
[~, idx_closest] = min(abs(f_values - f_time));
U = U_values(idx_closest);
phi = phi_values(idx_closest);

% Time vector: show N cycles
N_cycles = 20;
T = 1/f_time;
Fs = max(5000, 2000*f_time);          % sampling rate (Hz)
dt = 1/Fs;
t = 0:dt:N_cycles*T;

% Base motion
z_b = W*cos(omega*t);
a_b = -omega^2*W*cos(omega*t);

% Relative motion
u = U*cos(omega*t - phi);
a_u = -omega^2*U*cos(omega*t - phi);

% Absolute top response
z_top = u + z_b;
a_top = a_u + a_b;

% Plot acceleration vs time
figure('Position', [120, 140, 950, 520]);
subplot(2,1,1);
plot(t, a_top, 'b-', 'LineWidth', 1.8); hold on;
plot(t, a_b, 'k--', 'LineWidth', 1.2);
xlabel('Time (s)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Acceleration (m/s^2)', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('Steady-state Acceleration Response at %.2f Hz', f_time), 'FontSize', 16, 'FontWeight', 'bold');
legend('Top acceleration a_{top}', 'Base acceleration a_{base}', 'Location', 'best');
grid on; grid minor;
set(gca, 'FontSize', 12);

subplot(2,1,2);
plot(t, a_top/g, 'b-', 'LineWidth', 1.8); hold on;
plot(t, a_b/g, 'k--', 'LineWidth', 1.2);
xlabel('Time (s)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Acceleration (g)', 'FontSize', 14, 'FontWeight', 'bold');
title('Same Response in g', 'FontSize', 14, 'FontWeight', 'bold');
legend('a_{top}/g', 'a_{base}/g', 'Location', 'best');
grid on; grid minor;
set(gca, 'FontSize', 12);

%% ====================== 6. Plot 2: Steady-state Acceleration Response (Multi-Frequency, Centered Figure) ======================
% Define frequencies to analyze (resonance peak, 10 Hz, 20 Hz)
freq_list = [f_res, 10, 20];  % f_res: resonance frequency from transmissibility curve
freq_labels = {
    sprintf('Resonance Frequency (%.2f Hz)', f_res), ...
    '10 Hz', ...
    '20 Hz'
};

% ---------------------- 核心：计算屏幕中央的窗口位置 ----------------------
screen_size = get(0, 'ScreenSize');  % 获取屏幕尺寸 [left, bottom, screen_width, screen_height]
fig_width = 950;                     % 窗口宽度
fig_height = 800;                    % 窗口高度
% 计算居中位置（MATLAB坐标原点在屏幕左下角）
fig_left = (screen_size(3) - fig_width) / 2;
fig_bottom = (screen_size(4) - fig_height) / 2;

% 创建居中显示的图窗
fig = figure('Position', [fig_left, fig_bottom, fig_width, fig_height]);  % 显式获取图窗句柄
colors = {'#2980B9', '#E67E22'};  % Blue (top) / Orange (base)
default_color = [0 0 0];          % 默认黑色（X轴/边框等）

% Loop through each frequency to calculate and plot response
for i = 1:length(freq_list)
    f_time = freq_list(i);
    omega = 2*pi*f_time;

    % Get U and phi at the closest frequency point (reuse precomputed data)
    [~, idx_closest] = min(abs(f_values - f_time));
    U = U_values(idx_closest);
    phi = phi_values(idx_closest);

    % Time vector: show 20 cycles (consistent with original logic)
    N_cycles = 20;
    T = 1/f_time;
    Fs = max(5000, 2000*f_time);  % High sampling rate to ensure smooth curves
    dt = 1/Fs;
    t = 0:dt:N_cycles*T;

    % Calculate base and top acceleration
    z_b = W*cos(omega*t);          % Base displacement
    a_b = -omega^2 * z_b;          % Base acceleration (m/s²)
    u = U*cos(omega*t - phi);      % Relative displacement
    a_u = -omega^2 * u;            % Relative acceleration
    a_top = a_u + a_b;             % Top mass acceleration (m/s²)

    % Create subplot for current frequency
    hax = subplot(3, 1, i, 'Parent', fig);  % 显式指定父图窗
    hold on; box on; grid on; grid minor;

    % Plot acceleration curves (m/s²)
    plot(t, a_top, 'Color', colors{1}, 'LineWidth', 1.8, 'DisplayName', 'Top Acceleration (a_{top})');
    plot(t, a_b, 'Color', colors{2}, 'LineStyle', '--', 'LineWidth', 1.2, 'DisplayName', 'Base Acceleration (a_{base})');

    % Set subplot labels and title
    xlabel(hax, 'Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(hax, 'Acceleration (m/s²)', 'FontSize', 12, 'FontWeight', 'bold');
    title(hax, freq_labels{i}, 'FontSize', 14, 'FontWeight', 'bold');

    % Add legend and adjust font
    legend(hax, 'Location', 'best', 'FontSize', 10);
    set(hax, 'FontSize', 11, 'Color', 'none');

    % ---------------------- 步骤1：强制左侧轴包含0点（避免0在轴范围外） ----------------------
    ylim_auto = ylim(hax);
    if ylim_auto(1) > 0
        ylim(hax, [0, ylim_auto(2)]);  % 若轴最小值>0，强制左边界为0
    elseif ylim_auto(2) < 0
        ylim(hax, [ylim_auto(1), 0]);  % 若轴最大值<0，强制右边界为0
    end
    ylim_left = ylim(hax);  % 确保包含0的左侧轴范围

    % ---------------------- 步骤2：Y轴刻度/标签颜色匹配（无版本限制） ----------------------
    set(hax, 'YColor', colors{1});  % 左侧Y轴刻度线→蓝色
    % 筛选左侧Y轴刻度标签→蓝色
    y_tick_text = findobj(hax, 'Type', 'Text', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    if ~isempty(y_tick_text), set(y_tick_text, 'Color', colors{1}); end
    % 筛选X轴刻度标签→黑色
    x_tick_text = findobj(hax, 'Type', 'Text', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    if ~isempty(x_tick_text), set(x_tick_text, 'Color', default_color); end

    % ---------------------- 步骤3：像素级精准对齐左右Y轴0点（核心修正） ----------------------
    % 临时切换轴单位为像素，便于计算
    original_units = get(hax, 'Units');
    set(hax, 'Units', 'pixels');
    ax_pos = get(hax, 'Position');  % 轴像素位置 [left, bottom, width, height]
    ax_height_px = ax_pos(4);       % 轴的像素高度
    
    % 计算0点在左侧轴的像素位置（从轴底部开始的像素数）
    zero_val = 0;
    zero_px = ax_height_px * (zero_val - ylim_left(1)) / (ylim_left(2) - ylim_left(1));
    
    % 创建右侧轴（g单位）
    ax2 = axes('Parent', fig, 'Position', ax_pos, 'YAxisLocation', 'right', 'Color', 'none', 'XTick', [], 'Units', 'pixels');
    ylabel(ax2, 'Acceleration (g)', 'FontSize', 12, 'FontWeight', 'bold');
    set(ax2, 'FontSize', 11);
    
    % 计算右侧轴的数值范围，强制0点对应相同像素位置
    ylim_right_min = zero_val - (zero_px / ax_height_px) * ( (zero_val/g) - (ylim_left(1)/g) );
    ylim_right_max = zero_val + ( (ax_height_px - zero_px) / ax_height_px ) * ( (ylim_left(2)/g) - zero_val );
    ylim(ax2, [ylim_right_min, ylim_right_max]);  % 精准对齐0点
    
    % 恢复轴单位为原始值
    set(hax, 'Units', original_units);
    set(ax2, 'Units', original_units);

    % ---------------------- 步骤4：右侧Y轴样式设置 ----------------------
    set(ax2, 'YColor', colors{2});  % 右侧Y轴刻度线→橙色
    % 筛选右侧Y轴刻度标签→橙色
    y2_tick_text = findobj(ax2, 'Type', 'Text', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    if ~isempty(y2_tick_text), set(y2_tick_text, 'Color', colors{2}); end
end

% Adjust subplot spacing for better readability
sgtitle('Steady-State Acceleration Response at Different Frequencies', 'FontSize', 16, 'FontWeight', 'bold');
% 恢复默认轴单位（避免影响后续绘图）
set(findall(fig, 'Type', 'Axes'), 'Units', 'normalized');

%%
%% ====================== 6. Plot 2: Displacement + Acceleration vs Time (Multi-Frequency, Centered Figure) ======================
% Define frequencies to analyze (resonance peak, 10 Hz, 20 Hz)
freq_list = [f_res, 10, 20];  % f_res: resonance frequency from transmissibility curve
freq_labels = {
    sprintf('Resonance Frequency (%.2f Hz)', f_res), ...
    '10 Hz', ...
    '20 Hz'
};

% ---------------------- 核心：计算屏幕中央的窗口位置 ----------------------
screen_size = get(0, 'ScreenSize');  % 获取屏幕尺寸 [left, bottom, screen_width, screen_height]
fig_width = 1200;                    % 加宽窗口适配双列子图
fig_height = 800;                    % 窗口高度
% 计算居中位置（MATLAB坐标原点在屏幕左下角）
fig_left = (screen_size(3) - fig_width) / 2;
fig_bottom = (screen_size(4) - fig_height) / 2;

% 创建居中显示的图窗
fig = figure('Position', [fig_left, fig_bottom, fig_width, fig_height]);  % 显式获取图窗句柄
colors = {'#2980B9', '#E67E22'};  % Blue (top) / Orange (base)
default_color = [0 0 0];          % 默认黑色（X轴/边框等）

% Loop through each frequency to calculate and plot response
for i = 1:length(freq_list)
    f_time = freq_list(i);
    omega = 2*pi*f_time;

    % Get U and phi at the closest frequency point (reuse precomputed data)
    [~, idx_closest] = min(abs(f_values - f_time));
    U = U_values(idx_closest);
    phi = phi_values(idx_closest);

    % Time vector: show 20 cycles (consistent with original logic)
    N_cycles = 20;
    T = 1/f_time;
    Fs = max(5000, 2000*f_time);  % High sampling rate to ensure smooth curves
    dt = 1/Fs;
    t = 0:dt:N_cycles*T;

    % ---------------------- 1. 计算位移和加速度 ----------------------
    % 位移（Displacement）：单位转换为mm（更直观）
    z_b = W * cos(omega*t) * 1000;          % Base displacement (mm)
    u = U * cos(omega*t - phi) * 1000;      % Relative displacement (mm)
    z_top = u + z_b;                        % Top mass displacement (mm)
    
    % 加速度（Acceleration）：单位m/s²
    a_b = -omega^2 * (z_b / 1000);          % Base acceleration (m/s²)
    a_u = -omega^2 * (u / 1000);            % Relative acceleration (m/s²)
    a_top = a_u + a_b;                      % Top mass acceleration (m/s²)

    % ---------------------- 2. 绘制位移-时间子图（左列） ----------------------
    hax_disp = subplot(3, 2, 2*i-1, 'Parent', fig);  % 3行2列，第i行第1列
    hold on; box on; grid on; grid minor;

    % 绘制位移曲线（Top=蓝色实线，Base=橙色虚线）
    plot(t, z_top, 'Color', colors{1}, 'LineWidth', 1.8, 'DisplayName', 'Top Displacement (z_{top})');
    plot(t, z_b, 'Color', colors{2}, 'LineStyle', '--', 'LineWidth', 1.2, 'DisplayName', 'Base Displacement (z_{base})');

    % 子图标签和标题
    xlabel(hax_disp, 'Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(hax_disp, 'Displacement (mm)', 'FontSize', 12, 'FontWeight', 'bold');
    title(hax_disp, [freq_labels{i}, ' - Displacement'], 'FontSize', 14, 'FontWeight', 'bold');

    % 图例和字体
    legend(hax_disp, 'Location', 'best', 'FontSize', 10);
    set(hax_disp, 'FontSize', 11, 'Color', 'none');

    % ---------------------- 位移子图：Y轴优化（0点对齐+颜色匹配） ----------------------
    % 强制包含0点
    ylim_auto = ylim(hax_disp);
    if ylim_auto(1) > 0, ylim(hax_disp, [0, ylim_auto(2)]);
    elseif ylim_auto(2) < 0, ylim(hax_disp, [ylim_auto(1), 0]); end
    ylim_left = ylim(hax_disp);

    % Y轴刻度线+标签颜色匹配
    set(hax_disp, 'YColor', colors{1});
    y_tick_text = findobj(hax_disp, 'Type', 'Text', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    if ~isempty(y_tick_text), set(y_tick_text, 'Color', colors{1}); end
    x_tick_text = findobj(hax_disp, 'Type', 'Text', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    if ~isempty(x_tick_text), set(x_tick_text, 'Color', default_color); end

    % ---------------------- 3. 绘制加速度-时间子图（右列） ----------------------
    hax_acc = subplot(3, 2, 2*i, 'Parent', fig);  % 3行2列，第i行第2列
    hold on; box on; grid on; grid minor;

    % 绘制加速度曲线（复用原有逻辑）
    plot(t, a_top, 'Color', colors{1}, 'LineWidth', 1.8, 'DisplayName', 'Top Acceleration (a_{top})');
    plot(t, a_b, 'Color', colors{2}, 'LineStyle', '--', 'LineWidth', 1.2, 'DisplayName', 'Base Acceleration (a_{base})');

    % 子图标签和标题
    xlabel(hax_acc, 'Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(hax_acc, 'Acceleration (m/s²)', 'FontSize', 12, 'FontWeight', 'bold');
    title(hax_acc, [freq_labels{i}, ' - Acceleration'], 'FontSize', 14, 'FontWeight', 'bold');

    % 图例和字体
    legend(hax_acc, 'Location', 'best', 'FontSize', 10);
    set(hax_acc, 'FontSize', 11, 'Color', 'none');

    % ---------------------- 加速度子图：原有优化逻辑（0点对齐+颜色匹配） ----------------------
    % 步骤1：强制左侧轴包含0点
    ylim_auto = ylim(hax_acc);
    if ylim_auto(1) > 0
        ylim(hax_acc, [0, ylim_auto(2)]);
    elseif ylim_auto(2) < 0
        ylim(hax_acc, [ylim_auto(1), 0]);
    end
    ylim_left = ylim(hax_acc);

    % 步骤2：Y轴刻度/标签颜色匹配
    set(hax_acc, 'YColor', colors{1});
    y_tick_text = findobj(hax_acc, 'Type', 'Text', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    if ~isempty(y_tick_text), set(y_tick_text, 'Color', colors{1}); end
    x_tick_text = findobj(hax_acc, 'Type', 'Text', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    if ~isempty(x_tick_text), set(x_tick_text, 'Color', default_color); end

    % 步骤3：像素级精准对齐左右Y轴0点（加速度双Y轴）
    original_units = get(hax_acc, 'Units');
    set(hax_acc, 'Units', 'pixels');
    ax_pos = get(hax_acc, 'Position');
    ax_height_px = ax_pos(4);
    
    zero_val = 0;
    zero_px = ax_height_px * (zero_val - ylim_left(1)) / (ylim_left(2) - ylim_left(1));
    
    % 创建右侧g单位轴
    ax2 = axes('Parent', fig, 'Position', ax_pos, 'YAxisLocation', 'right', 'Color', 'none', 'XTick', [], 'Units', 'pixels');
    ylabel(ax2, 'Acceleration (g)', 'FontSize', 12, 'FontWeight', 'bold');
    set(ax2, 'FontSize', 11);
    
    ylim_right_min = zero_val - (zero_px / ax_height_px) * ( (zero_val/g) - (ylim_left(1)/g) );
    ylim_right_max = zero_val + ( (ax_height_px - zero_px) / ax_height_px ) * ( (ylim_left(2)/g) - zero_val );
    ylim(ax2, [ylim_right_min, ylim_right_max]);
    
    % 恢复单位
    set(hax_acc, 'Units', original_units);
    set(ax2, 'Units', original_units);

    % 步骤4：右侧Y轴样式
    set(ax2, 'YColor', colors{2});
    y2_tick_text = findobj(ax2, 'Type', 'Text', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    if ~isempty(y2_tick_text), set(y2_tick_text, 'Color', colors{2}); end
end

% ---------------------- 全局设置 ----------------------
% 总标题
sgtitle('Steady-State Displacement & Acceleration Response at Different Frequencies', 'FontSize', 16, 'FontWeight', 'bold');
% 调整子图间距
set(gcf, 'DefaultAxesPosition', get(gcf, 'DefaultAxesPosition'));
% 恢复默认轴单位
set(findall(fig, 'Type', 'Axes'), 'Units', 'normalized');

%% ====================== 7. Plot 3: Theory vs Experimental Comparison (1) ======================
% 请确保CSV文件路径正确，第一列为频率(Hz)，第二列为传递率(dB)
% exp_data_path = 'Vib_data_2026\Data_0121.csv';  % 实验数据路径
exp_data_path = 'Vib_data_2026\Data_0129-1.csv';  % 实验数据路径

try
    exp_data = csvread(exp_data_path, 0, 0);  % 跳过表头（如果有），从第2行开始读取
    freq_exp = exp_data(:, 1)-1;  % 原始实验频率
    Tr_exp = exp_data(:, 2);    % 原始实验传递率
    fprintf('实验数据读取成功，共 %d 个原始数据点\n', length(freq_exp));
catch ME
    error('实验数据读取失败：%s，请检查文件路径和格式', ME.message);
end

% ---------------------- 新增：实验频率自动取整+去重（核心功能） ----------------------
% 步骤1：将实验频率四舍五入到临近整数
freq_exp_round = round(freq_exp);

% 步骤2：处理重复整数频率（相同整数频率仅保留一个，取对应Tr的平均值）
% 获取唯一的整数频率，并按升序排列
[unique_freq, ~, idx] = unique(freq_exp_round, 'sorted');
% 对每个唯一整数频率，计算对应Tr的平均值（保留实验数据统计特征）
unique_Tr = arrayfun(@(x) mean(Tr_exp(idx==x)), 1:length(unique_freq));

% 步骤3：更新实验数据为“整数频率+平均Tr”（无重复，相邻数值不同）
freq_exp = unique_freq;    % 去重后的整数频率
Tr_exp = unique_Tr;        % 对应平均传递率
fprintf('频率取整+去重后，剩余 %d 个有效数据点\n', length(freq_exp));

% ---------------------- 数据预处理（补充） ----------------------
% 过滤超出理论频率范围的数据点（0-90Hz）
valid_idx = freq_exp >= 0 & freq_exp <= 90;
freq_exp = freq_exp(valid_idx);
Tr_exp = Tr_exp(valid_idx);
fprintf('过滤0-90Hz范围后，剩余 %d 个数据点\n', length(freq_exp));

% ---------------------- 创建对比图 ----------------------
figure('Position', [150, 150, 950, 550]);
hold on; box on; grid on; grid minor;

% 1. 绘制理论传递率曲线（实线）
plot(f_values, Tr_values_dB, 'b-', 'LineWidth', 2.2, 'DisplayName', 'Theory (HB, disp. Tr)');

% 2. 绘制处理后的实验数据（点连线）
plot(freq_exp, Tr_exp, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'r', ...
     'DisplayName', 'Experimental Data');

% 3. 添加参考线和标注
yline(0, '--k', 'LineWidth', 1, 'DisplayName', '0 dB');  % 0 dB参考线

% 标注理论共振峰
[Tr_max_dB, idx_res] = max(Tr_values_dB);
f_res = f_values(idx_res);
plot(f_res, Tr_max_dB, 'bs', 'MarkerSize', 8, 'MarkerFaceColor', 'b', ...
     'DisplayName', sprintf('Theory Peak: %.1f Hz', f_res));

% 标注理论隔振起始频率
idx_isolation = find(Tr_values_dB < 0, 1, 'first');
if ~isempty(idx_isolation)
    f_isolation = f_values(idx_isolation);
    xline(f_isolation, '--g', 'LineWidth', 1.5, ...
          'DisplayName', sprintf('Isolation Onset: %.1f Hz', f_isolation));
end

% ---------------------- 图形样式设置 ----------------------
xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Transmissibility (dB)', 'FontSize', 14, 'FontWeight', 'bold');
title('Experimental vs Theoretical Transmissibility Comparison of QZS Vibration Isolator', 'FontSize', 16, 'FontWeight', 'bold');
xlim([0, 90]);
ylim([-60, 20]);
legend('Location', 'northeast', 'FontSize', 11);
set(gca, 'FontSize', 12);

hold off;


 

%% ====================== 8. (Optional) Static plots kept for reference ======================
% Static force-displacement curve
figure('Position', [200, 200, 850, 520]);
hold on; box on;
plot(z_disp*1000, F_Ver, 'b-', 'LineWidth', 2, 'DisplayName', 'Force-Displacement');
plot(z_eq*1000, F_max, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Equilibrium');
xlabel('Displacement (mm)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Restoring Force (N)', 'FontSize', 14, 'FontWeight', 'bold');
title('Static Force-Displacement Characteristic', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
grid on; grid minor;
set(gca, 'FontSize', 12);

% Local stiffness
dz = diff(z_disp);
dF = diff(F_Ver);
k_local = dF ./ dz;
z_mid = (z_disp(1:end-1) + z_disp(2:end)) / 2;

figure('Position', [250, 250, 850, 420]);
plot(z_mid*1000, k_local, 'b-', 'LineWidth', 1.5); hold on;
QZS_threshold = 0.1*k;
yline(0, '--k', 'LineWidth', 1);
yline(QZS_threshold, ':r', 'LineWidth', 1);
yline(-QZS_threshold, ':r', 'LineWidth', 1);
xlabel('Displacement (mm)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Local Stiffness (N/m)', 'FontSize', 14, 'FontWeight', 'bold');
title('Local Stiffness vs Displacement', 'FontSize', 16, 'FontWeight', 'bold');
grid on; grid minor;
set(gca, 'FontSize', 12);

%% 3 Cases: theorey + experiment


%% ====================== 9. Output Summary ======================
fprintf('\n========== ANALYSIS SUMMARY ==========\n');
fprintf('Geometric Parameters:\n');
fprintf('  b = %.1f mm, l0 = %.1f mm\n', b*1000, l0*1000);
fprintf('  Initial angle φ_0 = %.2f°\n', phi_0*180/pi);
fprintf('\nDynamic Properties:\n');
fprintf('  Supported mass M = %.3f kg\n', M);
fprintf('  Linear stiffness k1 = %.2f N/m\n', k1);
fprintf('  Nonlinear stiffness k3 = %.2e N/m^3\n', k3);
fprintf('  Natural frequency = %.2f Hz\n', f_n);
fprintf('  Damping ratio ζ = %.4f\n', zeta);
fprintf('\nIsolation Performance (Disp. Tr):\n');
fprintf('  Peak transmissibility = %.2f dB at %.1f Hz\n', Tr_max_dB, f_res);
if ~isempty(idx_isolation)
    fprintf('  Isolation onset frequency = %.1f Hz\n', f_isolation);
end
fprintf('======================================\n');