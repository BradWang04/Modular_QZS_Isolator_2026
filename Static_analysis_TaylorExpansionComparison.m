%% QZS Vibration Isolator - Enhanced Force-Displacement Approximation
% 对比原始力-位移曲线与多种高阶近似方法
% 方法1：高阶泰勒展开（数值微分）
% 方法2：多项式最小二乘拟合（推荐）
clear all; clc; close all;
format long g;

%% ====================== 1. 系统基础参数 ======================
b = 80e-3;        % 连杆长度 (m)
k = 300;          % 单根弹簧刚度 (N/m)
x0 = 5e-3;        % 初始偏移量 (m)
g = 9.81;         % 重力加速度 (m/s^2)
phi_lim = asin((b - 50e-3) / b);
n_phi = 500;      % 采样点数

%% ====================== 2. 计算原始力-位移曲线 ======================
phi0 = asin((b - x0) / b);
phi_array = linspace(phi0, phi_lim, n_phi);
z_disp = 2*b*(sin(phi0) - sin(phi_array));
F_Ver = (9*k*b/2) .* sin(phi_array) .* (cos(phi_array) - cos(phi0)) ./ cos(phi_array);

% 找到静平衡位置
[F_eq, idx_eq] = max(F_Ver);
phi_eq = phi_array(idx_eq);
z_eq = z_disp(idx_eq);
M_eq = F_eq / g;

fprintf('====================== 静平衡关键参数 ======================\n');
fprintf('平衡位置位移 z_eq = %.4f mm\n', z_eq*1000);
fprintf('平衡位置承载力 F_eq = %.4f N\n', F_eq);
fprintf('平衡位置支撑质量 M_eq = %.4f kg\n', M_eq);
fprintf('============================================================\n\n');

%% ====================== 3. 相对位移与恢复力 ======================
u = z_disp - z_eq;  % 相对平衡位置的位移 (m)
F_restore_original = F_Ver - F_eq;  % 恢复力 (N)

%% ====================== 4. 方法1：高阶泰勒展开（数值微分） ======================
fprintf('====================== 方法1：数值泰勒展开 ======================\n');

% 使用中心差分法计算各阶导数（在平衡位置处）
% 选择合适的步长（相对于位移范围）
h = (max(z_disp) - min(z_disp)) / 1000;  % 步长

% 定义力函数（插值原始数据）
F_interp = @(z) interp1(z_disp, F_Ver, z, 'spline', 'extrap');

% 计算1-7阶导数（使用中心差分）
taylor_order = 7;  % 泰勒展开阶数
k_taylor = zeros(1, taylor_order);

% 一阶导数（线性刚度）
k_taylor(1) = (F_interp(z_eq + h) - F_interp(z_eq - h)) / (2*h);

% 二阶导数
k_taylor(2) = (F_interp(z_eq + h) - 2*F_interp(z_eq) + F_interp(z_eq - h)) / h^2;

% 三阶导数
k_taylor(3) = (F_interp(z_eq + 2*h) - 2*F_interp(z_eq + h) + 2*F_interp(z_eq - h) - F_interp(z_eq - 2*h)) / (2*h^3);

% 四阶导数
k_taylor(4) = (F_interp(z_eq + 2*h) - 4*F_interp(z_eq + h) + 6*F_interp(z_eq) - 4*F_interp(z_eq - h) + F_interp(z_eq - 2*h)) / h^4;

% 五阶及以上使用数值微分（简化处理）
for n = 5:taylor_order
    % 使用有限差分近似高阶导数
    delta = h *2;
    z_sample = z_eq + linspace(-delta, delta, 2*n+1);
    F_sample = arrayfun(F_interp, z_sample);
    % 使用多项式拟合后求导
    p = polyfit(z_sample - z_eq, F_sample - F_eq, n);
    k_taylor(n) = p(end-n+1) * factorial(n);
end

% 构建泰勒展开恢复力
F_restore_taylor = zeros(size(u));
for n = 1:taylor_order
    F_restore_taylor = F_restore_taylor + k_taylor(n) .* (u.^n) / factorial(n);
end

% 输出泰勒系数
fprintf('泰勒展开系数（阶数=%d）：\n', taylor_order);
for n = 1:taylor_order
    fprintf('  k%d / %d! = %.6e N/m^%d\n', n, n, k_taylor(n)/factorial(n), n);
end

% 计算泰勒展开误差
rmse_taylor = sqrt(mean((F_restore_original - F_restore_taylor).^2));
r2_taylor = 1 - sum((F_restore_original - F_restore_taylor).^2) / sum((F_restore_original - mean(F_restore_original)).^2);
fprintf('泰勒展开RMSE = %.6f N, R² = %.6f\n', rmse_taylor, r2_taylor);
fprintf('============================================================\n\n');

%% ====================== 5. 方法2：多项式最小二乘拟合（推荐） ======================
fprintf('====================== 方法2：多项式拟合（推荐） ======================\n');

% 尝试不同阶数的多项式拟合
poly_orders = [3, 5, 7, 9];
colors = lines(length(poly_orders));
rmse_poly = zeros(size(poly_orders));
r2_poly = zeros(size(poly_orders));
poly_coeffs = cell(size(poly_orders));

for i = 1:length(poly_orders)
    order = poly_orders(i);
    % polyfit拟合（注意：polyfit返回降序系数，p(1)*x^n + ... + p(n+1)）
    p = polyfit(u, F_restore_original, order);
    poly_coeffs{i} = p;
    
    % 计算拟合值
    F_restore_poly = polyval(p, u);
    
    % 计算误差
    rmse_poly(i) = sqrt(mean((F_restore_original - F_restore_poly).^2));
    r2_poly(i) = 1 - sum((F_restore_original - F_restore_poly).^2) / sum((F_restore_original - mean(F_restore_original)).^2);
    
    fprintf('%d阶多项式拟合：\n', order);
    fprintf('  RMSE = %.6f N, R² = %.8f\n', rmse_poly(i), r2_poly(i));
    fprintf('  系数（降序）：');
    fprintf('%.6e ', p);
    fprintf('\n');
end

% 选择最佳拟合阶数（R²最高且RMSE最小）
[~, best_idx] = max(r2_poly);
best_order = poly_orders(best_idx);
best_coeffs = poly_coeffs{best_idx};
F_restore_best = polyval(best_coeffs, u);

fprintf('\n推荐使用：%d阶多项式拟合\n', best_order);
fprintf('  RMSE = %.6f N, R² = %.8f\n', rmse_poly(best_idx), r2_poly(best_idx));
fprintf('============================================================\n\n');

%% ====================== 6. 提取物理刚度系数（从最佳拟合） ======================
fprintf('====================== 物理刚度系数（从%d阶拟合提取） ======================\n', best_order);
% polyfit返回降序系数：p(1)*u^n + p(2)*u^(n-1) + ... + p(n+1)
% 转换为升序：c1*u + c2*u^2 + c3*u^3 + ...
c = fliplr(best_coeffs);% 转为升序
c = c(2:end);  % 去掉常数项（应该接近0）

fprintf('恢复力表达式：F_restore = ');
for n = 1:length(c)
    if n == 1
        fprintf('%.4e*u', c(n));
    else
        fprintf(' + %.4e*u^%d', c(n), n);
    end
end
fprintf(' (u单位:m, F单位:N)\n\n');

% 关键刚度系数
k1_fit = c(1);  % 线性刚度
if length(c) >= 3
    k3_fit = c(3);  % 三次非线性刚度
    fprintf('线性刚度 k1 = %.4f N/m\n', k1_fit);
    fprintf('三次非线性刚度 k3 = %.4e N/m³\n', k3_fit);
end
fprintf('============================================================\n\n');

%% ====================== 7. 绘制专业对比图 ======================
figure('Position', [100, 100, 1400, 900], 'Color', 'w');

% --- 子图1：绝对位移-总力对比 ---
subplot(2, 2, 1); hold on; box on; grid on;
plot(z_disp*1000, F_Ver, 'k-', 'LineWidth', 3, 'DisplayName', 'Original');
plot(z_disp*1000, F_eq + F_restore_best, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Poly Fit (order=%d)', best_order));
plot(z_eq*1000, F_eq, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'y', 'DisplayName', 'Equilibrium');
xlabel('Vertical Displacement z (mm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Vertical Force F (N)', 'FontSize', 12, 'FontWeight', 'bold');
title('(a) Absolute F-z Curve', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11, 'LineWidth', 1.2);

% --- 子图2：相对位移-恢复力对比（全范围） ---
subplot(2, 2, 2); hold on; box on; grid on;
plot(u*1000, F_restore_original, 'k-', 'LineWidth', 3, 'DisplayName', 'Original');
plot(u*1000, F_restore_taylor, 'b--', 'LineWidth', 2, 'DisplayName', sprintf('Taylor (order=%d)', taylor_order));
plot(u*1000, F_restore_best, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Poly Fit (order=%d)', best_order));
plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
xlabel('Relative Displacement u (mm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Restoring Force F_{restore} (N)', 'FontSize', 12, 'FontWeight', 'bold');
title('(b) Restoring Force (Full Range)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11, 'LineWidth', 1.2);

% --- 子图3：小位移范围细节（±2mm） ---
subplot(2, 2, 3); hold on; box on; grid on;
u_range = 2;  % mm
idx_small = abs(u*1000) <= u_range;
plot(u(idx_small)*1000, F_restore_original(idx_small), 'k-', 'LineWidth', 3, 'DisplayName', 'Original');
plot(u(idx_small)*1000, F_restore_taylor(idx_small), 'b--', 'LineWidth', 2, 'DisplayName', 'Taylor');
plot(u(idx_small)*1000, F_restore_best(idx_small), 'r--', 'LineWidth', 2, 'DisplayName', 'Poly Fit');
plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
xlabel('Relative Displacement u (mm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Restoring Force F_{restore} (N)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('(c) Small Displacement Detail (±%dmm)', u_range), 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11, 'LineWidth', 1.2);
xlim([-u_range, u_range]);

% --- 子图4：拟合误差对比 ---
subplot(2, 2, 4); hold on; box on; grid on;
error_taylor = F_restore_original - F_restore_taylor;
error_poly = F_restore_original - F_restore_best;
plot(u*1000, error_taylor, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Taylor Error (RMSE=%.3fN)', rmse_taylor));
plot(u*1000, error_poly, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('Poly Fit Error (RMSE=%.3fN)', rmse_poly(best_idx)));
yline(0, 'k--', 'LineWidth', 1);
xlabel('Relative Displacement u (mm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Approximation Error (N)', 'FontSize', 12, 'FontWeight', 'bold');
title('(d) Fitting Error Comparison', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11, 'LineWidth', 1.2);

sgtitle('QZS Isolator: Enhanced Force-Displacement Approximation', 'FontSize', 16, 'FontWeight', 'bold');

%% ====================== 8. 绘制多阶多项式拟合对比 ======================
figure('Position', [150, 150, 1200, 500], 'Color', 'w');

subplot(1, 2, 1); hold on; box on; grid on;
plot(u*1000, F_restore_original, 'k-', 'LineWidth', 3, 'DisplayName', 'Original');
for i = 1:length(poly_orders)
    F_poly_i = polyval(poly_coeffs{i}, u);
    plot(u*1000, F_poly_i, '--', 'LineWidth', 2, 'Color', colors(i,:), 'DisplayName', sprintf('Order %d (R²=%.6f)', poly_orders(i), r2_poly(i)));
end
xlabel('Relative Displacement u (mm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Restoring Force F_{restore} (N)', 'FontSize', 12, 'FontWeight', 'bold');
title('Polynomial Fitting: Different Orders', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11, 'LineWidth', 1.2);

subplot(1, 2, 2); hold on; box on; grid on;
yyaxis left
bar(poly_orders, rmse_poly, 'FaceColor', [0.3, 0.6, 0.9], 'EdgeColor', 'k', 'LineWidth', 1.5);
ylabel('RMSE (N)', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'YColor', 'k');

yyaxis right
plot(poly_orders, r2_poly, 'ro-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
ylabel('R² (Goodness of Fit)', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'YColor', 'r');
ylim([min(r2_poly)*0.999, 1.001]);

xlabel('Polynomial Order', 'FontSize', 12, 'FontWeight', 'bold');
title('Fitting Accuracy vs Polynomial Order', 'FontSize', 13, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11, 'LineWidth', 1.2);

%% ====================== 9. 保存结果 ======================
saveas(gcf, 'QZS_Polynomial_Fitting_Comparison.png');
save('QZS_Enhanced_Approximation.mat', 'best_order', 'best_coeffs', 'k1_fit', 'k3_fit', ...
    'rmse_poly', 'r2_poly', 'z_eq', 'F_eq', 'M_eq');

fprintf('\n====================== 计算完成 ======================\n');
fprintf('对比图已保存\n');
fprintf('推荐使用 %d 阶多项式拟合，R² = %.8f\n', best_order, r2_poly(best_idx));
fprintf('====================================================\n');