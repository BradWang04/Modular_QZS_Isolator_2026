%% 10/26 New Formula

clear all; clc;
%% QZS Figure - 50*80 configurations * 3
close all;

% 设置下采样因子
downsample_factor = 10; % 例如，减少到原来的 1/10

% Data_0429_1 = csvread("Static_Testing/0429_QZS_1_50mm.csv",1);
Data_0429_1 = csvread("Static_Testing/Stat_QZS_1.csv",1);
deform_1    = Data_0429_1(:,1);
force_1     = Data_0429_1(:,2);
% 对变形数据进行下采样
% deform_1_downsamplved = downsample(deform_1, downsample_factor);
% % 对力数据进行下采样
% force_1_downsampled = downsample(force_1, downsample_factor);

Data_0429_2 = csvread("Static_Testing/Stat_QZS_2.csv",1);
deform_2    = Data_0429_2(:,1);
% deform_2_downsampled = downsample(deform_2, downsample_factor);
force_2     = Data_0429_2(:,2);
% force_2_downsampled = downsample(force_2, downsample_factor);


Data_0429_3 = csvread("Static_Testing/Stat_QZS_3.csv",1);
deform_3    = Data_0429_3(:,1);
% deform_3_downsampled = downsample(deform_3, downsample_factor);
force_3     = Data_0429_3(:,2);
% force_3_downsampled = downsample(force_3, downsample_factor);


Data_0429_4 = csvread("Static_Testing/Stat_QZS_4.csv",1);
deform_4    = Data_0429_4(:,1);
% deform_4_downsampled = downsample(deform_4, downsample_factor);
force_4     = Data_0429_4(:,2);
% force_4_downsampled = downsample(force_4, downsample_factor);

% ------- Spring Structure Maths Derivation, static -----
% 第一段
a_1 = 50; %mm
b_1 = 80; %mm
l_0_1 = 20;
k_1 = (140/1E3*9.8)/(80); % N/mm
z_t_1 = [0:1:60]; % Displacement
l_1 = ( l_0_1^2 - 9/4*z_t_1.^2 + 3*z_t_1.*(3*b_1^2 - l_0_1^2)^0.5 ).^0.5; % Eq. 3-6
F_ver_1 = k_1*l_1.*(l_1 - l_0_1).*sqrt(3*b_1^2 - l_1.^2)/b_1^2;

% 第二段
a_2 = 50; %mm
b_2 = 80; %mm
l_0_2 = 30;
k_2 = (460/1E3*9.8)/(80); % N/mm
z_t_2 = [0:1:60]; % Displacement
l_2 = ( l_0_2^2 - 9/4*z_t_2.^2 + 3*z_t_2.*(3*b_2^2 - l_0_2^2)^0.5 ).^0.5; % Eq. 3-6
F_ver_2 = k_2*l_2.*(l_2 - l_0_2).*sqrt(3*b_2^2 - l_2.^2)/b_2^2;

% 第三段
a_3 = 50; %mm
b_3 = 80; %mm
l_0_3 = 30;
k_3 = (950/1E3*9.8)/(80); % N/mm
z_t_3 = [0:1:60]; % Displacement
l_3 = ( l_0_3^2 - 9/4*z_t_3.^2 + 3*z_t_3.*(3*b_3^2 - l_0_3^2)^0.5 ).^0.5; % Eq. 3-6
F_ver_3 = k_3*l_3.*(l_3 - l_0_3).*sqrt(3*b_3^2 - l_3.^2)/b_3^2;


%------ Plotting
color_full = ["#2e974e" ,  "#e25508" , "#2e7ebb" ,"#5f5f5f" , "#7262ac" , "#d92523"];
colors_2 = ["#b8e3b2","#fdc38d" ,"#b7d4ea", "#cecece", "#cfcfe5",   "#fcab8f", "#adadae", "#c2e6f7", "#e7c6db" , "#b3c7c8"]; 
rgb_array = cellfun(@(hex) sscanf(hex(2:end), '%2x%2x%2x', [1 3]), colors_2, 'UniformOutput', false);
rgb_array = cell2mat(rgb_array);
% disp(rgb_array)

figure()
hold on
plot(deform_4,force_4,'^','Color',colors_2(3)) % 04x05x50
plot(z_t_1,F_ver_1,'Color',color_full(3),'LineWidth',1.5)

plot(deform_1,force_1,'o','Color',colors_2(1)) % 05x04x35
plot(z_t_2,F_ver_2,'Color',color_full(1),'LineWidth',1.5)

plot(deform_3,force_3,'s','Color',colors_2(2)) % 05x05x35
plot(z_t_3,F_ver_3,'Color',color_full(2),'LineWidth',1.5)

% legend('05x04x35','05x05x35','04x05x50','Maths',fontsize=12);
legend('Spring 1 - Test','Spring 1 - Model','Spring 2 - Test','Spring 2 - Model','Spring 3 - Test','Spring 3 - Model')
axis([0 60 0 18])
grid on
box on;

% title('Force-Deflection Curve for Horizontal Spring QZS Structure',fontsize=12)
xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')

%% QZS characterization

%% QZS Figure - 50*80 configurations * 3
close all;

% 设置下采样因子
downsample_factor = 10; % 例如，减少到原来的 1/10

% Data_0429_1 = csvread("Static_Testing/0429_QZS_1_50mm.csv",1);
Data_0429_1 = csvread("Static_Testing/Stat_QZS_1.csv",1);
deform_1    = Data_0429_1(:,1);
force_1     = Data_0429_1(:,2);

Data_0429_2 = csvread("Static_Testing/Stat_QZS_2.csv",1);
deform_2    = Data_0429_2(:,1);
force_2     = Data_0429_2(:,2);

Data_0429_3 = csvread("Static_Testing/Stat_QZS_3.csv",1);
deform_3    = Data_0429_3(:,1);
force_3     = Data_0429_3(:,2);

Data_0429_4 = csvread("Static_Testing/Stat_QZS_4.csv",1);
deform_4    = Data_0429_4(:,1);
force_4     = Data_0429_4(:,2);

% ------- Spring Structure Maths Derivation, static -----
% 第一段
a_1 = 50; %mm
b_1 = 80; %mm
l_0_1 = 20;
k_1 = (140/1E3*9.8)/(80); % N/mm
z_t_1 = [0:0.1:60]; % 使用更细的步长以便更准确地计算导数
l_1 = ( l_0_1^2 - 9/4*z_t_1.^2 + 3*z_t_1.*(3*b_1^2 - l_0_1^2)^0.5 ).^0.5; % Eq. 3-6
F_ver_1 = k_1*l_1.*(l_1 - l_0_1).*sqrt(3*b_1^2 - l_1.^2)/b_1^2;

% 第二段
a_2 = 50; %mm
b_2 = 80; %mm
l_0_2 = 30;
k_2 = (460/1E3*9.8)/(80); % N/mm
z_t_2 = [0:0.1:60]; % 使用更细的步长
l_2 = ( l_0_2^2 - 9/4*z_t_2.^2 + 3*z_t_2.*(3*b_2^2 - l_0_2^2)^0.5 ).^0.5; % Eq. 3-6
F_ver_2 = k_2*l_2.*(l_2 - l_0_2).*sqrt(3*b_2^2 - l_2.^2)/b_2^2;

% 第三段
a_3 = 50; %mm
b_3 = 80; %mm
l_0_3 = 30;
k_3 = (950/1E3*9.8)/(80); % N/mm
z_t_3 = [0:0.1:60]; % 使用更细的步长
l_3 = ( l_0_3^2 - 9/4*z_t_3.^2 + 3*z_t_3.*(3*b_3^2 - l_0_3^2)^0.5 ).^0.5; % Eq. 3-6
F_ver_3 = k_3*l_3.*(l_3 - l_0_3).*sqrt(3*b_3^2 - l_3.^2)/b_3^2;

% 计算刚度（力对位移的导数）
stiffness_1 = gradient(F_ver_1) ./ gradient(z_t_1);
stiffness_2 = gradient(F_ver_2) ./ gradient(z_t_2);
stiffness_3 = gradient(F_ver_3) ./ gradient(z_t_3);

% 找到QZS区间（刚度绝对值小于0.005 N/mm）
qzs_threshold = 0.035; % N/mm
qzs_threshold_1 = 0.1 * k_1;
qzs_threshold_2 = 0.1 * k_2;
qzs_threshold_3 = 0.1 * k_3;

qzs_indices_1 = find(abs(stiffness_1) < qzs_threshold_1);
qzs_indices_2 = find(abs(stiffness_2) < qzs_threshold_2);
qzs_indices_3 = find(abs(stiffness_3) < qzs_threshold_3);


[qzs_start_1, qzs_end_1] = find_continuous_regions(qzs_indices_1);
[qzs_start_2, qzs_end_2] = find_continuous_regions(qzs_indices_2);
[qzs_start_3, qzs_end_3] = find_continuous_regions(qzs_indices_3);

%------ Plotting
color_full = ["#2e974e" ,  "#e25508" , "#2e7ebb" ,"#5f5f5f" , "#7262ac" , "#d92523"];
colors_2 = ["#b8e3b2","#fdc38d" ,"#b7d4ea", "#cecece", "#cfcfe5",   "#fcab8f", "#adadae", "#c2e6f7", "#e7c6db" , "#b3c7c8"]; 
rgb_array = cellfun(@(hex) sscanf(hex(2:end), '%2x%2x%2x', [1 3]), colors_2, 'UniformOutput', false);
rgb_array = cell2mat(rgb_array);

figure()
hold on

% 绘制QZS区间的半透明矩形
% 对于Spring 1
for i = 1:length(qzs_start_1)
    x_start = z_t_1(qzs_start_1(i));
    x_end = z_t_1(qzs_end_1(i));
    y_start = F_ver_1(qzs_start_1(i));
    y_end = max(F_ver_1);
    rectangle('Position', [x_start, 0.8*y_start, x_end-x_start, 6*(y_end-y_start)], ...
              'FaceColor', [0.5, 0.7, 1, 0.2], 'EdgeColor', 'k');
end

% 对于Spring 2
for i = 1:length(qzs_start_2)
    x_start = z_t_2(qzs_start_2(i));
    x_end = z_t_2(qzs_end_2(i));
    y_start = F_ver_2(qzs_start_2(i));
    y_end = max(F_ver_2);
    rectangle('Position', [x_start, 0.95*y_start, x_end-x_start, 8*(y_end-y_start)], ...
              'FaceColor', [0.5, 1, 0.5, 0.2], 'EdgeColor', 'k');
end

% 对于Spring 3
for i = 1:length(qzs_start_3)
    x_start = z_t_3(qzs_start_3(i));
    x_end = z_t_3(qzs_end_3(i));
    y_start = F_ver_3(qzs_start_3(i));
    y_end = max(F_ver_3);
    rectangle('Position', [x_start, 0.95*y_start, x_end-x_start, 25*(y_end-y_start)], ...
              'FaceColor', [1, 0.7, 0.5, 0.2], 'EdgeColor', 'k');
end

% 绘制原始数据和曲线
plot(deform_4,force_4,'^','Color',colors_2(3)) % 04x05x50
plot(z_t_1,F_ver_1,'Color',color_full(3),'LineWidth',1.5)

plot(deform_1,force_1,'o','Color',colors_2(1)) % 05x04x35
plot(z_t_2,F_ver_2,'Color',color_full(1),'LineWidth',1.5)

plot(deform_3,force_3,'s','Color',colors_2(2)) % 05x05x35
plot(z_t_3,F_ver_3,'Color',color_full(2),'LineWidth',1.5)

legend('Spring 1 - Test','Spring 1 - Model','Spring 2 - Test','Spring 2 - Model','Spring 3 - Test','Spring 3 - Model')
axis([0 60 0 18])
grid on
box on;

xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')

% 添加子图显示刚度曲线
figure()
subplot(3,1,1)
plot(z_t_1, abs(stiffness_1), 'Color', color_full(3), 'LineWidth', 1.5)
hold on
plot([0 60], [qzs_threshold qzs_threshold], 'r--', 'LineWidth', 1)
% ylabel('|Stiffness| [N/mm]')
title('Spring 1 Stiffness','FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
grid on
axis([0 60 0 max(abs(stiffness_1))*1.1])

subplot(3,1,2)
plot(z_t_2, abs(stiffness_2), 'Color', color_full(1), 'LineWidth', 1.5)
hold on
plot([0 60], [qzs_threshold qzs_threshold], 'r--', 'LineWidth', 1)
ylabel('|Stiffness| [N/mm]','FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
title('Spring 2 Stiffness','FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
grid on
axis([0 60 0 max(abs(stiffness_2))*1.1])

subplot(3,1,3)
plot(z_t_3, abs(stiffness_3), 'Color', color_full(2), 'LineWidth', 1.5)
hold on
plot([0 60], [qzs_threshold qzs_threshold], 'r--', 'LineWidth', 1)
xlabel('Displacement [mm]','FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
% ylabel('|Stiffness| [N/mm]','FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
title('Spring 3 Stiffness','FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
grid on
axis([0 60 0 max(abs(stiffness_3))*1.1])

% 输出QZS区间信息
fprintf('\nQZS Regions (|stiffness| < %.3f N/mm):\n', qzs_threshold);
fprintf('\nSpring 1:\n');
for i = 1:length(qzs_start_1)
    fprintf('  Region %d: %.2f mm - %.2f mm\n', i, z_t_1(qzs_start_1(i)), z_t_1(qzs_end_1(i)));
end
fprintf('\nSpring 2:\n');
for i = 1:length(qzs_start_2)
    fprintf('  Region %d: %.2f mm - %.2f mm\n', i, z_t_2(qzs_start_2(i)), z_t_2(qzs_end_2(i)));
end
fprintf('\nSpring 3:\n');
for i = 1:length(qzs_start_3)
    fprintf('  Region %d: %.2f mm - %.2f mm\n', i, z_t_3(qzs_start_3(i)), z_t_3(qzs_end_3(i)));
end




%% CH 5: Geometry schematic

%------- Spring Structure Maths Derivation, static 1 -----
a_1 = 50; %mm
b_1 = 80; %mm
% x0_1 = 20;

phi = pi/2:-pi/100:0;

y = 2*b_1*sin(phi);
x =  2*b_1*cos(phi);

figure()
plot(x,y,'LineWidth', 1.2);
grid on

xlabel('Horizontal displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('Vertical displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
% title('QZS Configuration: b=50mm', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
box on


%% Vertical force and displacement relation
% ------- Spring Structure Maths Derivation, static -----
% 第一段
a_1 = 50; %mm
b_1 = 80; %mm
l_0_1 = 20;
k_1 = (140/1E3*9.8)/(80); % N/mm
z_t_1 = [0:1:60]; % Displacement
l_1 = ( l_0_1^2 - 9/4*z_t_1.^2 + 3*z_t_1.*(3*b_1^2 - l_0_1^2)^0.5 ).^0.5; % Eq. 3-6
F_ver_1 = k_1*l_1.*(l_1 - l_0_1).*sqrt(3*b_1^2 - l_1.^2)/b_1^2;

% 1. Vertical displacement w.r.t. horizontal & vert. force
figure();
subplot(1,2,1)
[ax, h1, h2] = plotyy(z_t_1, l_1, z_t_1, F_ver_1);

% Label the axes
xlabel('Vertical displacement (mm)', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
ylabel(ax(1), 'Horizontal displacement (mm)', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
ylabel(ax(2), 'Vertical force (N)', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
set(h1, 'LineWidth', 2);
set(h2, 'LineWidth', 2);
grid on
% Add titles and legends
title('Vertical displacement relations', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
legend('Horizontal displacement', 'Vertical force (N)', 'FontName', 'Calibri', 'FontSize', 12, 'FontWeight', 'bold');

% 2. horizontal displacement w.r.t. Vertical & vert. force
% figure()
subplot(1,2,2)
[ax, h1, h2] = plotyy(l_1, z_t_1, l_1, F_ver_1);

% Labeling the axes
xlabel('Horizontal displacement (mm)', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
ylabel(ax(1), 'Vertical displacement (mm)', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
ylabel(ax(2), 'Vertical force (N)', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');

% Adding titles and legends
title('Horizontal displacement relations', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
legend('Vertical displacement', 'Vertical force', 'FontName', 'Calibri', 'FontSize', 12, 'FontWeight', 'bold');

% Adjusting the appearance of the plot
set(h1, 'LineWidth', 2);
set(h2, 'LineWidth', 2);
grid on;


%% Case 1: same k, various l
filename_b50_1 = 'Multiple_QZS/Geometry_50_60/b50_s04-5-20.xlsx';
data_b50_1 = readtable(filename_b50_1, 'Sheet', 'b50_s04-5-20');
load_b50_1 = data_b50_1.LoadValue;
defl_b50_1 = data_b50_1.PositionValue;

filename_b50_2 = 'Multiple_QZS/Geometry_50_60/b50_s05-6-35.xlsx';
data_b50_2 = readtable(filename_b50_2, 'Sheet', 'b50_s05-6-35');
load_b50_2 = data_b50_2.LoadValue;
defl_b50_2 = data_b50_2.PositionValue;


filename_b60_1 = 'Multiple_QZS/Geometry_50_60/b60_s04-5-20.xlsx';
data_b60_1 = readtable(filename_b60_1, 'Sheet', 'b60_s04-5-20');
load_b60_1 = data_b60_1.LoadValue;
defl_b60_1 = data_b60_1.PositionValue;

filename_b60_2 = 'Multiple_QZS/Geometry_50_60/b60_s05-6-35.xlsx';
data_b60_2 = readtable(filename_b60_2, 'Sheet', 'b60_s05-6-35');
load_b60_2 = data_b60_2.LoadValue;
defl_b60_2 = data_b60_2.PositionValue;

figure()
subplot(1,2,1)
hold on
plot(defl_b50_1,load_b50_1, 'LineWidth', 1.2);
plot(defl_b50_2,load_b50_2, 'LineWidth', 1.2);
xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
legend('0.4-5-20','0.5-6-35','FontName', 'Calibri', 'FontSize', 12, 'FontWeight', 'bold')
title('Same k, b=50mm', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
hold off
grid on
box on

subplot(1,2,2)
hold on
plot(defl_b60_1,load_b60_1, 'LineWidth', 1.2);
plot(defl_b60_2,load_b60_2, 'LineWidth', 1.2);
xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
legend('0.4-5-20','0.5-6-35','FontName', 'Calibri', 'FontSize', 12, 'FontWeight', 'bold')
title('Same k, b=60mm', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
hold off
grid on
box on

%% Case 1: same l, various k
clear all;clc;
filename_b50_1 =    'Multiple_QZS/Geometry_50_60/b50_s04-5-25.xlsx';
data_b50_1 = readtable(filename_b50_1, 'Sheet', 'b50_s04-5-25');
load_b50_1 = data_b50_1.LoadValue;
defl_b50_1 = data_b50_1.PositionValue;

filename_b50_2 =    'Multiple_QZS/Geometry_50_60/b50_s05-4-25.xlsx';
data_b50_2 = readtable(filename_b50_2, 'Sheet', 'b50_s05-4-25');
load_b50_2 = data_b50_2.LoadValue;
defl_b50_2 = data_b50_2.PositionValue;

filename_b50_3 =    'Multiple_QZS/Geometry_50_60/b50_s05-5-25.xlsx';
data_b50_3 = readtable(filename_b50_3, 'Sheet', 'b50_s05-5-25');
load_b50_3 = data_b50_3.LoadValue;
defl_b50_3 = data_b50_3.PositionValue;

filename_b60_1 =    'Multiple_QZS/Geometry_50_60/b60_s04-5-25.xlsx';
data_b60_1 = readtable(filename_b60_1, 'Sheet', 'b60_s04-5-25');
load_b60_1 = data_b60_1.LoadValue;
defl_b60_1 = data_b60_1.PositionValue;

filename_b60_2 =    'Multiple_QZS/Geometry_50_60/b60_s05-4-25.xlsx';
data_b60_2 = readtable(filename_b60_2, 'Sheet', 'b60_s05-4-25');
load_b60_2 = data_b60_2.LoadValue;
defl_b60_2 = data_b60_2.PositionValue;

filename_b60_3 =    'Multiple_QZS/Geometry_50_60/b60_s05-5-25.xlsx';
data_b60_3 = readtable(filename_b60_3, 'Sheet', 'b60_s05-5-25');
load_b60_3 = data_b60_3.LoadValue;
defl_b60_3 = data_b60_3.PositionValue;

figure()
subplot(1,2,1)
hold on
plot(defl_b50_1,load_b50_1, 'LineWidth', 1.2);
plot(defl_b50_2,load_b50_2, 'LineWidth', 1.2);
plot(defl_b50_3,load_b50_3, 'LineWidth', 1.2);
xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
legend('0.4-5-25','0.5-4-25','0.5-5-25','FontName', 'Calibri', 'FontSize', 12, 'FontWeight', 'bold')
title('Same initial length, b=50mm', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
hold off
grid on
box on

subplot(1,2,2)
hold on
plot(defl_b60_1,load_b60_1, 'LineWidth', 1.2);
plot(defl_b60_2,load_b60_2, 'LineWidth', 1.2);
plot(defl_b60_3,load_b60_3, 'LineWidth', 1.2);
xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
legend('0.4-5-25','0.5-4-25','0.5-5-25','FontName', 'Calibri', 'FontSize', 12, 'FontWeight', 'bold')
title('Same initial length, b=60mm', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
hold off
grid on
box on



