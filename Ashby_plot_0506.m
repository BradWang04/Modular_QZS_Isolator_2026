close all; clear all;clc;
%%
% 示例数据
% materials = {'Metamaterial Isolator', 'Metamaterial B'};

% 散点图数据
frequency_1 =  [50, 35, 40, 60, 80];
loading_1 =    [50, 40, 70, 80, 65];
    
% frequency_2 =   [171.3195743	148.6132835	220.4860184	243.2914966	231.0380346	328.8978751	301.6289627	280.7080456	256.3460878	308.8162362];
% loading_2   =   [50.4759089	107.0795025	140.6124824	69.5647198	103.9751268	103.1358486	74.17121276	120.7797658	168.4898714	27.52451074];
frequency_2 =   [171.3195743	148.6132835	231.0380346	301.6289627	256.3460878	308.8162362];
loading_2   =   [50.4759089	    107.0795025	103.9751268	74.17121276	168.4898714	27.52451074];

% frequency_3 =   [128.0701754	105.2631579	56.14035088	63.15789474	73.68421053	207.0175439	170.1754386	280.7017544	138.5964912	100	170.1754386];
% loading_3   =   [13.63821612	60.0353206	88.66053408	51.06932153	144.0187859	184.2105263	113.3121798	135.7417327	181.8089194	100.9499689	55.79005589];
frequency_3 =   [128.0701754	207.0175439	170.1754386	138.5964912	100	170.1754386];
loading_3   =   [13.63821612	184.2105263	113.3121798	181.8089194	100.9499689	55.79005589];

% frequency_4 =   [71.92982456	140.3508772	77.19298246	66.66666667	133.3333333	205.2631579	185.9649123	231.5789474];
% loading_4   =   [207.0670315	180.7075765	131.8603478	159.4861046	106.5731253	71.37381618	118.8868188	127.8625213];
frequency_4 =   [71.92982456	77.19298246	133.3333333	185.9649123	231.5789474];
loading_4   =   [207.0670315	131.8603478	106.5731253	118.8868188	127.8625213];

frequency_5 =   [364.9122807	470.1754386	635.0877193	817.5438596	452.6315789	603.5087719	817.5438596	729.8245614	947.3684211];
loading_5   =   [257.6560317	178.3011178	134.509393	175.9431765	441.5269368	230.6609998	295.4122031	399.1519174	102.1871604];

% 定义颜色 -- 半透明
colors = ["#cecece", "#cfcfe5", "#b7d4ea", "#b8e3b2", "#fdc38d" , "#fcab8f", "#adadae", "#c2e6f7", "#e7c6db" , "#b3c7c8"]; 
color_full = ["#5f5f5f" , "#7262ac" , "#2e7ebb" , "2e974e" , "#e25508" , "#d92523"];

% 创建散点图
figure(1);
hold on
scatter(frequency_1, loading_1, 'filled','b');
scatter(frequency_2, loading_2, 'filled','r');
scatter(frequency_3, loading_3, 'filled');
scatter(frequency_4, loading_4, 'filled');
scatter(frequency_5, loading_5, 'filled');


% 添加标签
% for i = 1:length(materials)
%     text(frequency_beamBuckle(i), loading_beamBuckle(i), materials{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
% end
% text(frequency_beamBuckle(1), loading_beamBuckle(1), materials{1}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');


% 设置轴和标题
xlabel('Starting Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 30, 'FontWeight', 'bold');
ylabel('Loading (N)', 'FontName', 'Times New Roman', 'FontSize', 30, 'FontWeight', 'bold');
title('Vibration Isolator Operation Range - Loading and Starting Frequency', 'FontName', 'Times New Roman', 'FontSize', 25, 'FontWeight', 'bold');

% 设置轴范围
xlim([0 1000]);
ylim([0 500]);

% 绘制渐近线
res_N = 5:5:500;
res_m = res_N/9.8;
res_f = 1/(2*pi)*sqrt(5E7./res_m);

plot(res_f,res_N,'--','Color','#d92523','LineWidth',2)
plot(res_f*sqrt(2),res_N,'g','Color','#2e974e','LineWidth',2)


% 获取当前图形的轴对象
ax = gca;

%=============  椭圆绘制  ================
% 计算椭圆的中心位置
% 1. Metamaterial, 2. Buckling beam, 3. Magnets, 4. QZS Isolators
centerX = [60, 160, 160, 650, 250]; 
centerY = [60, 130, 120, 280, 100]; 

% 设置椭圆的半长轴和半短轴
semiMajorAxis = [40, 140, 140, 400 , 140];  
semiMinorAxis = [50,  60, 100 , 200 , 80];  

% 定义椭圆的角度和点数
theta = linspace(0, 2*pi, 100);
rotationAngle = [0, 150, 30, 340 , 170];  


% 计算椭圆的点坐标
numEllipses = length(centerX);
x_ellipse = zeros(numEllipses, length(theta));
y_ellipse = zeros(numEllipses, length(theta));

for i = 1:numEllipses
    x_ellipse(i, :) = centerX(i) + semiMajorAxis(i) * cos(theta) * cosd(rotationAngle(i)) - semiMinorAxis(i) * sin(theta) * sind(rotationAngle(i));
    y_ellipse(i, :) = centerY(i) + semiMajorAxis(i) * cos(theta) * sind(rotationAngle(i)) + semiMinorAxis(i) * sin(theta) * cosd(rotationAngle(i));
end

% 使用patch函数绘制椭圆
for i = 1:numEllipses
    patch('XData', x_ellipse(i, :), 'YData', y_ellipse(i, :), 'FaceColor', colors(i), 'FaceAlpha', 0.4, 'EdgeColor', 'none');
end

% 添加网格线
grid on;

% 添加黑色实线框
box on;

% 标注文字内容
text(200, 350, 'Resonance frequency: $f=\frac{1}{2\pi}\sqrt{\frac{k}{m}}$','Interpreter','latex', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','FontSize', 30, 'Color', '#d92523', 'FontName', 'Times New Roman');
text(300,300, 'Cut-off frequency: $f=\frac{\sqrt{2}}{2\pi}\sqrt{\frac{k}{m}}$','Interpreter','latex', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','FontSize', 30, 'Color', '#2e974e', 'FontName', 'Times New Roman');

text(30,0, 'Nonlinear beam design', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','FontSize', 24, 'FontName', 'Times New Roman');
text(500, 250, 'Linear Isolators', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','FontSize', 24, 'FontName', 'Times New Roman');
text(300, 45, 'Negative Stiffness Mechanism', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','FontSize', 24, 'FontName', 'Times New Roman');
text(150, 200, 'QZS Isolators', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','FontSize', 24, 'FontName', 'Times New Roman');
text(100, 150, 'Magnets', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','FontSize', 24, 'FontName', 'Times New Roman');

text(250, 100, 'This work', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 30, 'FontWeight', 'bold', 'Color', 'blue', 'FontName', 'Times New Roman');

% 调整坐标轴字体
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
% 调整图表大小
set(gcf, 'Position', [100, 100, 800, 600]);
%% Ultra-low frequency case

figure(2)
% 设置轴和标题
% 原始 FontSize 15 -> 双倍 30, FontName 'Calibri' -> 'Times New Roman'
xlabel('Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 30, 'FontWeight', 'bold');
ylabel('Loading (N)', 'FontName', 'Times New Roman', 'FontSize', 30, 'FontWeight', 'bold');
title('Ultra-low Frequency Vibration Isolator Operation Range', 'FontName', 'Times New Roman', 'FontSize', 30, 'FontWeight', 'bold');

% 获取当前坐标轴对象，用于设置刻度标签字体
ax = gca;
% 假设原始轴刻度标签 FontSize 是 12 (根据您之前提供的代码)，双倍后是 24
set(ax, 'FontName', 'Times New Roman', 'FontSize', 24);


grid on
box on
% 设置轴范围
xlim([0 50]);
ylim([0 100]);

hold on
%---- Ref 1
x_center = [1.8 10];  % 横轴线段的x坐标
y_center = [80 80];   % 横轴线段的y坐标

x_left = [1.8 1.8 2 2];  % 左侧工字端的x坐标
y_left = [78 82 82 78];      % 左侧工字端的y坐标

x_right = [9.8 9.8 10 10];  % 右侧工字端的x坐标
y_right = [78 82 82 78];     % 右侧工字端的y坐标
% 原始 FontSize 12 -> 双倍 24, FontName 'Calibri' -> 'Times New Roman'
text(5, 81, 'Chai, Y., et al. (2022)','FontName', 'Times New Roman', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 24);
% 绘制横轴线段
line(x_center, y_center, 'LineWidth', 2, 'Color', 'r');
% 绘制左侧工字端
patch(x_left, y_left, 'r', 'EdgeColor', 'r', 'LineWidth', 2);
% 绘制右侧工字端
patch(x_right, y_right, 'r', 'EdgeColor', 'r', 'LineWidth', 2);

%---- Ref 2
x_center = [20 40];  % 横轴线段的x坐标
y_center = [26.7 26.7];   % 横轴线段的y坐标

x_left = [20 20 20.2 20.2];  % 左侧工字端的x坐标
y_left = [24.7 28.7 28.7 24.7];      % 左侧工字端的y坐标

x_right = [39.8 39.8 40 40];  % 右侧工字端的x坐标
y_right = y_left;     % 右侧工字端的y坐标
% 原始 FontSize 12 -> 双倍 24, FontName 'Calibri' -> 'Times New Roman'
text(22, 27, 'Zhang, Q., et al. (2021)','FontName', 'Times New Roman', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 24);
% 绘制横轴线段
line(x_center, y_center, 'LineWidth', 2, 'Color', 'r');
% 绘制左侧工字端
patch(x_left, y_left, 'r', 'EdgeColor', 'r', 'LineWidth', 2);
% 绘制右侧工字端
patch(x_right, y_right, 'r', 'EdgeColor', 'r', 'LineWidth', 2);

%---- Ref 3
x_center = [12 50];  % 横轴线段的x坐标
y_center = [15 15];   % 横轴线段的y坐标

x_left = [12 12 12.2 12.2];  % 左侧工字端的x坐标
y_left = [13 17 17 13];      % 左侧工字端的y坐标

x_right = [49.8 49.8 50 50];  % 右侧工字端的x坐标
y_right = y_left;     % 右侧工字端的y坐标
% 原始 FontSize 12 -> 双倍 24, FontName 'Calibri' -> 'Times New Roman'
text(25, 15, 'Jiahao Zhou, et al.(2024)','FontName', 'Times New Roman', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 24);
% 绘制横轴线段
line(x_center, y_center, 'LineWidth', 2, 'Color', 'r');
% 绘制左侧工字端
patch(x_left, y_left, 'r', 'EdgeColor', 'r', 'LineWidth', 2);
% 绘制右侧工字端
patch(x_right, y_right, 'r', 'EdgeColor', 'r', 'LineWidth', 2);

%---- Ref 4
x_center = [15 40];  % 横轴线段的x坐标
y_center = [5 5];   % 横轴线段的y坐标

x_left = [15 15 15.2 15.2];  % 左侧工字端的x坐标
y_left = [3 7 7 3];      % 左侧工字端的y坐标

x_right = [39.8 39.8 40 40];  % 右侧工字端的x坐标
y_right = y_left;     % 右侧工字端的y坐标
% 原始 FontSize 12 -> 双倍 24, FontName 'Calibri' -> 'Times New Roman'
text(20, 5, 'Banerjee, P., et al. (2023)','FontName', 'Times New Roman', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 24);
% 绘制横轴线段
line(x_center, y_center, 'LineWidth', 2, 'Color', 'r');
% 绘制左侧工字端
patch(x_left, y_left, 'r', 'EdgeColor', 'r', 'LineWidth', 2);
% 绘制右侧工字端
patch(x_right, y_right, 'r', 'EdgeColor', 'r', 'LineWidth', 2);

%---- Ref 5
x_center = [5 12];  % 横轴线段的x坐标
y_center = [50 50];   % 横轴线段的y坐标

x_left = [5 5 5.2 5.2];  % 左侧工字端的x坐标
y_left = [48 52 52 48];      % 左侧工字端的y坐标

x_right = [11.8 11.8 12 12];  % 右侧工字端的x坐标
y_right = y_left;     % 右侧工字端的y坐标
% 原始 FontSize 12 -> 双倍 24, FontName 'Calibri' -> 'Times New Roman'
text(6, 51, 'Ye, K., et al. (2020)','FontName', 'Times New Roman', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 24);
% 绘制横轴线段
line(x_center, y_center, 'LineWidth', 2, 'Color', 'r');
% 绘制左侧工字端
patch(x_left, y_left, 'r', 'EdgeColor', 'r', 'LineWidth', 2);
% 绘制右侧工字端
patch(x_right, y_right, 'r', 'EdgeColor', 'r', 'LineWidth', 2);

%---- Ref 6
x_center = [2 5];  % 横轴线段的x坐标
y_center = [11.2 11.2];   % 横轴线段的y坐标

x_left = [2 2 2.2 2.2];  % 左侧工字端的x坐标
y_left = [9.2 13.2 13.2 9.2];      % 左侧工字端的y坐标

x_right = [4.8 4.8 5 5];  % 右侧工字端的x坐标
y_right = y_left;     % 右侧工字端的y坐标
% 原始 FontSize 12 -> 双倍 24, FontName 'Calibri' -> 'Times New Roman'
text(2, 16, 'Yang, H. and H. Zhao (2024)','FontName', 'Times New Roman', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 24);
% 绘制横轴线段
line(x_center, y_center, 'LineWidth', 2, 'Color', 'r');
% 绘制左侧工字端
patch(x_left, y_left, 'r', 'EdgeColor', 'r', 'LineWidth', 2);
% 绘制右侧工字端
patch(x_right, y_right, 'r', 'EdgeColor', 'r', 'LineWidth', 2);

%---- Ref 7
x_center = [25 45];  % 横轴线段的x坐标
y_center = [65 65];   % 横轴线段的y坐标

x_left = [25 25 25.2 25.2];  % 左侧工字端的x坐标
y_left = [63 67 67 63];      % 左侧工字端的y坐标

x_right = [44.8 44.8 45 45];  % 右侧工字端的x坐标
y_right = y_left;     % 右侧工字端的y坐标
% 原始 FontSize 12 -> 双倍 24, FontName 'Calibri' -> 'Times New Roman'
text(27, 65, 'Ling, P., et al. (2022)','FontName', 'Times New Roman', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 24);
% 绘制横轴线段
line(x_center, y_center, 'LineWidth', 2, 'Color', 'r');
% 绘制左侧工字端
patch(x_left, y_left, 'r', 'EdgeColor', 'r', 'LineWidth', 2);
% 绘制右侧工字端
patch(x_right, y_right, 'r', 'EdgeColor', 'r', 'LineWidth', 2);

hold off; % 确保在段落结束时释放 hold
