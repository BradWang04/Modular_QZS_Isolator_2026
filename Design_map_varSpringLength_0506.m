%% Dimensionless format of design map
%% Design map section
% 1. a取值范围为[40:5:60]，x0取值范围为[4:1:6]，b取值[30:5:100]，k取值为([100:10:400]/1E3*g)/(80)，
%    从而扩展其余相关参数的维度，并且删除输出x_1长度部分
% 2. 将abs(k_2) <= 0.25*k改为：abs(k_2) <= 0.5*k
% 3. 对每个a、x0、b和k的取值情况，找出abs(k_2) <= 0.5*k时所对应的x_1值，
%    并用二维热力图的方式进行表达，横纵坐标分别为b和k的值，并且添加对应的图表标注

clear all; clc;

a = 50; % mm
b_values = 50:1:100; % mm
g = 9.8;
x0_values = 5:0.05:20; % mm

k_values = ([100:5:2000]/1E3*g)/(80); % N/mm

% 初始化热力图数据
num_x0 = length(x0_values);
heatmap_data_length = zeros(length(b_values), length(k_values), num_x0);
heatmap_data_loading = zeros(length(b_values), length(k_values), num_x0);

for x0_idx = 1:num_x0
    x0 = x0_values(x0_idx);
    for b_idx = 1:length(b_values)
        b = b_values(b_idx);
        for k_idx = 1:length(k_values)
            k = k_values(k_idx);
            
            phi_0 = asin((b-x0)/b);
            
            phi_1 = [phi_0:-0.05:asin((b-50)/b)];
            
            disp_QZS = 2*b*ones(size(phi_1))-2*b*sin(phi_1);
            F_Ver = 9/2*k*b*sin(phi_1).*(cos(phi_1) - cos(phi_0)*ones(size(phi_1)) )./cos(phi_1);
            
            x_1 = disp_QZS(end:-1:1)-2*x0;
            y_1 = F_Ver(end:-1:1);
            
            % 使用数值微分方法对 y_1 关于 x_1 求导
            dx = diff(x_1);
            dy = diff(y_1);
            k_2 = dy ./ dx;
            
            % 找出满足 abs(k_2) <= 0.5*k 的 x_1 值
            valid_indices = abs(k_2) <= 0.5*k;
            valid_x_1 = x_1(1:end-1); % 注意 k_2 的长度比 x_1 少 1
            valid_x_1 = valid_x_1(valid_indices);
            
            % 计算对应的 x_1 长度
            total_x_1_length = sum(diff(valid_x_1));
            
            % 存储热力图数据
            heatmap_data_length(b_idx, k_idx, x0_idx) = abs(total_x_1_length);
            heatmap_data_loading(b_idx, k_idx, x0_idx) = max(F_Ver);
        end
    end
end


%% 绘制二维热力图
figure('Position', [100, 100, 1600, 600])

% 第一个子图：Total x_1 Length 的最底层
subplot(1, 2, 1)
imagesc(k_values, b_values, heatmap_data_length(:,:,1))
hold on
contour(k_values, b_values, heatmap_data_length(:,:,1), 'ShowText', 'on', 'LineColor', 'black')
colormap(parula)
colorbar('Location', 'eastoutside', 'FontSize', 12)
title('QZS Range Design Map (x0 = 4)','FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
xlabel('k [N/mm]','FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
ylabel('b [mm]','FontName', 'Calibri', 'FontSize', 15  ,'FontWeight', 'bold')
set(gca, 'YDir', 'normal')

% 第二个子图：F_Ver 的最大值 的最底层
subplot(1, 2, 2)
imagesc(k_values, b_values, heatmap_data_loading(:,:,1))
hold on
contour(k_values, b_values, heatmap_data_loading(:,:,1), 'ShowText', 'on', 'LineColor', 'black')
colormap(turbo)
colorbar('Location', 'eastoutside', 'FontSize', 12)
title('QZS Loading Design Map (x0 = 4)','FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
xlabel('k [N/mm]','FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
ylabel('b [mm]', 'FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
set(gca, 'YDir', 'normal')

%% 绘制三维热力图 
figure('Position', [100, 100, 1600, 600])
% % 第一个子图：Total x_1 Length 的三维图像
% subplot(1, 2, 1)
% % 使用 surf 绘制三维曲面
% surf(b_values, x0_values, squeeze(heatmap_data_length(:,1,:))', 'EdgeColor', 'none')
% hold on
% % 在三维曲面上绘制等高线
% contour3(b_values, x0_values, squeeze(heatmap_data_length(:,1,:))', 'ShowText', 'off', 'LineColor', 'black')
% colormap(parula)
% colorbar('Location', 'eastoutside', 'FontSize', 12)
% title('QZS Range Design Map (Case: k = 1.225 N/mm)','FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
% xlabel('b [mm]','FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
% ylabel('x0 [mm]','FontName', 'Calibri', 'FontSize', 15  ,'FontWeight', 'bold')
% zlabel('QZS vertical range [mm]','FontName', 'Calibri', 'FontSize', 15  ,'FontWeight', 'bold')
% set(gca, 'YDir', 'normal')
% box on
% 第一个子图：Total x_1 Length 的二维热力图
subplot(1, 2, 1)
% 使用 imagesc 绘制热力图
imagesc(b_values, x0_values, squeeze(heatmap_data_length(:,1,:))')
hold on
% 在热力图上绘制等高线
% contour(b_values, x0_values, squeeze(heatmap_data_length(:,1,:))', 'ShowText', 'on', 'LineColor', 'black')
colormap(parula)
colorbar('Location', 'eastoutside', 'FontSize', 12)
title('QZS Range Design Map','FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
xlabel('b [mm]','FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
ylabel('x0 [mm]','FontName', 'Calibri', 'FontSize', 15  ,'FontWeight', 'bold')
set(gca, 'YDir', 'normal')
box on

% 第二个子图：F_Ver 的最大值
subplot(1, 2, 2)
for x0_idx = 4:4
    surf(k_values, b_values, heatmap_data_loading(:,:,x0_idx), 'EdgeColor', 'none')
    hold on
end
contour3(k_values, b_values, heatmap_data_loading(:,:,x0_idx), 'ShowText', 'off', 'LineColor', 'black')
colormap(turbo)
colorbar('Location', 'eastoutside', 'FontSize', 12)
title('QZS Loading Design Map (Case: x_0=12mm)','FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
xlabel('k [N/mm]','FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
ylabel('b [mm]', 'FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
zlabel('QZS Loading [N]','FontName', 'Calibri', 'FontSize', 15  ,'FontWeight', 'bold')
set(gca, 'YDir', 'normal')
box on

%% Variable x0 scales ———— adding dimension for design map

clear all; clc;

a = 50; % mm
b_values = 50:1:100; % mm
g = 9.8;
x0_values = 5:5:20; % mm

k_values = ([100:5:2000]/1E3*g)/(80); % N/mm

% 初始化热力图数据
num_x0 = length(x0_values);
heatmap_data_length = zeros(length(b_values), length(k_values), num_x0);
heatmap_data_loading = zeros(length(b_values), length(k_values), num_x0);

for x0_idx = 1:num_x0
    x0 = x0_values(x0_idx);
    for b_idx = 1:length(b_values)
        b = b_values(b_idx);
        for k_idx = 1:length(k_values)
            k = k_values(k_idx);
            
            phi_0 = asin((b-x0)/b);
            
            phi_1 = [phi_0:-0.05:asin((b-50)/b)];
            
            disp_QZS = 2*b*ones(size(phi_1))-2*b*sin(phi_1);
            F_Ver = 9/2*k*b*sin(phi_1).*(cos(phi_1) - cos(phi_0)*ones(size(phi_1)) )./cos(phi_1);
            
            x_1 = disp_QZS(end:-1:1)-2*x0;
            y_1 = F_Ver(end:-1:1);
            
            % 使用数值微分方法对 y_1 关于 x_1 求导
            dx = diff(x_1);
            dy = diff(y_1);
            k_2 = dy ./ dx;
            
%             % 找出满足 abs(k_2) <= 0.5*k 的 x_1 值
%             valid_indices = abs(k_2) <= 0.5*k;
%             valid_x_1 = x_1(1:end-1); % 注意 k_2 的长度比 x_1 少 1
%             valid_x_1 = valid_x_1(valid_indices);
            
            % 找出满足 F_Ver >= 0.9 * MAX(F_Ver) 的 x_1 值
            max_F_Ver = max(F_Ver);
            valid_indices = F_Ver >= 0.9 * max_F_Ver;
            valid_x_1 = x_1(valid_indices);
            
            % 计算对应的 x_1 长度
            total_x_1_length = sum(diff(valid_x_1));
            
            % 存储热力图数据
            heatmap_data_length(b_idx, k_idx, x0_idx) = abs(total_x_1_length);
            heatmap_data_loading(b_idx, k_idx, x0_idx) = max(F_Ver);
        end
    end
end
%% Combined 3D map
figure()
trans = [0.3 0.5 0.7 0.9];
for x0_idx = 1:num_x0
    % 调整透明度
    surf(k_values, b_values, heatmap_data_loading(:,:,x0_idx), 'EdgeColor', 'none', 'FaceAlpha', trans(x0_idx));
    hold on
end
colormap(turbo)
colorbar('Location', 'eastoutside', 'FontSize', 24)
title('QZS Loading Design Map','FontName', 'Times New Roman', 'FontSize', 24,'FontWeight', 'bold')
xlabel('k [N/mm]','FontName', 'Times New Roman', 'FontSize', 24,'FontWeight', 'bold')
ylabel('b [mm]', 'FontName', 'Times New Roman', 'FontSize', 24,'FontWeight', 'bold')
zlabel('QZS loading value [N]','FontName', 'Times New Roman', 'FontSize', 24,'FontWeight', 'bold')
set(gca, 'YDir', 'normal')
% 在颜色条右侧添加 "QZS range" 图例
c = colorbar('Location', 'eastoutside', 'FontSize', 24);
c.Label.String = 'QZS loading';
c.Label.FontName = 'Times New Roman';
c.Label.FontSize = 24;
c.Label.FontWeight = 'bold';
% 标注不同图层对应的 x0 值
for x0_idx = 1:num_x0
    text(max(k_values), max(b_values), max(max(heatmap_data_loading(:,:,x0_idx))), ...
         ['x0 = ' num2str(x0_values(x0_idx)) ' mm'], 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
end

box on;

% 设置坐标轴字体
set(gca, 'FontName', 'Times New Roman', 'FontSize', 24)


%% Separate 2D maps
% close all;
% % 生成三张二维热力图，并添加等高线
% figure()
% % figure('Position', [100, 100, 1600, 1600])
% for x0_idx = 1:num_x0
%     subplot(2,2,x0_idx)
%     
%     % 绘制热力图
%     imagesc(k_values, b_values, heatmap_data_loading(:,:,x0_idx))
%     set(gca, 'YDir', 'normal')
%     colormap(turbo)
%     colorbar('Location', 'eastoutside', 'FontName', 'Times New Roman', 'FontSize', 24)
%     title(['x_0 = ' num2str(x0_values(x0_idx)) ' mm'], ...
%           'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold')
%     
%     % 添加等高线
%     hold on
%     contour(k_values, b_values, heatmap_data_loading(:,:,x0_idx), 'ShowText', 'on', 'LineColor', 'black')
%     hold off
%     
%     % 在颜色条右侧添加 "QZS range" 图例
%     c = colorbar('Location', 'eastoutside');
%     switch x0_idx
%         case 1
%             ylabel('b [mm]', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold')
%         case 2
%             c.Label.String = 'QZS loading';
%             c.Label.FontName = 'Times New Roman';
%             c.Label.FontSize = 24;
%             c.Label.FontWeight = 'bold';
%         case 3
%             xlabel('k [N/mm]', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold')
%             ylabel('b [mm]', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold')
%         case 4
%             c.Label.String = 'QZS loading';
%             c.Label.FontName = 'Times New Roman';
%             c.Label.FontSize = 24;
%             c.Label.FontWeight = 'bold';
%             xlabel('k [N/mm]', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold')
%     end
% 
%     % 调整整个图形的字体大小
%     set(gca, 'FontName', 'Times New Roman', 'FontSize', 24)
% end

%% 0506 - New version

% Assuming num_x0 is defined elsewhere
% Assuming k_values, b_values, heatmap_data_loading, x0_values are defined elsewhere.

% Define font names and sizes (adjust other sizes as needed based on your overall script)
newFontName = 'Times New Roman';
% Example sizes - adjust based on your overall plot design goals
titleFontSize = 34; % Example: if original title was 20, double it
axisTickFontSize = 30; % Example: if original axis ticks were 12, double it
contourTextFontSize = 20; % Specifically requested size for contour text

figure(); % Create a new figure

for x0_idx = 1:num_x0
    subplot(2,2,x0_idx); % Select the current subplot position

    % 绘制热力图
    imagesc(k_values, b_values, heatmap_data_loading(:,:,x0_idx));
    set(gca, 'YDir', 'normal'); % Set Y-axis direction
    colormap(turbo); % Apply colormap

    % Get current axes handle
    ax = gca;
    % Set axis tick font and size (Applies to axis ticks and potentially other inherited text)
    set(ax, 'FontName', newFontName, 'FontSize', axisTickFontSize);

    % Add the title for the subplot
    title(['x_0 = ' num2str(x0_values(x0_idx)) ' mm'], ...
          'FontName', newFontName, 'FontSize', titleFontSize, 'FontWeight', 'bold');

    % --- 修改后的等高线部分 ---
    hold on;
    % 绘制等高线，但不直接显示文本
    [C, h_contour] = contour(k_values, b_values, heatmap_data_loading(:,:,x0_idx), 'LineColor', 'black');
    % 使用 clabel 函数添加等高线文本，并设置字体和字号
    clabel(C, h_contour, 'FontName', newFontName, 'FontSize', contourTextFontSize);
    hold off;
    % --- 等高线修改结束 ---

    % Add colorbar
    % Setting colorbar font size - usually matches axis ticks or contour text
    c = colorbar('Location', 'eastoutside', 'FontSize', axisTickFontSize);

    % Customize colorbar label and potentially axis labels based on subplot index (as in original code)
    switch x0_idx
        case 1
            ylabel('b [mm]', 'FontName', newFontName, 'FontSize', axisTickFontSize, 'FontWeight', 'bold');
        case 2
            c.Label.String = 'QZS loading';
            c.Label.FontName = newFontName;
            c.Label.FontSize = axisTickFontSize; % Or contourTextFontSize? Choose based on desired appearance
            c.Label.FontWeight = 'bold';
        case 3
            xlabel('k [N/mm]', 'FontName', newFontName, 'FontSize', axisTickFontSize, 'FontWeight', 'bold');
            ylabel('b [mm]', 'FontName', newFontName, 'FontSize', axisTickFontSize, 'FontWeight', 'bold');
        case 4
            c.Label.String = 'QZS loading';
            c.Label.FontName = newFontName;
            c.Label.FontSize = axisTickFontSize; % Or contourTextFontSize? Choose based on desired appearance
            c.Label.FontWeight = 'bold';
            xlabel('k [N/mm]', 'FontName', newFontName, 'FontSize', axisTickFontSize, 'FontWeight', 'bold');
    end

    % Add grid and box for each subplot
    grid on;
    box on;

end % End of subplot loop

% Optional: Adjust the figure size after creating all subplots if needed
% set(gcf, 'Position', [100, 100, 1600, 1600]);




%% 创建自定义颜色映射
custom_colormap = turbo(256);

% 绘制 heatmap_data_length 与 b 的关系
figure()
% imagesc(x0_values, b_values, heatmap_data_length)

trans = [0.2 0.6 0.9];
for x0_idx = 1:num_x0
    % 调整透明度
    surf(k_values, b_values, heatmap_data_length(:,:,x0_idx), 'EdgeColor', 'none', 'FaceAlpha', trans(x0_idx));
    hold on
end

colormap(custom_colormap)
colorbar('Location', 'eastoutside', 'FontSize', 12)
title('QZS Length Design Map','FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
xlabel('Stiffness k','FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
ylabel('b [mm]', 'FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
set(gca, 'YDir', 'normal')
% 在颜色条右侧添加 "QZS Length" 图例
c = colorbar('Location', 'eastoutside', 'FontSize', 12);
c.Label.String = 'QZS Length';
c.Label.FontSize = 15;
c.Label.FontWeight = 'bold';



