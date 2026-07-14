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
            
% 找出F_Ver位于0.95*max(F_Ver)到max(F_Ver)范围内的点
max_F = max(F_Ver);
valid_indices = F_Ver >= 0.95*max_F;
valid_x_1 = x_1(valid_indices);

% 计算对应的位移长度
if length(valid_x_1) >= 2
    total_x_1_length = max(valid_x_1) - min(valid_x_1);
else
    total_x_1_length = 0;
end

% 存储热力图数据
heatmap_data_length(b_idx, k_idx, x0_idx) = abs(total_x_1_length);
heatmap_data_loading(b_idx, k_idx, x0_idx) = max_F;

        end
    end
end

%% ver 2
clear all; clc;

a = 50; % mm
b_values = 50:0.2:100; % mm
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
            
            % Ensure phi_1 end point is valid and phi_0 is not less than it
            phi_1_end_val = asin(max(-1, min(1, (b-50)/b))); % Clamp argument to asin to [-1, 1]
            if phi_0 < phi_1_end_val
                phi_1 = []; % Or handle as an error/skip
            else
                phi_1 = [phi_0:-0.05:phi_1_end_val];
                % If phi_0 is very close to phi_1_end_val, phi_1 might contain only phi_0
                % or be [phi_0, phi_0 - small_step_that_is_still_>=_phi_1_end_val]
                % Ensure the last element is indeed phi_1_end_val if the step doesn't land on it exactly
                % and phi_0 is not equal to phi_1_end_val
                if ~isempty(phi_1) && phi_1(end) > phi_1_end_val && phi_1(end) - 0.05 < phi_1_end_val && phi_0 ~= phi_1_end_val
                    phi_1(end+1) = phi_1_end_val; % Add the exact end point if missed by a small margin
                elseif isempty(phi_1) && phi_0 == phi_1_end_val % Case where start and end are the same
                     phi_1 = phi_0;
                end

            end
            
            if isempty(phi_1) || length(phi_1) < 1 % Need at least one point for F_Ver
                total_x_1_length = 0;
                current_max_F_Ver = 0; % Or NaN, depending on how you want to represent this
            else
                disp_QZS = 2*b*ones(size(phi_1))-2*b*sin(phi_1);
                F_Ver = 9/2*k*b*sin(phi_1).*(cos(phi_1) - cos(phi_0)*ones(size(phi_1)) )./cos(phi_1);
                
                % Handle potential NaN or Inf in F_Ver if cos(phi_1) is near zero (phi_1 near pi/2)
                F_Ver(isinf(F_Ver) | isnan(F_Ver)) = []; % Remove them
                phi_1_valid_indices = ~(isinf(F_Ver) | isnan(F_Ver)); % if needed for disp_QZS
                disp_QZS = disp_QZS(phi_1_valid_indices);
                % Re-filter F_Ver just in case (should be redundant now)
                F_Ver = F_Ver(~(isinf(F_Ver) | isnan(F_Ver)));


                if isempty(F_Ver) || length(F_Ver) < 1
                    total_x_1_length = 0;
                    current_max_F_Ver = 0;
                else
%                     x_1 = disp_QZS(end:-1:1)-2*x0;
                    x_1 = disp_QZS(end:-1:1);
                    
                    y_1 = F_Ver(end:-1:1); % y_1 is the reversed F_Ver

                    current_max_F_Ver = max(y_1); % y_1 is the F_Ver used for this calculation

                    if current_max_F_Ver <= 0 || isempty(y_1) % If max force is not positive, or y_1 is empty
                        total_x_1_length = 0;
                    else
                        threshold = 0.95 * current_max_F_Ver;
                        
                        % Find indices where y_1 (F_Ver) is within the continuous range
                        % [max(F_Ver)*0.95, max(F_Ver)]
                        high_force_indices = find(y_1 >= threshold);
                        
                        if isempty(high_force_indices) || length(high_force_indices) < 2
                            % If no points or only one point is in the high force range,
                            % the length is considered 0.
                            total_x_1_length = 0;
                        else
                            % Get the x_1 values corresponding to these high forces
                            x_values_in_range = x_1(high_force_indices);
                            
                            % Calculate the length as the difference between the max and min x_1
                            % in this continuous high-force region
                            total_x_1_length = max(x_values_in_range) - min(x_values_in_range);
                        end
                    end
                end
            end
            
            % 存储热力图数据
            heatmap_data_length(b_idx, k_idx, x0_idx) = abs(total_x_1_length); % abs for safety
            if isempty(F_Ver) % If F_Ver was empty due to phi_1 issues or filtering
                 heatmap_data_loading(b_idx, k_idx, x0_idx) = 0; % Or NaN
            else
                 heatmap_data_loading(b_idx, k_idx, x0_idx) = max(F_Ver); % Original max F_Ver before reversal
            end
        end
    end
    disp(['Finished x0_idx = ', num2str(x0_idx), ' / ', num2str(num_x0)]); % Progress indicator
end

disp('Calculation complete.');

% Example of how to visualize one slice (optional)
% figure;
% x0_slice_idx = 1; % Choose an x0 slice to display
% imagesc(k_values, b_values, heatmap_data_length(:,:,x0_slice_idx));
% colorbar;
% xlabel('k (N/mm)');
% ylabel('b (mm)');
% title(['Heatmap of Displacement Length for x0 = ', num2str(x0_values(x0_slice_idx)), ' mm']);
% set(gca,'YDir','normal');


%% 绘制二维热力图
figure('Position', [100, 100, 1600, 600])

% 第一个子图：Total x_1 Length 的最底层
subplot(1, 2, 1)
imagesc(k_values, b_values, heatmap_data_length(:,:,1))
hold on
contour(k_values, b_values, heatmap_data_length(:,:,1), 'ShowText', 'on', 'LineColor', 'black')
colormap(parula)
colorbar('Location', 'eastoutside', 'FontSize', 12)
title('QZS Range Design Map','FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
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
title('QZS Loading Design Map','FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
xlabel('k [N/mm]','FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
ylabel('b [mm]', 'FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
set(gca, 'YDir', 'normal')
%%
% --- 绘图代码：为选定的 x0_values 绘制叠加的3D曲面 ---
% 假设 heatmap_data_length, k_values, b_values, x0_values 已经计算完毕

figure;
hold on; % 保持图形，以便在同一坐标轴上绘制多个曲面

% 创建 k_values 和 b_values 的网格作为 X 和 Y 坐标
[K_grid, B_grid] = meshgrid(k_values, b_values);

% --- 定义要绘制的 x0_values ---
% 您可以通过以下几种方式选择：

% 方式一: 选择特定索引 (例如：第一个, 中间的, 最后一个)
num_x0_total = length(x0_values);
% indices_to_plot = [1, round(num_x0_total/2), num_x0_total]; % 示例：选择3个
% indices_to_plot = [1, 50, 100, 150, 200, 250, num_x0_total]; % 示例：选择更多

% 方式二: 指定具体的 x0 值，然后找到它们对应的索引
target_x0_for_plot = [5, 12, 20]; % 示例：绘制 x0 = 5mm, 12mm, 20mm
% target_x0_for_plot = [x0_values(1), x0_values(round(num_x0_total*0.25)), x0_values(round(num_x0_total*0.5)), x0_values(round(num_x0_total*0.75)), x0_values(end)]; % 按比例选择

indices_to_plot = [];
for val_x0 = target_x0_for_plot
    [~, idx] = min(abs(x0_values - val_x0)); % 找到最接近的索引
    if ~isempty(idx)
        indices_to_plot = [indices_to_plot, idx(1)];
    end
end
indices_to_plot = unique(indices_to_plot); % 确保索引唯一且有序

if isempty(indices_to_plot)
    warning('没有选择任何 x0 值进行绘图，将默认选择第一个、中间和最后一个。');
    indices_to_plot = [1, round(num_x0_total/2), num_x0_total];
    indices_to_plot = unique(indices_to_plot); % 处理 num_x0_total < 2 的情况
end

% 为每个曲面定义颜色
% 您可以使用预定义的颜色列表或MATLAB的颜色图
num_surfaces_to_plot = length(indices_to_plot);
if num_surfaces_to_plot > 0
    % 使用 'turbo' colormap 中的一部分颜色，或者 'lines', 'jet', 'hsv' 等
    surface_colors = turbo(num_surfaces_to_plot + floor(num_surfaces_to_plot/2)); % 获取一些分散的颜色
    % 或者手动指定颜色:
    % custom_colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; ...];
    % if num_surfaces_to_plot <= size(custom_colors,1)
    %     surface_colors = custom_colors(1:num_surfaces_to_plot,:);
    % end
else
    disp('没有有效的 x0 索引可供绘图。');
    return; % 提前退出，因为没有东西可画
end


legend_entries = cell(1, num_surfaces_to_plot);
plot_handles = gobjects(1, num_surfaces_to_plot); % Preallocate graphics object array

for i = 1:num_surfaces_to_plot
    x0_idx = indices_to_plot(i);
    current_x0_value = x0_values(x0_idx);
    
    % 当前 x0 对应的 Z 数据 (曲面高度)
    Z_data_surface = heatmap_data_length(:,:,x0_idx);
    
    % 绘制曲面
    % 每个曲面使用一种独特的颜色，以便区分不同的 x0 值
    % 曲面的颜色是均匀的，其高度由 Z_data_surface 决定
    h_surf = surf(K_grid, B_grid, Z_data_surface, ...
                  'FaceColor', surface_colors(i,:), ...
                  'EdgeColor', 'none', ... % 可以设为 'black' 或其他颜色来看网格
                  'DisplayName', sprintf('x_0 = %.2f mm', current_x0_value)); % For legend
    
    % 可选：设置透明度，以便能看到被遮挡的曲面
    alpha(h_surf, 0.65); % 透明度范围 0 (完全透明) 到 1 (不透明)
    
    plot_handles(i) = h_surf; % 存储句柄用于图例
    % legend_entries{i} = sprintf('x_0 = %.2f mm', current_x0_value); % 也可以这样创建图例条目
end

hold off;

% --- 图形外观设置 ---
% colormap(turbo); % 此处的 colormap 主要影响 colorbar (如果添加的话)
                 % 因为我们为每个曲面直接指定了 FaceColor

% 添加图例来标识每个曲面对应的 x0 值
if num_surfaces_to_plot > 0 && all(isgraphics(plot_handles))
    legend(plot_handles, 'FontName', 'Calibri', 'Location', 'northwest', 'FontSize', 12);
    % 或者 legend(legend_entries, 'Location', ...); 如果 DisplayName 未在 surf 中设置
end

title_str = 'QZS Displacement Design Map';
% 如果只选择了一个x0值，可以调整标题
if num_surfaces_to_plot == 1
    title_str = sprintf('Displacement Length (x_0 = %.2f mm)', x0_values(indices_to_plot(1)));
end

title(title_str, 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
xlabel('k [N/mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
ylabel('b [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold');
zlabel('QZS Displacement Length [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold'); % Z轴是长度
% set(gca, 'YDir', 'normal', 'FontName', 'Calibri', 'FontSize', 12);
grid on;
view(3); % 设置为标准三维视角，例如 view([-37.5, 30])
% view([60 25]); % 可以调整视角参数 [azimuth, elevation]

axis tight; % 使坐标轴紧密贴合数据范围
box on; % 显示坐标区边框

% 关于 Colorbar:
% 由于每个曲面已经用不同颜色标示不同的 x0 值，
% 传统的 colorbar (将颜色映射到 Z 值) 在此处的直接意义不大。
% Z 轴本身已经表示了 Displacement Length。
% 如果确实需要 colorbar 来指示 Z 值的范围，可以这样做，但它不会改变曲面的颜色：
% cb = colorbar('Location', 'eastoutside', 'FontSize', 12);
% ylabel(cb, 'Displacement Length Scale [mm]', 'FontName', 'Calibri', 'FontSize', 12);
% colormap(gca, turbo); % 确保 colorbar 使用 turbo 颜色图





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

%% Default Font
figure()
trans = [0.3 0.5 0.7 0.9];
for x0_idx = 1:num_x0
    % 调整透明度
    surf(k_values, b_values, heatmap_data_loading(:,:,x0_idx), 'EdgeColor', 'none', 'FaceAlpha', trans(x0_idx));
    hold on
end
colormap(turbo)
colorbar('Location', 'eastoutside', 'FontSize', 12)
title('QZS Loading Design Map', 'FontSize', 12,'FontWeight', 'bold')
xlabel('k [N/mm]', 'FontSize', 24,'FontWeight', 'bold')
ylabel('b [mm]', 'FontSize', 24,'FontWeight', 'bold')
zlabel('QZS loading value [N]', 'FontSize', 24,'FontWeight', 'bold')
set(gca, 'YDir', 'normal')
% 在颜色条右侧添加 "QZS range" 图例
c = colorbar('Location', 'eastoutside', 'FontSize', 12);
c.Label.String = 'QZS loading';
c.Label.FontSize = 16;
c.Label.FontWeight = 'bold';
% 标注不同图层对应的 x0 值
for x0_idx = 1:num_x0
    text(max(k_values), max(b_values), max(max(heatmap_data_loading(:,:,x0_idx))), ...
         ['l_0 = ' num2str(x0_values(x0_idx)) ' mm'], 'FontSize', 12, 'FontWeight', 'bold')
end

box on;

% 设置坐标轴字体
set(gca, 'FontSize', 12)


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



