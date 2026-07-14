%% Design map section
% 1. a取值不变，b取值[30:5:100]，k取值为([100:10:400]/1E3*g)/(80)，
%    从而扩展其余相关参数的维度，并且删除输出x_1长度部分
% 2. 将abs(k_2) <= 0.25*k改为：abs(k_2) <= 0.5*k
% 3. 对每个b和k的取值情况，找出abs(k_2) <= 0.5*k时所对应的x_1值，
%    并用二维热力图的方式进行表达，横纵坐标分别为b和k的值，并且添加对应的图表标注

clear all; clc;

a = 50; % mm
b_values = 50:1:120; % mm
g = 9.8;
x0 = 5;

k_values = ([100:100:20000]/1E3*g)/(80); % N/mm

% 初始化热力图数据
heatmap_data_length = zeros(length(b_values), length(k_values));
heatmap_data_loading = zeros(length(b_values), length(k_values));

for b_idx = 1:length(b_values)
    b = b_values(b_idx);
    for k_idx = 1:length(k_values)
        k = k_values(k_idx);
        
        phi_0 = asin((b-x0)/b);
        
        phi_1 = [phi_0:-0.05:asin((b-50)/b)];
        
        disp = 2*b*ones(size(phi_1))-2*b*sin(phi_1);
        F_Ver = 9/2*k*b*sin(phi_1).*(cos(phi_1) - cos(phi_0)*ones(size(phi_1)) )./cos(phi_1);
        
        x_1 = disp(end:-1:1)-2*x0;
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
        heatmap_data_length(b_idx, k_idx) = abs(total_x_1_length);
        heatmap_data_loading(b_idx, k_idx) = max(F_Ver);
    end
end

%% 绘制热力图
figure('Position', [100, 100, 1600, 600])

% 第一个子图：Total x_1 Length
subplot(1, 2, 1)
imagesc(k_values, b_values, heatmap_data_length)
hold on
contour(k_values, b_values, heatmap_data_length, 'ShowText', 'on', 'LineColor', 'black')
colormap(parula)

cb1 = colorbar('Location', 'eastoutside', 'FontSize', 12);
cb1.Label.String = 'QZS range [mm]';
cb1.Label.FontName = 'Calibri';
cb1.Label.FontSize = 15;
cb1.Label.FontWeight = 'bold';
cb1.Label.Color = 'blue';
% title('QZS Range Design Map','FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
xlabel('k [N/mm]','FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
ylabel('b [mm]','FontName', 'Calibri', 'FontSize', 15  ,'FontWeight', 'bold')
set(gca, 'YDir', 'normal')


% 第二个子图：F_Ver 的最大值
subplot(1, 2, 2)
imagesc(k_values, b_values, heatmap_data_loading)
hold on
contour(k_values, b_values, heatmap_data_loading, 'ShowText', 'on', 'LineColor', 'black')
cb2 = colorbar('Location', 'eastoutside', 'FontSize', 12);
cb2.Label.String = 'QZS loading [N]';
cb2.Label.FontName = 'Calibri';
cb2.Label.FontSize = 15;
cb2.Label.FontWeight = 'bold';
cb2.Label.Color = 'red';
% title('QZS Loading Design Map','FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
xlabel('k [N/mm]','FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
ylabel('b [mm]', 'FontName', 'Calibri', 'FontSize', 15,'FontWeight', 'bold')
set(gca, 'YDir', 'normal')
colormap(turbo)


%===== 嵌入放大图像 =====
%===== 嵌入放大图像 =====
% 标记放大的区域
rectangle('Position', [0.05, 50, 0.5, 30], 'EdgeColor', 'y', 'LineStyle', '-','LineWidth',3)

% 创建可拖拽的放大区域
pan on;
axes('Position', [0.593, 0.151, 0.2, 0.2]); % 设置嵌入区域的位置和大小
imagesc(k_values(k_values >= 0 & k_values <= 0.5), b_values(b_values >= 50 & b_values <= 80), heatmap_data_loading(b_values >= 50 & b_values <= 80, k_values >= 0 & k_values <= 0.5))
hold on
contour(k_values(k_values >= 0 & k_values <= 0.5), b_values(b_values >= 50 & b_values <= 80), heatmap_data_loading(b_values >= 50 & b_values <= 80, k_values >= 0 & k_values <= 0.5), 'ShowText', 'on', 'LineColor', 'black')
set(gca, 'YDir', 'normal')
colormap(turbo)
axis tight
axis square

% 修改x和y轴的字体大小和颜色
set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'XColor', 'yellow', 'YColor', 'yellow');
% xlabel('k [N/mm]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'blue');
% ylabel('b [mm]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'blue');

%% 


