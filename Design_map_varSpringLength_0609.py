# 导入核心库
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import warnings
warnings.filterwarnings('ignore')  # 屏蔽数值计算的警告


# -------------------------- 全局设置（最终无无效参数，Calibri全局生效）--------------------------
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题
plt.rcParams['font.family'] = 'sans-serif'  # Calibri是无衬线字体，指定类别
plt.rcParams['font.sans-serif'] = ['Calibri', 'Microsoft YaHei', 'DejaVu Sans', 'Arial']  # 系统兼容备选
plt.rcParams['axes.titleweight'] = 'bold'   # 全局标题加粗
plt.rcParams['axes.labelweight'] = 'bold'   # 全局轴标签加粗
plt.rcParams['font.size'] = 12              # 全局基础字号
# 【修复点1】删除无效的contour.label.fontsize，等高线字号在clabel处单独设置

# -------------------------- 1. 按要求定义参数范围（核心：扩展参数维度）--------------------------
g = 9.8  # 重力加速度 m/s²
a_values = np.arange(40, 61, 5)  # a: [40,45,50,55,60] mm
x0_values = np.arange(4, 7, 1)   # x0: [4,5,6] mm （4:1:6）
b_values = np.arange(30, 101, 5) # b: [30,35,...,100] mm
# k: ([100:10:400]/1E3*g)/80 N/mm
k_numerator = np.arange(100, 401, 10)  # [100,110,...,400]
k_values = (k_numerator / 1e3 * g) / 80  # 计算k的最终值 N/mm

# 取a为中间值50（原MATLAB默认，也可循环a，此处贴合原逻辑）
a = 50  # mm

# -------------------------- 2. 初始化热力图数据（维度：b_len × k_len × x0_len）--------------------------
b_len = len(b_values)
k_len = len(k_values)
x0_len = len(x0_values)
# 存储满足条件的x1有效长度
heatmap_data_length = np.zeros((b_len, k_len, x0_len))
# 存储对应的最大F_Ver
heatmap_data_loading = np.zeros((b_len, k_len, x0_len))

# -------------------------- 3. 三重循环计算（x0 → b → k，核心：abs(k2)<=0.5k判定）--------------------------
for x0_idx in range(x0_len):
    x0 = x0_values[x0_idx]
    for b_idx in range(b_len):
        b = b_values[b_idx]
        # 限制asin参数在[-1,1]，避免NaN/报错
        phi0_arg = np.clip((b - x0) / b, -1, 1)
        phi0 = np.arcsin(phi0_arg)  # 初始角度
        
        # 生成phi1序列：从phi0递减，步长-0.05，终点合理范围（贴合原MATLAB）
        phi1_end_arg = np.clip((b - 50) / b, -1, 1)
        phi1_end = np.arcsin(phi1_end_arg)
        phi1 = np.arange(phi0, phi1_end - 0.01, -0.05)  # 左闭右开，减0.01确保包含终点
        if len(phi1) < 2:  # 点数不足，跳过计算
            continue
        
        # 计算位移和竖向力（原MATLAB公式1:1还原）
        disp_QZS = 2 * b - 2 * b * np.sin(phi1)
        F_Ver = (9/2) * k_values[:, np.newaxis] * b * np.sin(phi1) * (np.cos(phi1) - np.cos(phi0)) / np.cos(phi1)
        # 反转得到x1, y1（原MATLAB的end:-1:1）
        x1 = disp_QZS[::-1]
        y1 = F_Ver[:, ::-1]  # k为第一维，phi1为第二维
        
        for k_idx in range(k_len):
            k = k_values[k_idx]
            y1_k = y1[k_idx]
            # 计算刚度k2：y1对x1的数值微分（dy/dx，刚度定义）
            dx = np.diff(x1)
            dy = np.diff(y1_k)
            k2 = dy / dx  # k2为刚度，长度=len(x1)-1
            if len(k2) == 0:
                continue
            
            # 核心筛选：abs(k2) <= 0.5*k 筛选有效x1
            valid_k2_idx = np.abs(k2) <= 0.5 * k
            valid_x1 = x1[:-1][valid_k2_idx]  # 微分后x1长度减1，对应筛选
            
            # 计算有效x1的长度（删除了原x1长度输出，仅存储）
            if len(valid_x1) >= 2:
                total_x1_length = np.max(valid_x1) - np.min(valid_x1)
            else:
                total_x1_length = 0
            
            # 存储当前计算结果
            heatmap_data_length[b_idx, k_idx, x0_idx] = np.abs(total_x1_length)
            heatmap_data_loading[b_idx, k_idx, x0_idx] = np.max(y1_k) if np.max(y1_k) > 0 else 0
    
    # 打印进度（替代原MATLAB的disp）
    print(f'完成x0={x0_values[x0_idx]}mm / 共{x0_len}个x0值')

print('=== 计算完成，开始绘制设计图 ===')

# -------------------------- 4. 绘制核心二维热力图（横纵轴：b、k，贴合需求）--------------------------
fig1 = plt.figure(figsize=(16, 6), tight_layout=True)
# 子图1：QZS有效位移长度热力图（b-k轴）
ax1 = fig1.add_subplot(1, 2, 1)
im1 = ax1.imshow(heatmap_data_length[:, :, 0], extent=[k_values.min(), k_values.max(), b_values.min(), b_values.max()],
                 origin='lower', cmap='viridis', aspect='auto')
# 【修复点2】clabel处单独设置fontsize=10，替代无效的全局设置
cont1 = ax1.contour(k_values, b_values, heatmap_data_length[:, :, 0], colors='black', linewidths=1)
ax1.clabel(cont1, inline=True, fontsize=10)  # 直接指定等高线文本字号
cbar1 = fig1.colorbar(im1, ax=ax1, location='right')
cbar1.set_label('QZS Displacement Length [mm]')
ax1.set_title('QZS Range Design Map (x0=4mm)')
ax1.set_xlabel('k [N/mm]')
ax1.set_ylabel('b [mm]')

# 子图2：QZS最大荷载热力图（b-k轴）
ax2 = fig1.add_subplot(1, 2, 2)
im2 = ax2.imshow(heatmap_data_loading[:, :, 0], extent=[k_values.min(), k_values.max(), b_values.min(), b_values.max()],
                 origin='lower', cmap='turbo', aspect='auto')
cont2 = ax2.contour(k_values, b_values, heatmap_data_loading[:, :, 0], colors='black', linewidths=1)
ax2.clabel(cont2, inline=True, fontsize=10)  # 统一等高线文本字号
cbar2 = fig1.colorbar(im2, ax=ax2, location='right')
cbar2.set_label('QZS Max Loading [N]')
ax2.set_title('QZS Loading Design Map (x0=4mm)')
ax2.set_xlabel('k [N/mm]')
ax2.set_ylabel('b [mm]')

# -------------------------- 5. 绘制3D曲面图（多x0叠加，贴合原MATLAB的3D效果）--------------------------
fig2 = plt.figure(figsize=(10, 8), tight_layout=True)
ax3d = fig2.add_subplot(111, projection='3d')
# 生成b和k的网格（对应MATLAB的meshgrid）
K_grid, B_grid = np.meshgrid(k_values, b_values)
# 多x0叠加绘制，设置不同透明度
transparencies = [0.4, 0.6, 0.8]  # 对应x0[4,5,6]的透明度
for x0_idx in range(x0_len):
    Z = heatmap_data_loading[:, :, x0_idx]
    # 3D曲面：对应MATLAB的surf，edgecolor='none'消隐网格
    ax3d.plot_surface(K_grid, B_grid, Z, cmap='turbo', edgecolor='none',
                      alpha=transparencies[x0_idx], label=f'x0={x0_values[x0_idx]}mm')
    # 标注x0值（在曲面右上角）
    ax3d.text(k_values.max(), b_values.max(), Z.max(), f'x0={x0_values[x0_idx]}mm', fontsize=12)

# 3D图标签与样式
ax3d.set_title('QZS Loading Design Map (Multi-x0 Overlay)')
ax3d.set_xlabel('k [N/mm]')
ax3d.set_ylabel('b [mm]')
ax3d.set_zlabel('QZS Max Loading [N]')
ax3d.grid(True)
ax3d.legend(loc='upper left')
ax3d.view_init(elev=30, azim=-45)  # 调整3D视角

# -------------------------- 6. 绘制多x0子图热力图（2D，每个x0一个子图）--------------------------
fig3 = plt.figure(figsize=(15, 5), tight_layout=True)
for x0_idx in range(x0_len):
    ax = fig3.add_subplot(1, x0_len, x0_idx+1)
    im = ax.imshow(heatmap_data_loading[:, :, x0_idx], extent=[k_values.min(), k_values.max(), b_values.min(), b_values.max()],
                   origin='lower', cmap='turbo', aspect='auto')
    cont = ax.contour(k_values, b_values, heatmap_data_loading[:, :, x0_idx], colors='black', linewidths=1)
    ax.clabel(cont, inline=True, fontsize=10)  # 统一等高线字号
    cbar = fig3.colorbar(im, ax=ax, location='right')
    cbar.set_label('QZS Loading [N]')
    ax.set_title(f'x0={x0_values[x0_idx]}mm')
    ax.set_xlabel('k [N/mm]')
    if x0_idx == 0:  # 仅第一个子图显示y轴标签，避免重复
        ax.set_ylabel('b [mm]')

# 显示所有图
plt.show()

# 可选：保存高清无白边图片（论文/报告用，取消注释即可）
# fig1.savefig('QZS_2D_Heatmap.png', dpi=300, bbox_inches='tight', facecolor='white')
# fig2.savefig('QZS_3D_Surface.png', dpi=300, bbox_inches='tight', facecolor='white')
# fig3.savefig('QZS_Multi_x0_Heatmap.png', dpi=300, bbox_inches='tight', facecolor='white')