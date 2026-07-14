# -*- coding: utf-8 -*-
"""
QZS Vibration Isolator - Force-Displacement Curve: High-Order Taylor Expansion
最终版修改：1.子图右下角+高放大倍率 2.全局Calibri字体 3.所有字体放大1倍 4.修复所有报错
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# ====================== 全局设置（Calibri字体+绘图样式+负号正常显示） ======================
plt.rcParams['font.sans-serif'] = ['Calibri', 'DejaVu Sans']  # 全局Calibri字体
plt.rcParams['axes.unicode_minus'] = False  # 负号正常显示
# 论文级绘图基础样式
mpl.rcParams['lines.linewidth'] = 1.2
mpl.rcParams['axes.linewidth'] = 1.2
mpl.rcParams['grid.linestyle'] = '--'
mpl.rcParams['grid.alpha'] = 0.7
mpl.rcParams['font.family'] = 'sans-serif'  # 配合Calibri生效

# ====================== 1. 系统参数与原始曲线计算 ======================
b = 80e-3;        # 连杆长度 (m)
k = 300;          # 弹簧刚度 (N/m)
x0 = 5e-3;        # 初始偏移量 (m)
g = 9.81;         # 重力加速度 (m/s²)

# 几何角度限制
phi_lim = np.arcsin((b - 50e-3) / b);
n_phi = 1000;     # 高采样点保证微分/拟合精度
phi0 = np.arcsin((b - x0) / b);
phi_array = np.linspace(phi0, phi_lim, n_phi);  # 摆角数组

# 1.1 原始运动学与力学计算
z_disp = 2 * b * (np.sin(phi0) - np.sin(phi_array));  # 垂直位移 (m)
F_Ver = (9 * k * b / 2) * np.sin(phi_array) * (np.cos(phi_array) - np.cos(phi0)) / np.cos(phi_array);

# 1.2 寻找静平衡位置 (力最大处为准零刚度点)
F_eq = np.max(F_Ver);
idx_eq = np.argmax(F_Ver);
z_eq = z_disp[idx_eq];
u = z_disp - z_eq;        
F_restore_original = F_Ver - F_eq;  

# ====================== 2. 数值泰勒展开 (多项式拟合提取刚度系数) ======================
fit_window = 2e-3; 
mask = np.abs(u) < fit_window;
u_fit = u[mask];
F_fit = F_restore_original[mask];

# 3阶/5阶多项式拟合（降序系数解包）
p3 = np.polyfit(u_fit, F_fit, 3);
k3_num, k2_num, k1_num, _ = p3;
p5 = np.polyfit(u_fit, F_fit, 5);
k5_high, k4_high, k3_high, k2_high, k1_high, _ = p5;

# 格式化输出提取的刚度参数
print("====================== 数值泰勒展开系数 ======================")
print(f"拟合范围: ±{fit_window*1000:.1f} mm")
print("----------------------------------------------------------")
print(f"线性刚度 k1 (N/m)   : {k1_num:.4f} \t(理论应接近0)")
print(f"二次刚度 k2 (N/m^2) : {k2_num:.4f} \t(非对称项)")
print(f"三次刚度 k3 (N/m^3) : {k3_num:.4e} \t(主导非线性项)")
print(f"五次刚度 k5 (N/m^5) : {k5_high:.4e} \t(大变形修正项)")
print("============================================================\n")

# ====================== 3. 重构近似曲线 + 误差计算 ======================
# 3阶/5阶泰勒近似恢复力
F_taylor_3rd = k1_num * u + k2_num * (u ** 2) + k3_num * (u ** 3);
F_taylor_5th = (k1_high * u + k2_high * (u ** 2) + k3_high * (u ** 3) + 
                k4_high * (u ** 4) + k5_high * (u ** 5));
# 重构总力
F_total_3rd = F_eq + F_taylor_3rd;
F_total_5th = F_eq + F_taylor_5th;

# 计算RMSE拟合误差
rmse_3rd = np.sqrt(np.mean((F_restore_original - F_taylor_3rd) ** 2));
rmse_5th = np.sqrt(np.mean((F_restore_original - F_taylor_5th) ** 2));
print("拟合误差 (RMSE):")
print(f"3阶展开误差: {rmse_3rd:.4f} N")
print(f"5阶展开误差: {rmse_5th:.4f} N (更优)")

# ====================== 5. 绘图对比（核心：所有字体放大1倍 + 之前所有修改） ======================
fig = plt.figure(figsize=(10, 6), facecolor='w')  # 主图尺寸1000x600
ax_main = fig.add_subplot(111)

# 绘制主图曲线
ax_main.plot(z_disp*1000, F_Ver, 'b-', linewidth=2.5, label='Original Analytical Model')
ax_main.plot(z_disp*1000, F_total_3rd, 'r--', linewidth=2, label='Taylor Expansion (3rd Order)')
ax_main.plot(z_disp*1000, F_total_5th, 'g-.', linewidth=1.5, label='Taylor Expansion (5th Order)')
ax_main.plot(z_eq*1000, F_eq, 'ko', markerfacecolor='k', markersize=8, label='Equilibrium Point') # 标记点同步放大

# 主图标注：所有fontsize放大1倍（原12→24，原14→28，原11→22）
ax_main.set_xlabel('Vertical Displacement z (mm)', fontsize=24, fontweight='bold')  # 原12→24
ax_main.set_ylabel('Vertical Force F (N)', fontsize=24, fontweight='bold')          # 原12→24
ax_main.set_title('QZS Curve Approximation: Original vs Corrected Taylor Series', fontsize=20, fontweight='bold')  # 原14→28
ax_main.legend(loc='lower right', fontsize=16)  # 原11→22
ax_main.grid(True, linewidth=1.5)  # 网格线同步加粗，匹配大字体
# 纵轴范围适配
y_min, y_max = np.min(F_Ver)*1.1, np.max(F_Ver)*1.1
ax_main.set_ylim([y_min, y_max])
# 主图坐标轴刻度字体放大（新增，匹配整体大字体）
ax_main.tick_params(axis='both', labelsize=20)

# 子图：右下角 + 高放大倍率（位置/倍率不变，字体同步放大1倍）
ax_zoom = fig.add_axes([0.3, 0.25, 0.3, 0.3])  # 右下角归一化坐标
# 绘制放大图曲线
ax_zoom.plot(z_disp*1000, F_Ver, 'b-', linewidth=2)
ax_zoom.plot(z_disp*1000, F_total_3rd, 'r--', linewidth=1.5)
ax_zoom.plot(z_disp*1000, F_total_5th, 'g-.', linewidth=1.5)
# 高放大倍率范围（±5mm，F_eq±3）
x_zoom_min = (z_eq - 0.005) * 1000
x_zoom_max = (z_eq + 0.005) * 1000
ax_zoom.set_xlim([x_zoom_min, x_zoom_max])
ax_zoom.set_ylim([F_eq - 3, F_eq + 3])
# 放大图标注：字体放大1倍（原11→22）+ 刻度字体放大
ax_zoom.set_title('Zoom at Equilibrium', fontsize=22, fontweight='bold')  # 原11→22
ax_zoom.grid(True, linewidth=1.5)
ax_zoom.tick_params(axis='both', labelsize=20)  # 新增刻度大字体，匹配整体

# 防止大字体重叠，优化布局
plt.tight_layout(pad=2)
plt.show()