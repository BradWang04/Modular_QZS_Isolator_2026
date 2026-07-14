# -*- coding: utf-8 -*-
"""
合并版：QZS力-位移泰勒展开 + 自由振动衰减曲线（上下排布 + (a)(b)标注）
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# ====================== 全局字体/样式设置（参考第二组代码，统一Calibri风格） ======================
plt.rcParams['font.sans-serif'] = ['Calibri', 'DejaVu Sans']  # 全局Calibri字体
plt.rcParams['axes.unicode_minus'] = False  # 负号正常显示
plt.rcParams['axes.labelweight'] = 'bold'   # 轴标题加粗（核心参考第二组）
plt.rcParams['font.size'] = 16              # 基础字号（参考第二组）
# 通用绘图样式
mpl.rcParams['lines.linewidth'] = 1.2
mpl.rcParams['axes.linewidth'] = 1.2
mpl.rcParams['grid.linestyle'] = '--'
mpl.rcParams['grid.alpha'] = 0.5            # 参考第二组的网格透明度
mpl.rcParams['figure.facecolor'] = 'white'
mpl.rcParams['savefig.dpi'] = 300           # 高清保存

# ====================== 第一部分：QZS力-位移曲线（泰勒展开）计算 ======================
# 系统参数
b = 80e-3;        # 连杆长度 (m)
k_qzs = 300;      # 弹簧刚度 (N/m)
x0 = 5e-3;        # 初始偏移量 (m)
g = 9.81;         # 重力加速度 (m/s²)

# 几何角度限制与摆角数组
phi_lim = np.arcsin((b - 50e-3) / b);
n_phi = 1000;     # 高采样点保证精度
phi0 = np.arcsin((b - x0) / b);
phi_array = np.linspace(phi0, phi_lim, n_phi);

# 原始运动学与力学计算
z_disp = 2 * b * (np.sin(phi0) - np.sin(phi_array));  # 垂直位移 (m)
F_Ver = (9 * k_qzs * b / 2) * np.sin(phi_array) * (np.cos(phi_array) - np.cos(phi0)) / np.cos(phi_array);

# 静平衡位置（准零刚度点）
F_eq = np.max(F_Ver);
idx_eq = np.argmax(F_Ver);
z_eq = z_disp[idx_eq];
u = z_disp - z_eq;        
F_restore_original = F_Ver - F_eq;  

# 数值泰勒展开（3阶/5阶多项式拟合）
fit_window = 2e-3; 
mask = np.abs(u) < fit_window;
u_fit = u[mask];
F_fit = F_restore_original[mask];

p3 = np.polyfit(u_fit, F_fit, 3);
k3_num, k2_num, k1_num, _ = p3;
p5 = np.polyfit(u_fit, F_fit, 5);
k5_high, k4_high, k3_high, k2_high, k1_high, _ = p5;

# 输出刚度参数
print("====================== 数值泰勒展开系数 ======================")
print(f"拟合范围: ±{fit_window*1000:.1f} mm")
print(f"线性刚度 k1 (N/m)   : {k1_num:.4f} \t(理论应接近0)")
print(f"二次刚度 k2 (N/m^2) : {k2_num:.4f} \t(非对称项)")
print(f"三次刚度 k3 (N/m^3) : {k3_num:.4e} \t(主导非线性项)")
print(f"五次刚度 k5 (N/m^5) : {k5_high:.4e} \t(大变形修正项)")

# 重构近似曲线 + 误差计算
F_taylor_3rd = k1_num * u + k2_num * (u ** 2) + k3_num * (u ** 3);
F_taylor_5th = (k1_high * u + k2_high * (u ** 2) + k3_high * (u ** 3) + 
                k4_high * (u ** 4) + k5_high * (u ** 5));
F_total_3rd = F_eq + F_taylor_3rd;
F_total_5th = F_eq + F_taylor_5th;

rmse_3rd = np.sqrt(np.mean((F_restore_original - F_taylor_3rd) ** 2));
rmse_5th = np.sqrt(np.mean((F_restore_original - F_taylor_5th) ** 2));
print("----------------------------------------------------------")
print(f"3阶展开RMSE误差: {rmse_3rd:.4f} N")
print(f"5阶展开RMSE误差: {rmse_5th:.4f} N (更优)")
print("============================================================\n")

# ====================== 第二部分：自由振动衰减曲线计算 ======================
# 系统参数（与论文一致）
m = 50.0                # 质量 kg
k_vib = 78957.0         # 刚度 N/m
zeta = 0.0801           # 阻尼比
omega_n = np.sqrt(k_vib / m)
omega_d = omega_n * np.sqrt(1 - zeta**2)
T_d = 2 * np.pi / omega_d

A = 5.0                 # 初始位移 mm
fs = 1000               # 采样率 Hz
t_total = 1.6           # 总时间 s


df = pd.read_csv('free_vibration_decay_data.csv')

t       = df['Time_s'].to_numpy()
x_clean = df['Disp_clean_mm'].to_numpy()
x_noisy = df['Disp_noisy_mm'].to_numpy()
cycle   = df['Cycle'].to_numpy()

fs      = 1.0 / np.mean(np.diff(t))   # 采样率 (Hz)
t_total = t[-1] + 1.0/fs              # 总时长 (s)
A       = x_clean[0]                  # 初始位移 (mm)

# ====================== 合并绘图（2行1列上下布局） ======================
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), facecolor='w')  # 2行1列，竖版画布

# -------- 子图1（上）：力-位移泰勒展开曲线 --------
# 主曲线绘制
ax1.plot(z_disp*1000, F_Ver, 'b-', linewidth=2.5, label='Original Analytical Model')
ax1.plot(z_disp*1000, F_total_3rd, 'r--', linewidth=2, label='Taylor Expansion (3rd Order)')
ax1.plot(z_disp*1000, F_total_5th, 'g-.', linewidth=1.5, label='Taylor Expansion (5th Order)')
ax1.plot(z_eq*1000, F_eq, 'ko', markerfacecolor='k', markersize=8, label='Equilibrium Point')

# 样式设置（参考第二组：轴标题加粗、Calibri字体）
ax1.set_xlabel('Vertical Displacement z (mm)', fontweight='bold', fontsize=16)
ax1.set_ylabel('Vertical Force F (N)', fontweight='bold', fontsize=16)
ax1.set_title('QZS Curve Approximation: Original vs Taylor Series', fontweight='bold', fontsize=18)
ax1.legend(loc='lower right', fontsize=14)
ax1.grid(True, alpha=0.5)
ax1.tick_params(axis='both', labelsize=14)
# 纵轴范围适配
y_min, y_max = np.min(F_Ver)*1.1, np.max(F_Ver)*1.1
ax1.set_ylim([y_min, y_max])

# 子图1的局部放大框（按要求修改位置：[0.2, 0.42, 0.2, 0.2]）
ax_zoom = fig.add_axes([0.4, 0.6, 0.2, 0.2])  # 归一化坐标，适配上下布局
ax_zoom.plot(z_disp*1000, F_Ver, 'b-', linewidth=2)
ax_zoom.plot(z_disp*1000, F_total_3rd, 'r--', linewidth=1.5)
ax_zoom.plot(z_disp*1000, F_total_5th, 'g-.', linewidth=1.5)
# 放大范围
x_zoom_min = (z_eq - 0.002) * 1000
x_zoom_max = (z_eq + 0.002) * 1000
ax_zoom.set_xlim([x_zoom_min, x_zoom_max])
ax_zoom.set_ylim([F_eq - 1, F_eq + 1])
ax_zoom.set_title('Zoom at Equilibrium', fontweight='bold', fontsize=14)
ax_zoom.grid(True)
ax_zoom.tick_params(axis='both', labelsize=12)

# -------- 子图2（下）：自由振动衰减曲线 --------
ax2.plot(t, x_clean, 'b-', linewidth=1.5, label='Ideal decay model')
ax2.plot(t, x_noisy, 'r.', markersize=2.5, alpha=0.7, label='Measured signal')

# 样式设置（参考第二组，修正重复的xlabel/ylabel）
ax2.set_xlabel('Time (s)', fontweight='bold', fontsize=16)
ax2.set_ylabel('Displacement (mm)', fontweight='bold', fontsize=16)
ax2.set_title('Free Vibration Displacement Decay\n'
              r'$\zeta=%.3f,\ f_n\approx6.30\ \mathrm{Hz},\ T_d\approx%.4f\ \mathrm{s}$'
              % (zeta, T_d), fontweight='bold', fontsize=18)
ax2.grid(True, alpha=0.3)
ax2.legend(fontsize=14)
ax2.tick_params(axis='both', labelsize=14)

# -------- 添加(a)(b)标注（外部左上角，Calibri+加粗） --------
# (a)标注：对应上侧力-位移图
fig.text(0.02, 0.92, 'a', fontsize=26, fontweight='bold', fontfamily='Calibri')
# (b)标注：对应下侧振动衰减图
fig.text(0.02, 0.42, 'b', fontsize=26, fontweight='bold', fontfamily='Calibri')

# 整体布局优化 + 保存图片
plt.tight_layout(pad=2)
plt.subplots_adjust(left=0.12, top=0.9, bottom=0.08)  # 调整左侧边距，避免标注被遮挡
plt.savefig('combined_qzs_vibration_plot.png', dpi=300, bbox_inches='tight')
plt.show()
print("✅ 合并图已保存：combined_qzs_vibration_plot.png")