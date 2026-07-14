# ================================= 核心库导入 =================================
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from scipy.optimize import root
import warnings
warnings.filterwarnings('ignore')

# ================================= 全局绘图设置（Calibri字体+样式统一）=================================
plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Calibri', 'Microsoft YaHei', 'DejaVu Sans', 'Arial']
plt.rcParams['axes.titleweight'] = 'bold'   # 标题加粗
plt.rcParams['axes.labelweight'] = 'bold'   # 轴标签加粗
plt.rcParams['figure.facecolor'] = 'white'  # 画布白底
plt.rcParams['grid.alpha'] = 0.5            # 网格透明度
plt.rcParams['grid.linestyle'] = '--'       # 网格线型
plt.rcParams['legend.framealpha'] = 0.9     # 图例透明度
plt.rcParams['savefig.dpi'] = 300           # 保存图片高清DPI

# ================================= 1. 系统参数定义（与MATLAB1:1，SI单位）=================================
# 几何参数
b = 80e-3;          # 基圆半径 (m)
l0 = 80e-3;         # 弹簧自然长度 (m)，QZS满足l0=b
x0 = 5e-3;          # 初始偏移 (m)
g = 9.81;           # 重力加速度 (m/s²)

# 弹簧刚度
k = 300;            # 单个弹簧刚度 (N/m)

# 阻尼参数
zeta = 0.0311;      # 阻尼比（实验值）

# 激励参数
W = 2e-3;           # 基础激励振幅 (m)
freq_range = [0, 90];# 频率范围 (Hz)
n_freq = 200;       # 频率点数

# ================================= 2. 静态平衡分析 =================================
# 初始角度（无载荷位置）
phi_0 = np.arcsin((b - x0) / b)
# 角度范围
phi_array = np.linspace(phi_0, np.arcsin((b - 50e-3) / b), 200)
# 从初始位置的位移
z_disp = 2 * b * (np.sin(phi_0) - np.sin(phi_array))
# 竖向恢复力（手稿Eq.18）
F_Ver = (9 * k * b / 2) * np.sin(phi_array) * (np.cos(phi_array) - np.cos(phi_0)) / np.cos(phi_array)

# 有效载荷对应的平衡位置（取最大恢复力点）
F_max = np.max(F_Ver)
idx_eq = np.argmax(F_Ver)
phi_eq = phi_array[idx_eq]
z_eq = z_disp[idx_eq]
# 等效质量
M = F_max / g

# 打印系统参数
print('='*40 + ' SYSTEM PARAMETERS ' + '='*40)
print(f'  平衡角度 φ_eq = {phi_eq:.4f} rad ({phi_eq*180/np.pi:.2f}°)')
print(f'  支撑质量 M = {M:.4f} kg')
print(f'  平衡位置 z_eq = {z_eq*1000:.4f} mm')

# ================================= 3. 线性化刚度系数（手稿Eq.24-25）=================================
k1 = (9 * k / 2) * (1 - np.cos(phi_0) / np.cos(phi_eq))
k3 = (9 * k * np.cos(phi_0)) / (4 * b**2 * np.cos(phi_eq)**3)
# 固有频率（abs(k1)避免k1<0）
omega_n = np.sqrt(np.abs(k1) / M)
f_n = omega_n / (2 * np.pi)
# 阻尼系数（围绕k1线性化）
c = 2 * zeta * np.sqrt(np.abs(k1) * M)

# 打印线性化参数
print(f'  线性刚度 k1 = {k1:.4f} N/m')
print(f'  非线性刚度 k3 = {k3:.4e} N/m^3')
print(f'  固有频率 f_n = {f_n:.4f} Hz')
print(f'  阻尼系数 c = {c:.4f} N·s/m')
print('='*90)

# ================================= 4. HB法传递率计算 =================================
print('\nCalculating transmissibility via Harmonic Balance Method...')
# 频率数组（避免0Hz除零）
f_values = np.linspace(max(freq_range[0], 1e-6), freq_range[1], n_freq)
omega_values = 2 * np.pi * f_values
# 初始化数组
Tr_values_dB = np.zeros_like(f_values)  # 位移传递率(dB)
Tr_force_dB = np.zeros_like(f_values)   # 力传递率(dB)
U_values = np.zeros_like(f_values)      # 相对振幅
phi_values = np.zeros_like(f_values)    # 相位滞后

# Root求解器配置（替代MATLAB fsolve，彻底避坑）
root_options = {
    'method': 'hybr',  # 与fsolve默认算法一致
    'tol': 1e-12,      # 全局精度
    'options': {'disp': False}
}
U_prev = W  # 连续法初始猜测

# 频率循环
for i, omega in enumerate(omega_values):
    # 定义HB振幅方程
    def HB_eq(U):
        return ((k1 - M*omega**2 + 0.75*k3*U**2)**2 + (c*omega)**2) * U**2 - (M*omega**2*W)**2
    
    # 初始猜测：连续值+线性解
    U_linear = np.abs(M*omega**2*W) / np.sqrt((k1 - M*omega**2)**2 + (c*omega)**2)
    U0 = max([U_prev, U_linear, 1e-9])
    
    # Root求解
    res = root(HB_eq, x0=U0, **root_options)
    U_sol = res.x[0] if res.success else np.nan
    # 求解失败则尝试备选猜测
    if not res.success or not np.isfinite(U_sol):
        guess_set = [U_linear, W, 0.1*W, 10*W, U0]
        U_sol = np.nan
        for gg in guess_set:
            res_tmp = root(HB_eq, x0=max(gg, 1e-9), **root_options)
            if res_tmp.success and np.isfinite(res_tmp.x[0]):
                U_sol = res_tmp.x[0]
                break
        if np.isnan(U_sol):
            U_sol = U0  # 兜底避免崩溃
    
    # 保存结果
    U = np.abs(U_sol)
    U_values[i] = U
    U_prev = U
    # 有效刚度
    k_eff = k1 + 0.75 * k3 * U**2
    # 相位：u(t)=U cos(ωt-φ)
    phi = np.arctan2(c*omega, (k_eff - M*omega**2))
    phi_values[i] = phi
    # 力传递率（可选）
    Tr_force = np.sqrt(k_eff**2 + (c*omega)**2) / np.sqrt((k_eff - M*omega**2)**2 + (c*omega)**2)
    Tr_force_dB[i] = 20 * np.log10(Tr_force)
    # 位移传递率（dB）
    Tr_disp = np.sqrt(1 + (U/W)**2 + 2*(U/W)*np.cos(phi))
    Tr_values_dB[i] = 20 * np.log10(Tr_disp)

print('Transmissibility calculation complete!')

# ================================= 5. 绘图1：0-90Hz 传递率曲线 =================================
fig1, ax1 = plt.subplots(1,1,figsize=(9,5.2), dpi=100)
# 绘制理论曲线
ax1.plot(f_values, Tr_values_dB, 'b-', linewidth=2, label='Theory (HB, disp. Tr)')
# 0dB参考线
ax1.axhline(y=0, linestyle='--', color='k', linewidth=1, label='0 dB')
# 共振峰标记
Tr_max_dB, idx_res = np.max(Tr_values_dB), np.argmax(Tr_values_dB)
f_res = f_values[idx_res]
ax1.plot(f_res, Tr_max_dB, 'ro', markersize=9, markerfacecolor='r', 
         label=f'Peak: {f_res:.1f} Hz')
# 隔振起始频率（首次<0dB）
idx_isolation = np.argmax(Tr_values_dB < 0) if np.any(Tr_values_dB < 0) else None
if idx_isolation is not None and idx_isolation < len(f_values):
    f_isolation = f_values[idx_isolation]
    ax1.axvline(x=f_isolation, linestyle='--', color='g', linewidth=1.5, 
                label=f'Isolation onset: {f_isolation:.1f} Hz')

# 样式设置（字体+4）
ax1.set_xlabel('Frequency (Hz)', fontsize=18, fontweight='bold')
ax1.set_ylabel('Transmissibility (dB)', fontsize=18, fontweight='bold')
ax1.set_title('Transmissibility Curve of QZS Vibration Isolator (0–90 Hz)', fontsize=20, fontweight='bold')
ax1.set_xlim(0,90)
ax1.set_ylim(-60,20)
ax1.grid(True, which='both')
ax1.legend(loc='upper right', fontsize=15)
ax1.tick_params(axis='both', labelsize=16)
plt.tight_layout()

# ================================= 6. 绘图2：三参数组合传递率对比 =================================
# 定义三参数集（与MATLAB1:1）
param_sets = [
    {'b':60e-3, 'k':80,  'label':'b=60 mm, k=80 N/m',  'color':'#E74C3C', 'linestyle':'-'},
    {'b':80e-3, 'k':100, 'label':'b=80 mm, k=100 N/m', 'color':'#2ECC71', 'linestyle':'--'},
    {'b':90e-3, 'k':120, 'label':'b=90 mm, k=120 N/m', 'color':'#9B59B6', 'linestyle':':'}
]
# 共享参数
x0_share = 5e-3; zeta_share = 0.0311; W_share = 2e-3; g_share = 9.81
f_vals_share = np.linspace(max(freq_range[0],1e-6), freq_range[1], n_freq)
omega_vals_share = 2*np.pi*f_vals_share

# 创建画布
fig2, ax2 = plt.subplots(1,1,figsize=(9,5.2), dpi=100)
ax2.grid(True, which='both')

# 循环参数集
for params in param_sets:
    b_p, k_p, label_p, color_p, ls_p = params['b'], params['k'], params['label'], params['color'], params['linestyle']
    # 静态平衡
    phi0_p = np.arcsin((b_p - x0_share)/b_p)
    phi_arr_p = np.linspace(phi0_p, np.arcsin((b_p -50e-3)/b_p), 200)
    F_ver_p = (9*k_p*b_p/2)*np.sin(phi_arr_p)*(np.cos(phi_arr_p)-np.cos(phi0_p))/np.cos(phi_arr_p)
    F_max_p = np.max(F_ver_p)
    phi_eq_p = phi_arr_p[np.argmax(F_ver_p)]
    M_p = F_max_p / g_share
    # 线性化刚度
    k1_p = (9*k_p/2)*(1 - np.cos(phi0_p)/np.cos(phi_eq_p))
    k3_p = (9*k_p*np.cos(phi0_p))/(4*b_p**2*np.cos(phi_eq_p)**3)
    c_p = 2*zeta_share*np.sqrt(np.abs(k1_p)*M_p)
    # HB传递率计算
    Tr_dB_p = np.zeros_like(f_vals_share)
    U_prev_p = W_share
    for i, omega in enumerate(omega_vals_share):
        def HB_eq_p(U):
            return ((k1_p - M_p*omega**2 +0.75*k3_p*U**2)**2 + (c_p*omega)**2)*U**2 - (M_p*omega**2*W_share)**2
        U_lin_p = np.abs(M_p*omega**2*W_share)/np.sqrt((k1_p-M_p*omega**2)**2 + (c_p*omega)**2)
        U0_p = max([U_prev_p, U_lin_p, 1e-9])
        res_p = root(HB_eq_p, x0=U0_p, **root_options)
        U_sol_p = res_p.x[0] if res_p.success else U0_p
        U_p = np.abs(U_sol_p)
        U_prev_p = U_p
        k_eff_p = k1_p +0.75*k3_p*U_p**2
        phi_p = np.arctan2(c_p*omega, k_eff_p - M_p*omega**2)
        Tr_disp_p = np.sqrt(1 + (U_p/W_share)**2 + 2*(U_p/W_share)*np.cos(phi_p))
        Tr_dB_p[i] = 20*np.log10(Tr_disp_p)
    # 绘制曲线
    ax2.plot(f_vals_share, Tr_dB_p, color=color_p, linestyle=ls_p, linewidth=2, label=label_p)
    # 标记共振峰
    Tr_max_p, idx_res_p = np.max(Tr_dB_p), np.argmax(Tr_dB_p)
    f_res_p = f_vals_share[idx_res_p]
    ax2.plot(f_res_p, Tr_max_p, 'o', color=color_p, markersize=8, markerfacecolor=color_p)

# 0dB参考线
ax2.axhline(y=0, linestyle='--', color='k', linewidth=1, label='0 dB (Isolation Threshold)')
# 样式（字体+4）
ax2.set_xlabel('Frequency (Hz)', fontsize=18, fontweight='bold')
ax2.set_ylabel('Transmissibility (dB)', fontsize=18, fontweight='bold')
ax2.set_title('Transmissibility Curves of QZS Vibration Isolator with Different Parameters', fontsize=20, fontweight='bold')
ax2.set_xlim(0,90); ax2.set_ylim(-60,20)
ax2.legend(loc='upper right', fontsize=14)
ax2.tick_params(axis='both', labelsize=16)
plt.tight_layout()


# ================================= 7. 绘图3：多频率稳态响应（位移+加速度 3行2列）=================================
# 分析频率：共振峰+10Hz+20Hz
freq_list = [f_res, 10, 20]
freq_labels = [f'Resonance Frequency ({f_res:.2f} Hz)', '10 Hz', '20 Hz']
colors = ['#2980B9', '#E67E22']  # 蓝色(top)/橙色(base)
default_color = 'k'

# 全局统一字体大小（与之前+4设置一致，便于维护）
title_size = 18
label_size = 16
tick_size = 15
legend_size = 14

# 创建3行2列画布，加宽适配双列
fig3, axs3 = plt.subplots(3,2,figsize=(12,8), dpi=100)
# fig3.suptitle('Steady-State Displacement & Acceleration Response at Different Frequencies', 
#              fontsize=20, fontweight='bold', y=0.98)

# 循环频率
for i, f_time in enumerate(freq_list):
    omega = 2 * np.pi * f_time
    # 获取对应U和phi
    idx_closest = np.argmin(np.abs(f_values - f_time))
    U = U_values[idx_closest]
    phi = phi_values[idx_closest]
    # 时间向量：20个周期，高采样率
    N_cycles = 20; T = 1/f_time; Fs = max(5000, 2000*f_time); dt = 1/Fs
    t = np.linspace(0, N_cycles*T, int(N_cycles*T*Fs))
    
    # 计算位移（mm）和加速度（m/s²）
    z_b = W * np.cos(omega*t) * 1000  # 基础位移(mm)
    u = U * np.cos(omega*t - phi) * 1000  # 相对位移(mm)
    z_top = u + z_b  # 顶部位移(mm)
    a_b = -omega**2 * (z_b / 1000)  # 基础加速度(m/s²)
    a_u = -omega**2 * (u / 1000)  # 相对加速度(m/s²)
    a_top = a_u + a_b  # 顶部加速度(m/s²)
    
    # ---------------------- 左列：位移-时间子图 (3行2列，第i行1列) ----------------------
    ax_disp = axs3[i,0]
    ax_disp.plot(t, z_top, color=colors[0], linewidth=1.8, label='Top Displacement (z_{top})')
    ax_disp.plot(t, z_b, color=colors[1], linestyle='--', linewidth=1.2, label='Base Displacement (z_{base})')
    # 标题/轴标签（字体大小全局统一）
    ax_disp.set_title(f'{freq_labels[i]} - Displacement', fontsize=title_size, fontweight='bold')
    ax_disp.set_xlabel('Time (s)', fontsize=label_size, fontweight='bold')
    ax_disp.set_ylabel('Displacement (mm)', fontsize=label_size, fontweight='bold')
    # 图例（右上角+统一字体）
    ax_disp.legend(loc='upper right', fontsize=legend_size)
    ax_disp.grid(True, which='both')
    
    # 1. 强制位移轴包含0点，锁定Y轴范围
    ylim_d = ax_disp.get_ylim()
    if ylim_d[0] > 0:
        ylim_unify = (0, ylim_d[1])
    elif ylim_d[1] < 0:
        ylim_unify = (ylim_d[0], 0)
    else:
        ylim_unify = ylim_d  # 本身包含0点，直接用原范围
    ax_disp.set_ylim(ylim_unify)  # 锁定Y轴范围
    
    # 2. 左侧Y轴样式（蓝色+字体统一）
    ax_disp.spines['left'].set_color(colors[0])
    ax_disp.tick_params(axis='both', labelsize=tick_size)  # 刻度标注字体统一
    ax_disp.tick_params(axis='y', colors=colors[0])       # 左侧Y轴刻度颜色匹配
    ax_disp.yaxis.label.set_color(colors[0])              # 左侧Y轴标签颜色匹配
    
    # 3. 右侧Y轴 - 核心：共用刻度+仅改标签颜色
    ax_disp2 = ax_disp.twinx()  # 创建右侧双Y轴
    # 强制右侧Y轴范围与左侧完全一致（共用刻度）
    ax_disp2.set_ylim(ylim_unify)
    # 右侧Y轴标签：内容与左侧相同，仅颜色改为橙色
    ax_disp2.set_ylabel('Displacement (mm)', fontsize=label_size, fontweight='bold')
    # 右侧Y轴样式：仅显示轴线+标签，隐藏刻度数字（共用左侧刻度）
    ax_disp2.spines['right'].set_color(colors[1])
    ax_disp2.tick_params(axis='y', colors=colors[1], labelsize=tick_size)  
    ax_disp2.yaxis.label.set_color(colors[1])                                                # 右侧Y轴标签颜色匹配


    # ---------------------- 右列：加速度-时间子图 (3行2列，第i行2列) ----------------------
    ax_acc = axs3[i,1]
    # 直接绘图（样式参数与左列完全一致）
    ax_acc.plot(t, a_top, color=colors[0], linewidth=1.8, label='Top Acceleration (a_{top})')
    ax_acc.plot(t, a_b, color=colors[1], linestyle='--', linewidth=1.2, label='Base Acceleration (a_{base})')
    # 标题/轴标签（字体大小全局统一，与左列参数完全一致）
    ax_acc.set_title(f'{freq_labels[i]} - Acceleration', fontsize=title_size, fontweight='bold')
    ax_acc.set_xlabel('Time (s)', fontsize=label_size, fontweight='bold')
    ax_acc.set_ylabel('Acceleration (m/s²)', fontsize=label_size, fontweight='bold')
    # 图例：完全参考左列 - loc/字体大小1:1复刻
    ax_acc.legend(loc='upper right', fontsize=legend_size)
    ax_acc.grid(True, which='both')
    
    # 1. 强制加速度轴包含0点，锁定Y轴范围
    ylim_a = ax_acc.get_ylim()
    if ylim_a[0]>0:
        ylim_unify_a = (0, ylim_a[1])
    elif ylim_a[1]<0:
        ylim_unify_a = (ylim_a[0], 0)
    else:
        ylim_unify_a = ylim_a
    ax_acc.set_ylim(ylim_unify_a)  # 锁定Y轴范围
    
    # 2. 左侧Y轴样式（蓝色+字体统一，与左列逻辑一致）
    ax_acc.spines['left'].set_color(colors[0])
    ax_acc.tick_params(axis='both', labelsize=tick_size)
    ax_acc.tick_params(axis='y', colors=colors[0])
    ax_acc.yaxis.label.set_color(colors[0])
    
    # 3. 右侧Y轴 - 核心：共用刻度（取消g换算）+仅改标签颜色
    ax_acc2 = ax_acc.twinx()
    # 强制右侧Y轴范围与左侧完全一致（共用刻度，无g换算）
    ax_acc2.set_ylim(ylim_unify_a)
    # 右侧Y轴标签：内容与左侧相同，仅颜色改为橙色
    ax_acc2.set_ylabel('Acceleration (m/s²)', fontsize=label_size, fontweight='bold')
    # 右侧Y轴样式：仅显示轴线+标签，隐藏刻度数字（共用左侧刻度）
    ax_acc2.spines['right'].set_color(colors[1])
    ax_acc2.tick_params(axis='y', colors=colors[1], labelsize=tick_size)  
    ax_acc2.yaxis.label.set_color(colors[1])                                                # 右侧Y轴标签颜色匹配

# 调整子图间距，避免标签重叠
plt.tight_layout()
plt.subplots_adjust(top=0.93, hspace=0.3, wspace=0.35)

# ================================= 8. 绘图4：理论VS实验传递率对比 =================================
# ---------------------- 实验数据配置（修改为你的CSV实际路径！）----------------------
data_folder = 'Vib_data_2026'  # CSV文件夹
exp_file = 'Data_0129-1.csv'   # 实验CSV文件
exp_path = os.path.join(data_folder, exp_file)

# 读取并预处理实验数据
if os.path.exists(exp_path):
    exp_data = pd.read_csv(exp_path, header=None).values
    exp_data = exp_data[np.all(np.isfinite(exp_data), axis=1)]  # 过滤NaN
    freq_exp = exp_data[:,0] - 1  # 保留MATLAB的-1修正
    Tr_exp = exp_data[:,1]        # 取第二列（解决索引越界）
    print(f'\n实验数据读取成功：{len(freq_exp)} 个原始点')
    
    # 频率取整+去重平均
    freq_exp_round = np.round(freq_exp).astype(int)
    unique_freq, idxu = np.unique(freq_exp_round, return_inverse=True)
    unique_Tr = np.array([np.mean(Tr_exp[idxu==i]) for i in range(len(unique_freq))])
    # 过滤0-90Hz
    valid = (unique_freq >=0) & (unique_freq <=90)
    freq_exp, Tr_exp = unique_freq[valid], unique_Tr[valid]
    print(f'实验数据预处理后：{len(freq_exp)} 个有效点 (0-90Hz)')
else:
    raise FileNotFoundError(f'实验文件未找到：{exp_path}，请检查路径！')

# 创建对比画布
fig4, ax4 = plt.subplots(1,1,figsize=(9.5,5.5), dpi=100)
ax4.grid(True, which='both')
# 绘制理论曲线
ax4.plot(f_values, Tr_values_dB, 'b-', linewidth=2.2, label='Theory (HB, disp. Tr)')
# 绘制实验数据（红色点线图，与理论色完全区分）
ax4.plot(freq_exp, Tr_exp, 'ro-', linewidth=1.5, markersize=6, markerfacecolor='r', label='Experimental Data')
# 0dB参考线
ax4.axhline(y=0, linestyle='--', color='k', linewidth=1, label='0 dB')
# 理论共振峰（蓝色方形）
ax4.plot(f_res, Tr_max_dB, 'bs', markersize=8, markerfacecolor='b', label=f'Theory Peak: {f_res:.1f} Hz')
# 隔振起始频率
if idx_isolation is not None and idx_isolation < len(f_values):
    ax4.axvline(x=f_isolation, linestyle='--', color='g', linewidth=1.5, label=f'Isolation Onset: {f_isolation:.1f} Hz')

# 样式（字体+4）
ax4.set_xlabel('Frequency (Hz)', fontsize=18, fontweight='bold')
ax4.set_ylabel('Transmissibility (dB)', fontsize=18, fontweight='bold')
ax4.set_title('Experimental vs Theoretical Transmissibility Comparison of QZS Vibration Isolator', fontsize=20, fontweight='bold')
ax4.set_xlim(0,90); ax4.set_ylim(-60,20)
ax4.legend(loc='upper right', fontsize=15)
ax4.tick_params(axis='both', labelsize=16)
plt.tight_layout()

# ================================= 9. 绘图5：静态参考图（力-位移+局部刚度）=================================
# 9.1 静态力-位移曲线
fig5_1, ax5_1 = plt.subplots(1,1,figsize=(8.5,5.2), dpi=100)
ax5_1.plot(z_disp*1000, F_Ver, 'b-', linewidth=2, label='Force-Displacement')
ax5_1.plot(z_eq*1000, F_max, 'ro', markersize=10, markerfacecolor='r', label='Equilibrium')
ax5_1.set_xlabel('Displacement (mm)', fontsize=18, fontweight='bold')
ax5_1.set_ylabel('Restoring Force (N)', fontsize=18, fontweight='bold')
ax5_1.set_title('Static Force-Displacement Characteristic', fontsize=20, fontweight='bold')
ax5_1.legend(loc='upper right', fontsize=16)
ax5_1.grid(True, which='both')
ax5_1.tick_params(axis='both', labelsize=16)
plt.tight_layout()

# 9.2 局部刚度曲线
dz = np.diff(z_disp)
dF = np.diff(F_Ver)
k_local = dF / dz
z_mid = (z_disp[:-1] + z_disp[1:]) / 2
fig5_2, ax5_2 = plt.subplots(1,1,figsize=(8.5,4.2), dpi=100)
ax5_2.plot(z_mid*1000, k_local, 'b-', linewidth=1.5)
ax5_2.axhline(y=0, linestyle='--', color='k', linewidth=1)
QZS_threshold = 0.1*k
ax5_2.axhline(y=QZS_threshold, linestyle=':', color='r', linewidth=1)
ax5_2.axhline(y=-QZS_threshold, linestyle=':', color='r', linewidth=1)
ax5_2.set_xlabel('Displacement (mm)', fontsize=18, fontweight='bold')
ax5_2.set_ylabel('Local Stiffness (N/m)', fontsize=18, fontweight='bold')
ax5_2.set_title('Local Stiffness vs Displacement', fontsize=20, fontweight='bold')
ax5_2.grid(True, which='both')
ax5_2.tick_params(axis='both', labelsize=16)
plt.tight_layout()

# ================================= 10. 分析总结输出 =================================
print('\n' + '='*40 + ' ANALYSIS SUMMARY ' + '='*40)
print('Geometric Parameters:')
print(f'  b = {b*1000:.1f} mm, l0 = {l0*1000:.1f} mm')
print(f'  Initial angle φ_0 = {phi_0*180/np.pi:.2f}°')
print('\nDynamic Properties:')
print(f'  Supported mass M = {M:.3f} kg')
print(f'  Linear stiffness k1 = {k1:.2f} N/m')
print(f'  Nonlinear stiffness k3 = {k3:.2e} N/m^3')
print(f'  Natural frequency = {f_n:.2f} Hz')
print(f'  Damping ratio ζ = {zeta:.4f}')
print('\nIsolation Performance (Disp. Tr):')
print(f'  Peak transmissibility = {Tr_max_dB:.2f} dB at {f_res:.1f} Hz')
if idx_isolation is not None and idx_isolation < len(f_values):
    print(f'  Isolation onset frequency = {f_isolation:.1f} Hz')
print('='*90)

# 显示所有图片
plt.show()



# 保存所有图片（可选，取消注释即可）
# fig1.savefig('QZS_Tr_Curve_0-90Hz.png', bbox_inches='tight')
# fig2.savefig('QZS_Tr_3Params_Comparison.png', bbox_inches='tight')
# fig3.savefig('QZS_SteadyState_Disp_Acc.png', bbox_inches='tight')
# fig4.savefig('QZS_Theory_vs_Experiment.png', bbox_inches='tight')
# fig5_1.savefig('QZS_Static_Force_Displacement.png', bbox_inches='tight')
# fig5_2.savefig('QZS_Local_Stiffness.png', bbox_inches='tight')
# print('\n所有图片保存成功！')