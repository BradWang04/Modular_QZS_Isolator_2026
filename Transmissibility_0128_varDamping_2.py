#!/usr/bin/env python3
import numpy as np
from scipy.optimize import root
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')  # 屏蔽求解警告

# ====================== 1. Common System Parameters 公共系统参数（与MATLAB一致）=====================
M = 0.548               # 质量 (kg)
zeta = 0.03             # 阻尼比
W = 2e-3                # 基础激励幅值 (m)
freq_range = [0, 10]    # 频率范围
n_freq = 500            # 频率点数

# ====================== 2. Define System Properties 定义三类系统参数（替代MATLAB struct）=====================
# 用字典列表存储每个系统的参数，对应MATLAB的struct数组
systems = [
    # Case 1: 线性系统
    {
        'name': 'Linear System',
        'k1': M * (2 * np.pi * 3.0)**2,  # 线性固有频率3Hz
        'k3': 0,
        'color': (0.5, 0.5, 0.5),        # 灰色
        'linestyle': ':'
    },
    # Case 2: 硬化QZS系统（k3>0）
    {
        'name': 'Hardening QZS (k3 > 0)',
        'k1': M * (2 * np.pi * 0.5)**2,  # QZS低固有频率0.5Hz
        'k3': 2.7e5,                     # 正非线性刚度
        'color': 'b',
        'linestyle': '-'
    },
    # Case 3: 软化QZS系统（k3<0）
    {
        'name': 'Softening QZS (k3 < 0)',
        'k1': M * (2 * np.pi * 0.5)**2,  # 同硬化QZS的k1
        'k3': -2.7e5,                    # 负非线性刚度
        'color': 'r',
        'linestyle': '--'
    }
]

# ====================== 3. Dynamic Analysis Loop 动力学分析循环=====================
# 创建画布（匹配MATLAB的Position [100,100,900,520] → figsize=(9,5.2)）
plt.figure(figsize=(9, 5.2))

# 循环每个系统计算传递率并绘图
for s_idx, current_system in enumerate(systems):
    print(f'\n--- Calculating for: {current_system["name"]} ---')
    
    k1 = current_system['k1']
    k3 = current_system['k3']
    
    # 阻尼系数（基于刚度线性部分）
    c = 2 * zeta * np.sqrt(np.abs(k1) * M)
    
    # 频率数组（避开0Hz奇点，对应MATLAB linspace(1e-6, ...)）
    f_values = np.linspace(1e-6, freq_range[1], n_freq)
    omega_values = 2 * np.pi * f_values
    Tr_values_dB = np.zeros_like(f_values)  # 传递率dB数组
    
    U_prev = W  # 迭代初值，保证收敛
    
    # 定义谐波平衡方程（移到循环外提升效率）
    def HB_eq(U, omega):
        return ((k1 - M*omega**2 + 0.75*k3*U**2)**2 + (c*omega)**2) * U**2 - (M*omega**2*W)**2
    
    # 频率循环求解传递率
    for i, omega in enumerate(omega_values):
        # 线性解作为兜底初值
        U_linear = np.abs(M*omega**2 * W) / np.sqrt((k1 - M*omega**2)**2 + (c*omega)**2 + 1e-12)
        U0 = max(U_prev, U_linear, 1e-9)  # 优化初值
        
        # 用root替代MATLAB fsolve（method='hybr'匹配fsolve默认方法）
        res = root(HB_eq, U0, args=(omega,), method='hybr', options={
            'xtol': 1e-12,
            'maxfev': 200,
            'disp': False
        })
        
        # 收敛性判断，不收敛则用线性解
        if not res.success or not np.isfinite(res.x[0]) or res.x[0] < 0:
            U_sol = U_linear
        else:
            U_sol = res.x[0]
        
        U = np.abs(U_sol)
        U_prev = U
        
        # 计算有效刚度、相位差、传递率dB
        k_eff = k1 + 0.75 * k3 * U**2
        phi = np.arctan2(c * omega, k_eff - M*omega**2)
        Tr_disp = np.sqrt(1 + (U/W)**2 + 2*(U/W)*np.cos(phi))
        Tr_values_dB[i] = 20 * np.log10(Tr_disp)
    
    # 绘制当前系统的传递率曲线（匹配MATLAB样式）
    plt.plot(f_values, Tr_values_dB, 
             color=current_system['color'],
             linestyle=current_system['linestyle'],
             linewidth=2.5,
             label=current_system['name'])

# ====================== 4. Finalize Plot 绘图收尾（与MATLAB完全一致）=====================
plt.axhline(y=0, color='k', linestyle='-', linewidth=1)  # 0dB水平线
plt.xlabel('Frequency (Hz)', fontsize=14, fontweight='bold')
plt.ylabel('Transmissibility (dB)', fontsize=14, fontweight='bold')
plt.title('Comparison of Linear, Hardening, and Softening Isolators', 
          fontsize=16, fontweight='bold')
plt.xlim(freq_range)
plt.ylim([-40, 30])
plt.grid(True, which='both')  # 开启主/次网格（对应MATLAB grid on; grid minor）
plt.legend(loc='upper right', fontsize=12)
plt.box(True)  # 对应MATLAB box on

# 显示图像
plt.show()