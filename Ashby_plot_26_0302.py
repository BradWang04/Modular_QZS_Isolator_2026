import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Calibri', 'Microsoft YaHei', 'DejaVu Sans', 'Arial']
plt.close('all')

studies = [
    ("Huang et al. [45]",  4.99,  166.7,  2.0,  15.0),
    ("Lan et al. [9]",     2.3,    58.9,  2.3,   7.0),
    ("Xu et al. [20]",     1.5,    57.6,  0.5,  20.0),
    ("Le & Ahn [42]",      0.6,    49.0,  0.1,  10.0),
    ("Zheng et al. [38]",  5.0,    36.3,  5.0,  25.0),  # capped at xlim in plot
    ("Yan et al. [2]",     8.76,    6.1,  None, None),
    ("Wang et al. [44]",   2.5,    98.0,  2.5,  10.0),
    ("Hu & Zhou [3]",      3.9,   148.0,  0.1,  20.0),
]

# ── FIX 3: unique colour for every study ──
colors_map = {
    "Huang et al. [45]":  '#d92523',
    "Lan et al. [9]":     '#2e7ebb',
    "Xu et al. [20]":     '#e25508',
    "Le & Ahn [42]":      '#2e974e',
    "Zheng et al. [38]":  '#7262ac',
    "Yan et al. [2]":     '#888888',
    "Wang et al. [44]":   '#c2880a',
    "Hu & Zhou [3]":      '#c4469e',   # FIX 3: was '#d92523', now distinct magenta
}

# ── FIX 1 & 5: (dx, dy, ha) — ha set correctly for the sign of dx ──
label_offsets = {
    "Huang et al. [45]":  ( 0.3,   6,  'left'),
    "Lan et al. [9]":     ( 0.3, -18,  'left'),
    "Xu et al. [20]":     ( 0.3,   6,  'left'),
    "Le & Ahn [42]":      (-0.3, -18,  'right'),  # FIX 1 & 5: ha='right' for negative dx
    "Zheng et al. [38]":  ( 0.3, -18,  'left'),
    "Yan et al. [2]":     (-0.3,   4,  'right'),  # FIX 1: ha='right' for negative dx
    "Wang et al. [44]":   ( 0.3,   6,  'left'),
    "Hu & Zhou [3]":      (-0.3, -18,  'right'),  # FIX 1: ha='right' for negative dx
}

def draw_ibeam(ax, fx, fy, col, lw=2.5, ms=9):
    """
    Draw a vertical I-beam marker scaled to the current axis limits.
    FIX 2: cap dimensions derived from axis range, not hardcoded.
    """
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    xrange = xlim[1] - xlim[0]
    yrange = ylim[1] - ylim[0]

    hw = xrange * 0.012      # horizontal cap half-width = 1.2 % of x-range
    hh = yrange * 0.025      # vertical stem half-height  = 2.5 % of y-range
    cap_h = yrange * 0.008   # cap rectangle half-height  = 0.8 % of y-range

    # Vertical stem
    ax.plot([fx, fx], [fy - hh, fy + hh], linewidth=lw, color=col)

    # Top cap (full width, centred on fx)
    ax.fill([fx - hw, fx - hw, fx + hw, fx + hw],
            [fy + hh - cap_h, fy + hh + cap_h,
             fy + hh + cap_h, fy + hh - cap_h],
            color=col, edgecolor=col, linewidth=lw)

    # Bottom cap
    ax.fill([fx - hw, fx - hw, fx + hw, fx + hw],
            [fy - hh - cap_h, fy - hh + cap_h,
             fy - hh + cap_h, fy - hh - cap_h],
            color=col, edgecolor=col, linewidth=lw)

    # Centre dot
    ax.plot(fx, fy, marker='o', markersize=ms, color=col, linestyle='none')


# ── Figure ──
fig2, ax2 = plt.subplots(figsize=(13, 8.5))

ax2.set_xlabel('Minimum Isolation Frequency (Hz)', fontsize=26, fontweight='bold')
ax2.set_ylabel('Loading (N)', fontsize=26, fontweight='bold')
ax2.set_title('QZS Vibration Isolator — Verified Experimental Data\n(8 Studies)',
              fontsize=22, fontweight='bold', pad=14)

ax2.set_xlim([0, 22])
ax2.set_ylim([0, 200])
ax2.grid(True, linestyle='-', alpha=0.7)
ax2.spines['top'].set_visible(True)
ax2.spines['right'].set_visible(True)

for (label, f_iso, loading, f_low, f_high) in studies:
    col = colors_map[label]

    # Dashed horizontal range bar
    if f_low is not None and f_high is not None:
        f_high_clipped = min(f_high, ax2.get_xlim()[1])   # clip to xlim
        ax2.plot([f_low, f_high_clipped], [loading, loading],
                 linewidth=2, color=col, alpha=0.35, linestyle='--')

    # I-beam marker (FIX 2: axis-scaled caps)
    draw_ibeam(ax2, f_iso, loading, col=col, lw=2.5, ms=9)

    # Label
    dx, dy, ha = label_offsets[label]
    ax2.text(f_iso + dx, loading + dy, label,
             ha=ha, va='bottom', fontsize=16, fontweight='bold', color=col)

# Annotation for shock isolator
ax2.annotate('(shock isolator,\nno f_iso defined)',
             xy=(8.76, 6.1), xytext=(10.5, 22),
             fontsize=13, color='#888888', style='italic',
             arrowprops=dict(arrowstyle='->', color='#888888', lw=1.2))

# Exclusion note (FIX 4: manual line break instead of wrap=True)
ax2.text(0.4, 196,
         'Note: Carrella et al. (2007) and Ye et al. (2021) excluded\n'
         '(purely theoretical; no experimental loading or frequency data).',
         ha='left', va='top', fontsize=12, color='#555555', style='italic',
         bbox=dict(boxstyle='round,pad=0.4', facecolor='#f5f5f5',
                   edgecolor='#cccccc', alpha=0.85))

# Legend
legend_entries = [
    mpatches.Patch(color='#d92523', label='Euler buckled beam QZS (Huang et al.)'),
    mpatches.Patch(color='#2e7ebb', label='Cam-roller-spring QZS (Lan et al.)'),
    mpatches.Patch(color='#e25508', label='Magnetic spring QZS (Xu et al.)'),
    mpatches.Patch(color='#2e974e', label='Negative stiffness seat isolator (Le & Ahn)'),
    mpatches.Patch(color='#7262ac', label='Multi-layer QZS (Zheng et al.)'),
    mpatches.Patch(color='#888888', label='Bistable shock isolator (Yan et al.) — no f_iso'),
    mpatches.Patch(color='#c2880a', label='Dual QZS ultra-low-frequency (Wang et al.)'),
    mpatches.Patch(color='#c4469e', label='QZS with quadratic stiffness (Hu & Zhou)'),  # FIX 3
]
ax2.legend(handles=legend_entries, loc='upper right', fontsize=11,
           framealpha=0.92, edgecolor='#aaaaaa',
           title='Isolator type', title_fontsize=12)

# Caption (FIX 4: explicit \n breaks instead of wrap=True)
caption = (
    "Figure: Verified operating points of eight QZS vibration isolators, plotted as minimum isolation\n"
    "frequency vs. static loading. I-beam markers show isolation onset; dashed bars show tested frequency range.\n"
    "Carrella et al. (2007) and Ye et al. (2021) omitted (theory only). Yan et al. (2022) flagged as shock isolator."
)
fig2.text(0.5, -0.03, caption, ha='center', va='top', fontsize=11,
          style='italic', color='#333333')

ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)

plt.tight_layout()
fig2.savefig('qzs_verified_plot.png', dpi=200, bbox_inches='tight', facecolor='white')
print("Plot saved.")
plt.show()