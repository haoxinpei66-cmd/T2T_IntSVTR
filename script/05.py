#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ==================== AI兼容性设置 ====================
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'SimHei']

plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 10,
    'figure.titlesize': 12,
    'figure.dpi': 300,
})

# 颜色定义（四种区域）
region_colors = {
    'SD': '#1f77b4',
    'SR': '#ff7f0e',
    'RM': '#2ca02c',
    'Unique': '#d62728'
}

# 区域顺序（固定）
region_order = ['SD', 'SR', 'RM', 'Unique']
# SV类型顺序
sv_types = ['DEL', 'INS', 'DUP', 'INV', 'Alu', 'LINE', 'LTR', 'SVA']

def load_data(stat_file):
    """读取统计文件"""
    df = pd.read_csv(stat_file, sep='\t')
    return df

def prepare_plot_data(df):
    """准备绘图数据：生成每个SV_TYPE下四个REGION的索引及数值"""
    x_labels = []
    x_positions = []
    idx = 0
    for sv in sv_types:
        for region in region_order:
            x_labels.append(f'{sv}\n{region}')
            x_positions.append(idx)
            idx += 1
    
    nw_total = []
    nw_subset = []
    kg_total = []
    kg_subset = []
    max_ratio = []
    
    for sv in sv_types:
        for region in region_order:
            row = df[(df['SV_TYPE'] == sv) & (df['REGION'] == region)]
            if len(row) == 0:
                nw_total.append(0); nw_subset.append(0)
                kg_total.append(0); kg_subset.append(0)
                max_ratio.append(0)
            else:
                row = row.iloc[0]
                nw_total.append(row['NyuWa_ALL'])
                nw_subset.append(row['NyuWa_SUBSET'])
                kg_total.append(row['1KG_ALL'])
                kg_subset.append(row['1KG_SUBSET'])
                max_ratio.append(max(row['NyuWa_RATIO'], row['1KG_RATIO']))
    
    return x_positions, x_labels, nw_total, nw_subset, kg_total, kg_subset, max_ratio

def create_three_panel_plot(x_positions, x_labels, nw_total, nw_subset, kg_total, kg_subset, max_ratio, output_file):
    """创建三合一图（只显示总数标签，不显示子集标签）"""
    all_counts = nw_total + kg_total
    y_max_count = max(all_counts) * 1.1 if max(all_counts) > 0 else 10
    
    fig, axes = plt.subplots(3, 1, figsize=(20, 12), sharex=True)
    ax1, ax2, ax3 = axes
    width = 0.8
    
    # ---------- 上子图：NyuWa绝对值 ----------
    for i, (total, subset) in enumerate(zip(nw_total, nw_subset)):
        region_name = x_labels[i].split('\n')[1]
        color = region_colors[region_name]
        ax1.bar(x_positions[i], total, width, color=color, alpha=0.6, edgecolor='none')
        ax1.bar(x_positions[i], subset, width, color=color, alpha=1.0, edgecolor='none')
        if total > 0:
            ax1.text(x_positions[i], total + y_max_count*0.01, f'{int(total):,}',
                     ha='center', va='bottom', fontsize=6, fontweight='bold', color='black')
    ax1.set_ylabel('Count', fontsize=10)
    ax1.set_ylim(0, y_max_count)
    ax1.set_title('')
    ax1.grid(True, alpha=0.3, axis='y')
    
    # ---------- 中子图：1KG绝对值 ----------
    for i, (total, subset) in enumerate(zip(kg_total, kg_subset)):
        region_name = x_labels[i].split('\n')[1]
        color = region_colors[region_name]
        ax2.bar(x_positions[i], total, width, color=color, alpha=0.6, edgecolor='none')
        ax2.bar(x_positions[i], subset, width, color=color, alpha=1.0, edgecolor='none')
        if total > 0:
            ax2.text(x_positions[i], total + y_max_count*0.01, f'{int(total):,}',
                     ha='center', va='bottom', fontsize=6, fontweight='bold', color='black')
    ax2.set_ylabel('Count', fontsize=10)
    ax2.set_ylim(0, y_max_count)
    ax2.set_title('')
    ax2.grid(True, alpha=0.3, axis='y')
    
    # ---------- 下子图：比例（取较大值），颜色改为浅绿色 lightgreen ----------
    for i, ratio in enumerate(max_ratio):
        ax3.bar(x_positions[i], ratio, width, color='lightgreen', alpha=0.8, edgecolor='none')
        if ratio > 0:
            ax3.text(x_positions[i], ratio + 0.02, f'{ratio:.2f}',
                     ha='center', va='bottom', fontsize=6, fontweight='bold', color='black')
    ax3.set_ylabel('Ratio (max of two)', fontsize=10)
    ax3.set_ylim(0, 1.05)
    ax3.set_title('')
    ax3.grid(True, alpha=0.3, axis='y')
    
    # 设置共同的x轴刻度
    ax3.set_xticks(x_positions)
    ax3.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=8)
    
    plt.subplots_adjust(hspace=0.1)
    plt.savefig(output_file, bbox_inches='tight', facecolor='white', format='pdf')
    print(f"三合一图已保存至: {output_file}")

def main():
    stat_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/1019SV_NWSV_NWME/1019sv_1019inv_nwsv_nwme_allsvtype_region_ratio_detail_over05_svlenge50.rm05mesv.rm05trsv_splitdup_new.txt"
    output_stat_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/1019SV_NWSV_NWME/1019sv_1019inv_nwsv_nwme_allsvtype_region_ratio_detail_over05_svlenge50.rm05mesv.rm05trsv_splitdup_new_.txt"
    output_plot_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/1019SV_NWSV_NWME/1019sv_1019inv_nwsv_nwme_allsvtype_region_ratio_distribution_over05_svlenge50.rm05mesv.rm05trsv_splitdup_new_addratio.pdf"
    
    print("读取统计文件...")
    df = load_data(stat_file)
    df.to_csv(output_stat_file, sep='\t', index=False)
    print(f"统计文件已复制至: {output_stat_file}")
    
    print("准备绘图数据...")
    x_positions, x_labels, nw_total, nw_subset, kg_total, kg_subset, max_ratio = prepare_plot_data(df)
    
    print("生成三合一图...")
    create_three_panel_plot(x_positions, x_labels, nw_total, nw_subset, kg_total, kg_subset, max_ratio, output_plot_file)
    
    print("完成！")

if __name__ == "__main__":
    main()
