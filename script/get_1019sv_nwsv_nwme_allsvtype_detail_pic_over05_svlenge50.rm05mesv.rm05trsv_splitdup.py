#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, Tuple, Set

# ==================== AI兼容性设置 ====================
def set_custom_style():
    plt.rcParams.update({
        'font.size': 24,
        'axes.titlesize': 20,
        'axes.labelsize': 28,
        'xtick.labelsize': 26,
        'ytick.labelsize': 26,
        'legend.fontsize': 24,
        'figure.titlesize': 22,
        'figure.dpi': 300,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'font.family': 'Arial',
        'axes.unicode_minus': False
    })

def read_1kg_all_variants(file_path: str) -> Dict[str, int]:
    """读取千人所有变异文件，统计变异类型分布"""
    print("正在读取千人所有变异文件...")
    svtype_count = {}
    try:
        df = pd.read_csv(file_path, sep='\t', comment='#', header=None,
                        names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'SV_TYPE', 'SVLEN'])
        svtype_counts = df['SV_TYPE'].value_counts()
        svtype_count = svtype_counts.to_dict()
        print(f"千人所有变异文件读取完成，共 {len(df)} 个变异")
        for sv_type, count in sorted(svtype_count.items(), key=lambda x: x[1], reverse=True):
            print(f"  {sv_type}: {count}")
    except Exception as e:
        print(f"错误：读取千人所有变异文件失败 - {e}")
    return svtype_count

def read_nyuwa_all_variants(file_path: str) -> Dict[str, int]:
    """读取NyuWa所有变异文件，统计变异类型分布"""
    print("正在读取NyuWa所有变异文件...")
    svtype_count = {}
    try:
        df = pd.read_csv(file_path, sep='\t', comment='#', header=None,
                        names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'SV_TYPE', 'SVLEN'])
        svtype_counts = df['SV_TYPE'].value_counts()
        svtype_count = svtype_counts.to_dict()
        print(f"NyuWa所有变异文件读取完成，共 {len(df)} 个变异")
        for sv_type, count in sorted(svtype_count.items(), key=lambda x: x[1], reverse=True):
            print(f"  {sv_type}: {count}")
    except Exception as e:
        print(f"错误：读取NyuWa所有变异文件失败 - {e}")
    return svtype_count

def read_svcmp_file(file_path: str) -> Tuple[Dict[str, int], Dict[str, int]]:
    """读取svcmp文件，统计千人和NyuWa在子集中的变异类型分布"""
    print("正在读取svcmp文件...")
    kg_svtype_count = {}
    nw_svtype_count = {}
    try:
        df = pd.read_csv(file_path, sep='\t', comment='#', header=None)
        if len(df.columns) < 19:
            print(f"错误：svcmp文件列数不足，只有 {len(df.columns)} 列")
            return kg_svtype_count, nw_svtype_count

        kg_df = df[[2, 17]]
        kg_df.columns = ['ID', 'SV_TYPE']
        kg_df['SV_TYPE'] = kg_df['SV_TYPE'].replace({'L1': 'LINE'})
        kg_unique = kg_df.drop_duplicates(subset=['ID'])
        kg_svtype_counts = kg_unique['SV_TYPE'].value_counts()
        kg_svtype_count = kg_svtype_counts.to_dict()

        nw_df = df[[8, 9, 18]]
        nw_df.columns = ['CHROM', 'POS', 'SV_TYPE']
        nw_unique = nw_df.drop_duplicates(subset=['CHROM', 'POS', 'SV_TYPE'])
        nw_svtype_counts = nw_unique['SV_TYPE'].value_counts()
        nw_svtype_count = nw_svtype_counts.to_dict()

        print(f"svcmp文件读取完成，共 {len(df)} 行")
        print("svcmp中千人变异类型分布:")
        for sv_type, count in sorted(kg_svtype_count.items(), key=lambda x: x[1], reverse=True):
            print(f"  {sv_type}: {count}")
        print("svcmp中NyuWa变异类型分布:")
        for sv_type, count in sorted(nw_svtype_count.items(), key=lambda x: x[1], reverse=True):
            print(f"  {sv_type}: {count}")
    except Exception as e:
        print(f"错误：读取svcmp文件失败 - {e}")
    return kg_svtype_count, nw_svtype_count

def calculate_ratios(kg_all: Dict[str, int], kg_subset: Dict[str, int],
                    nw_all: Dict[str, int], nw_subset: Dict[str, int]) -> pd.DataFrame:
    """计算比例并生成统计表格"""
    print("正在计算比例...")
    target_types = ['DEL', 'INS', 'DUP', 'INV', 'Alu', 'LINE', 'LTR', 'SVA']
    results = []
    for sv_type in target_types:
        total = kg_all.get(sv_type, 0)
        subset = kg_subset.get(sv_type, 0)
        ratio = subset / total if total > 0 else 0.0
        results.append({
            'SV_TYPE': sv_type,
            'DATASET': '1KG',
            'TOTAL_COUNT': total,
            'SUBSET_COUNT': subset,
            'RATIO': ratio
        })
    for sv_type in target_types:
        total = nw_all.get(sv_type, 0)
        subset = nw_subset.get(sv_type, 0)
        ratio = subset / total if total > 0 else 0.0
        results.append({
            'SV_TYPE': sv_type,
            'DATASET': 'NyuWa',
            'TOTAL_COUNT': total,
            'SUBSET_COUNT': subset,
            'RATIO': ratio
        })
    df = pd.DataFrame(results)
    return df

def save_statistics(df: pd.DataFrame, output_file: str):
    """保存统计结果到文件"""
    print("正在保存统计结果...")
    try:
        with open(output_file, 'w') as f:
            f.write("SV_TYPE\tDATASET\tTOTAL_COUNT\tSUBSET_COUNT\tRATIO\n")
            for _, row in df.iterrows():
                f.write(f"{row['SV_TYPE']}\t{row['DATASET']}\t{row['TOTAL_COUNT']}\t{row['SUBSET_COUNT']}\t{row['RATIO']:.6f}\n")
        print(f"统计结果已保存到: {output_file}")
    except Exception as e:
        print(f"错误：保存统计结果失败 - {e}")

def create_plot(df: pd.DataFrame, output_file: str):
    """
    创建三个子图的对比图（无数字标签，无标题，无纵轴标签，共享纵轴范围）
    - 上：NyuWa 堆叠柱状图
    - 中：1KG 堆叠柱状图
    - 下：一致性柱状图（max ratio），浅绿色
    """
    print("正在创建对比图...")
    try:
        # 准备数据
        kg_df = df[df['DATASET'] == '1KG'].sort_values('SV_TYPE')
        nw_df = df[df['DATASET'] == 'NyuWa'].sort_values('SV_TYPE')
        sv_types = ['DEL', 'INS', 'DUP', 'INV', 'Alu', 'LINE', 'LTR', 'SVA']
        kg_df = kg_df.set_index('SV_TYPE').reindex(sv_types).reset_index()
        nw_df = nw_df.set_index('SV_TYPE').reindex(sv_types).reset_index()

        # 计算一致性：每个SV类型取NyuWa和1KG比例的最大值
        consistency = [max(nw_df.loc[i, 'RATIO'], kg_df.loc[i, 'RATIO']) for i in range(len(sv_types))]

        # 确定统一的纵轴范围（基于两个数据集的总数最大值）
        max_total = max(nw_df['TOTAL_COUNT'].max(), kg_df['TOTAL_COUNT'].max())
        y_max = max_total * 1.1  # 留10%边距

        # 创建图形：3行1列，共享x轴，整体更矮（高度8）
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
        x = np.arange(len(sv_types))
        width = 0.6  # 柱子宽度

        # ----- 上子图：NyuWa 堆叠柱状图 -----
        total_nw = nw_df['TOTAL_COUNT'].values
        subset_nw = nw_df['SUBSET_COUNT'].values
        remainder_nw = total_nw - subset_nw

        ax1.bar(x, subset_nw, width, color='darkblue', label='Subset Count')
        ax1.bar(x, remainder_nw, width, bottom=subset_nw, color='lightblue', label='Remainder Count')
        ax1.set_ylim(0, y_max)
        ax1.set_ylabel('')  # 不显示纵轴标签
        ax1.legend(loc='upper right')
        ax1.grid(True, alpha=0.3, axis='y')
        ax1.tick_params(axis='x', labelbottom=False)  # 隐藏x轴刻度标签

        # ----- 中子图：1KG 堆叠柱状图 -----
        total_kg = kg_df['TOTAL_COUNT'].values
        subset_kg = kg_df['SUBSET_COUNT'].values
        remainder_kg = total_kg - subset_kg

        ax2.bar(x, subset_kg, width, color='darkred', label='Subset Count')
        ax2.bar(x, remainder_kg, width, bottom=subset_kg, color='lightcoral', label='Remainder Count')
        ax2.set_ylim(0, y_max)
        ax2.set_ylabel('')  # 不显示纵轴标签
        ax2.legend(loc='upper right')
        ax2.grid(True, alpha=0.3, axis='y')
        ax2.tick_params(axis='x', labelbottom=False)  # 隐藏x轴刻度标签

        # ----- 下子图：一致性柱状图（浅绿色）-----
        ax3.bar(x, consistency, width, color='lightgreen', label='Consistency (max ratio)')
        ax3.set_ylim(0, 1)
        ax3.set_ylabel('')  # 不显示纵轴标签
        ax3.set_xlabel('SV Type', fontsize=12, fontweight='bold')  # 明确横坐标标签
        ax3.legend(loc='upper right')
        ax3.grid(True, alpha=0.3, axis='y')

        # 设置x轴刻度标签（显示SV类型名称）
        ax3.set_xticks(x)
        ax3.set_xticklabels(sv_types, rotation=45, ha='right', fontsize=10)

        # 调整子图之间的间距
        plt.tight_layout()
        plt.savefig(output_file, bbox_inches='tight', facecolor='white', format='pdf')
        print(f"对比图已保存到: {output_file}")

    except Exception as e:
        print(f"错误：创建对比图失败 - {e}")

def main():
    # 文件路径配置（保持不变）
    kg_all_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/1019.svall.inv.addsvtype.rmsvtr.notduptoins_svlenge50.rm05mesv.rm05trsv_newsvtype_rmXY.txt"
    nw_all_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/nyuwa.svall.addsvtype.svlenge50.rm05mesv.rm05trsv_newsvtype_rmXY.txt"
    svcmp_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/svcmp.1k.sv.inv.over05.rmnan.addsvtype.rmdismatchINV.notduptoins_svlenge50.rm05mesv.rm05trsv.newsvtype.txt"
    output_stat_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/1019SV_NWSV_NWME/1019sv_1019inv_nwsv_nwme_allsvtype_detail_over05_svlenge50.rm05mesv.rm05trsv_splitdup.txt"
    output_plot_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/1019SV_NWSV_NWME/1019sv_1019inv_nwsv_nwme_allsvtype_distribution_over05_svlenge50.rm05mesv.rm05trsv_splitdup.pdf"

    print("开始统计SV变异类型分布...")
    set_custom_style()

    # 1. 读取所有变异文件
    kg_all = read_1kg_all_variants(kg_all_file)
    nw_all = read_nyuwa_all_variants(nw_all_file)

    # 2. 读取svcmp文件（子集）
    kg_subset, nw_subset = read_svcmp_file(svcmp_file)

    # 3. 计算比例
    stats_df = calculate_ratios(kg_all, kg_subset, nw_all, nw_subset)

    # 4. 保存统计结果
    save_statistics(stats_df, output_stat_file)

    # 5. 创建对比图
    create_plot(stats_df, output_plot_file)

    print("处理完成！")

if __name__ == "__main__":
    main()
