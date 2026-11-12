#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, Tuple, Set

def read_1kg_all_variants(file_path: str) -> Dict[str, int]:
    """
    读取千人所有变异文件，统计变异类型分布
    """
    print("正在读取千人所有变异文件...")
    svtype_count = {}

    try:
        # 使用pandas读取，跳过注释行
        df = pd.read_csv(file_path, sep='\t', comment='#', header=None,
                        names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'SV_TYPE', 'SVLEN'])

        # 统计变异类型
        svtype_counts = df['SV_TYPE'].value_counts()
        svtype_count = svtype_counts.to_dict()

        print(f"千人所有变异文件读取完成，共 {len(df)} 个变异")
        print("千人所有变异类型分布:")
        for sv_type, count in sorted(svtype_count.items(), key=lambda x: x[1], reverse=True):
            print(f"  {sv_type}: {count}")

    except Exception as e:
        print(f"错误：读取千人所有变异文件失败 - {e}")

    return svtype_count

def read_nyuwa_all_variants(file_path: str) -> Dict[str, int]:
    """
    读取NyuWa所有变异文件，统计变异类型分布
    """
    print("正在读取NyuWa所有变异文件...")
    svtype_count = {}

    try:
        # 使用pandas读取，跳过注释行
        df = pd.read_csv(file_path, sep='\t', comment='#', header=None,
                        names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'SV_TYPE', 'SVLEN'])

        # 统计变异类型
        svtype_counts = df['SV_TYPE'].value_counts()
        svtype_count = svtype_counts.to_dict()

        print(f"NyuWa所有变异文件读取完成，共 {len(df)} 个变异")
        print("NyuWa所有变异类型分布:")
        for sv_type, count in sorted(svtype_count.items(), key=lambda x: x[1], reverse=True):
            print(f"  {sv_type}: {count}")

    except Exception as e:
        print(f"错误：读取NyuWa所有变异文件失败 - {e}")

    return svtype_count

def read_svcmp_file(file_path: str) -> Tuple[Dict[str, int], Dict[str, int]]:
    """
    读取svcmp文件，统计千人和NyuWa在子集中的变异类型分布
    """
    print("正在读取svcmp文件...")
    kg_svtype_count = {}
    nw_svtype_count = {}

    try:
        # 读取文件
        df = pd.read_csv(file_path, sep='\t', comment='#', header=None)

        # 检查列数
        if len(df.columns) < 19:
            print(f"错误：svcmp文件列数不足，只有 {len(df.columns)} 列")
            return kg_svtype_count, nw_svtype_count

        # 千人变异：根据变异ID去重统计
        kg_df = df[[2, 17]]  # 第三列是变异ID，第十八列是千人变异类型
        kg_df.columns = ['ID', 'SV_TYPE']

        # 映射L1到LINE
        kg_df['SV_TYPE'] = kg_df['SV_TYPE'].replace({'L1': 'LINE'})

        kg_unique = kg_df.drop_duplicates(subset=['ID'])  # 根据变异ID去重
        kg_svtype_counts = kg_unique['SV_TYPE'].value_counts()
        kg_svtype_count = kg_svtype_counts.to_dict()

        # NyuWa变异：根据染色体、位置、变异类型去重统计
        nw_df = df[[8, 9, 18]]  # 第九列是染色体，第十列是位置，第十九列是NyuWa变异类型
        nw_df.columns = ['CHROM', 'POS', 'SV_TYPE']
        nw_unique = nw_df.drop_duplicates(subset=['CHROM', 'POS', 'SV_TYPE'])  # 根据三个条件去重
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
    """
    计算比例并生成统计表格
    """
    print("正在计算比例...")

    # 定义需要统计的7种变异类型（去掉DUP）
    target_types = ['DEL', 'INS', 'INV', 'Alu', 'LINE', 'LTR', 'SVA']

    results = []

    # 处理千人数据
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

    # 处理NyuWa数据
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
    """
    保存统计结果到文件
    """
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
    创建对比图
    """
    print("正在创建对比图...")

    try:
        # 设置字体 - 使用系统默认字体避免警告
        plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Liberation Sans']
        plt.rcParams['axes.unicode_minus'] = False

        # 分离千人和NyuWa数据
        kg_df = df[df['DATASET'] == '1KG'].sort_values('SV_TYPE')
        nw_df = df[df['DATASET'] == 'NyuWa'].sort_values('SV_TYPE')

        # 定义变异类型顺序（去掉DUP）
        sv_types = ['DEL', 'INS', 'INV', 'Alu', 'LINE', 'LTR', 'SVA']

        # 确保数据顺序一致
        kg_df = kg_df.set_index('SV_TYPE').reindex(sv_types).reset_index()
        nw_df = nw_df.set_index('SV_TYPE').reindex(sv_types).reset_index()

        # 创建图形
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))

        # 设置x轴位置
        x = np.arange(len(sv_types))
        width = 0.8

        # 计算统一的y轴范围
        max_count = max(kg_df['TOTAL_COUNT'].max(), nw_df['TOTAL_COUNT'].max())
        y_max_count = max_count * 1.15  # 留15%的边距

        # 绘制NyuWa图（上部）
        bars1 = ax1.bar(x, nw_df['TOTAL_COUNT'], width, color='lightblue', alpha=0.7, label='Total Count')
        ax1_twin = ax1.twinx()
        line1 = ax1_twin.plot(x, nw_df['RATIO'], 'ro-', linewidth=2, markersize=8, label='Ratio')

        # 设置NyuWa图的坐标轴
        ax1.set_ylabel('NyuWa Total Count', fontsize=12, fontweight='bold')
        ax1.set_ylim(0, y_max_count)
        ax1_twin.set_ylabel('Ratio', fontsize=12, fontweight='bold')
        ax1_twin.set_ylim(0, 1)
        ax1.set_title('NyuWa SV Type Distribution', fontsize=14, fontweight='bold')
        ax1.set_xticks(x)
        ax1.set_xticklabels(sv_types, rotation=45, ha='right', fontsize=10)

        # 在柱子上添加数量标签
        for i, bar in enumerate(bars1):
            height = bar.get_height()
            if height > 0:  # 只在有数据时显示标签
                ax1.text(bar.get_x() + bar.get_width()/2., height + y_max_count*0.01,
                        f'{int(height):,}', ha='center', va='bottom', fontsize=9, fontweight='bold')

        # 在折线点上添加比例标签
        for i, ratio in enumerate(nw_df['RATIO']):
            if ratio > 0:  # 只在有数据时显示标签
                ax1_twin.text(x[i], ratio + 0.03, f'{ratio:.2f}',
                            ha='center', va='bottom', fontsize=9, fontweight='bold', color='red')

        # 绘制千人图（下部）
        bars2 = ax2.bar(x, kg_df['TOTAL_COUNT'], width, color='lightcoral', alpha=0.7, label='Total Count')
        ax2_twin = ax2.twinx()
        line2 = ax2_twin.plot(x, kg_df['RATIO'], 'ro-', linewidth=2, markersize=8, label='Ratio')

        # 设置千人图的坐标轴
        ax2.set_ylabel('1KG Total Count', fontsize=12, fontweight='bold')
        ax2.set_ylim(0, y_max_count)
        ax2_twin.set_ylabel('Ratio', fontsize=12, fontweight='bold')
        ax2_twin.set_ylim(0, 1)
        ax2.set_title('1KG SV Type Distribution', fontsize=14, fontweight='bold')
        ax2.set_xticks(x)
        ax2.set_xticklabels(sv_types, rotation=45, ha='right', fontsize=10)

        # 在柱子上添加数量标签
        for i, bar in enumerate(bars2):
            height = bar.get_height()
            if height > 0:  # 只在有数据时显示标签
                ax2.text(bar.get_x() + bar.get_width()/2., height + y_max_count*0.01,
                        f'{int(height):,}', ha='center', va='bottom', fontsize=9, fontweight='bold')

        # 在折线点上添加比例标签
        for i, ratio in enumerate(kg_df['RATIO']):
            if ratio > 0:  # 只在有数据时显示标签
                ax2_twin.text(x[i], ratio + 0.03, f'{ratio:.2f}',
                            ha='center', va='bottom', fontsize=9, fontweight='bold', color='red')

        # 添加图例
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax1_twin.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

        lines1, labels1 = ax2.get_legend_handles_labels()
        lines2, labels2 = ax2_twin.get_legend_handles_labels()
        ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

        # 添加网格线
        ax1.grid(True, alpha=0.3, axis='y')
        ax2.grid(True, alpha=0.3, axis='y')

        # 调整布局
        plt.tight_layout()

        # 保存图片
#        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.savefig(output_file, bbox_inches='tight', facecolor='white', format='pdf')
        print(f"对比图已保存到: {output_file}")

    except Exception as e:
        print(f"错误：创建对比图失败 - {e}")

def main():
    # 文件路径配置
    kg_all_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/1019.svall.inv.addsvtype.rmsvtr_svlenge50.rm05mesv.rm05trsv_newsvtype_rmXY.txt"
    nw_all_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/nyuwa.svall.addsvtype.duptoins_svlenge50.rm05mesv.rm05trsv_newsvtype_rmXY.txt"
    svcmp_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/svcmp.1k.sv.inv.over05.rmnan.addsvtype.rmdismatchINV_svlenge50.rm05mesv.rm05trsv.newsvtype.txt"
    output_stat_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/1019SV_NWSV_NWME/1019sv_1019inv_nwsv_nwme_allsvtype_detail_over05_svlenge50.rm05mesv.rm05trsv.txt"
    output_plot_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/1019SV_NWSV_NWME/1019sv_1019inv_nwsv_nwme_allsvtype_distribution_over05_svlenge50.rm05mesv.rm05trsv.pdf"

    print("开始统计SV变异类型分布...")

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
