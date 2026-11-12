#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, Tuple, Set, List, Optional
import os

# ==================== AI兼容性设置 ====================
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'SimHei']


def set_custom_style():
    """设置自定义绘图样式"""
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


def normalize_sv_type(sv_type: str) -> str:
    """统一变异类型命名"""
    if pd.isna(sv_type):
        return sv_type
    sv_type = str(sv_type).strip()
    type_mapping = {
        'ALU': 'Alu',
        'alu': 'Alu',
        'L1': 'LINE',
        'l1': 'LINE'
    }
    return type_mapping.get(sv_type, sv_type)


def read_all_variants(file_path: str, dataset: str) -> Dict[str, int]:
    """读取所有变异文件，统计变异类型分布"""
    print(f"正在读取{dataset}所有变异文件...")
    svtype_count = {}
    try:
        df = pd.read_csv(file_path, sep='\t', comment='#', header=None,
                         names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'SV_TYPE', 'SVLEN'])
        df['SV_TYPE'] = df['SV_TYPE'].apply(normalize_sv_type)
        svtype_counts = df['SV_TYPE'].value_counts()
        svtype_count = svtype_counts.to_dict()
        print(f"{dataset}所有变异文件读取完成，共 {len(df)} 个变异")
        for sv_type, count in sorted(svtype_count.items(), key=lambda x: x[1], reverse=True):
            print(f"  {sv_type}: {count}")
    except Exception as e:
        print(f"错误：读取{dataset}所有变异文件失败 - {e}")
    return svtype_count


def read_cmp_1019_nw(file_path: str) -> Tuple[Dict[str, Set[Tuple]], Dict[str, Set[str]]]:
    """读取1019和NyuWa比较文件，返回按类型分组的变异集合"""
    print("正在读取1019-NyuWa比较文件...")
    nw_variants_by_type = {}
    kg_variants_by_type = {}
    target_types = ['DEL', 'INS', 'DUP', 'INV', 'Alu', 'LINE', 'LTR', 'SVA']
    for sv_type in target_types:
        nw_variants_by_type[sv_type] = set()
        kg_variants_by_type[sv_type] = set()
    try:
        df = pd.read_csv(file_path, sep='\t', comment='#', header=None)
        if len(df.columns) < 19:
            print(f"错误：1019-NyuWa比较文件列数不足，只有 {len(df.columns)} 列")
            return nw_variants_by_type, kg_variants_by_type
        nw_chrom = df[8]
        nw_pos = df[9]
        nw_svtype = df[18]
        nw_svlen = df[13]
        kg_id = df[2]
        kg_svtype = df[17]
        nw_svtype = nw_svtype.apply(normalize_sv_type)
        kg_svtype = kg_svtype.apply(normalize_sv_type)
        for i in range(len(df)):
            nw_type = nw_svtype.iloc[i] if i < len(nw_svtype) else None
            kg_type = kg_svtype.iloc[i] if i < len(kg_svtype) else None
            if pd.notna(nw_type) and nw_type in nw_variants_by_type:
                variant_id = (nw_chrom.iloc[i], nw_pos.iloc[i], nw_type, nw_svlen.iloc[i])
                nw_variants_by_type[nw_type].add(variant_id)
            if pd.notna(kg_type) and kg_type in kg_variants_by_type and pd.notna(kg_id.iloc[i]):
                kg_variants_by_type[kg_type].add(str(kg_id.iloc[i]))
        print(f"1019-NyuWa比较文件读取完成，共 {len(df)} 行")
        for sv_type in target_types:
            print(f"  {sv_type}: NyuWa={len(nw_variants_by_type[sv_type])}, 1019={len(kg_variants_by_type[sv_type])}")
    except Exception as e:
        print(f"错误：读取1019-NyuWa比较文件失败 - {e}")
    return nw_variants_by_type, kg_variants_by_type


def read_cmp_65_nw(file_path: str) -> Dict[str, Set[Tuple]]:
    """读取65和NyuWa比较文件，返回按类型分组的NyuWa变异集合"""
    print("正在读取65-NyuWa比较文件...")
    nw_variants_by_type = {}
    target_types = ['DEL', 'INS', 'DUP', 'INV', 'Alu', 'LINE', 'LTR', 'SVA']
    for sv_type in target_types:
        nw_variants_by_type[sv_type] = set()
    try:
        df = pd.read_csv(file_path, sep='\t', comment='#', header=None)
        if len(df.columns) < 19:
            print(f"错误：65-NyuWa比较文件列数不足，只有 {len(df.columns)} 列")
            return nw_variants_by_type
        nw_chrom = df[8]
        nw_pos = df[9]
        nw_svtype = df[18]
        nw_svlen = df[13]
        nw_svtype = nw_svtype.apply(normalize_sv_type)
        for i in range(len(df)):
            sv_type = nw_svtype.iloc[i] if i < len(nw_svtype) else None
            if pd.notna(sv_type) and sv_type in nw_variants_by_type:
                variant_id = (nw_chrom.iloc[i], nw_pos.iloc[i], sv_type, nw_svlen.iloc[i])
                nw_variants_by_type[sv_type].add(variant_id)
        print(f"65-NyuWa比较文件读取完成，共 {len(df)} 行")
        for sv_type in target_types:
            print(f"  {sv_type}: NyuWa={len(nw_variants_by_type[sv_type])}")
    except Exception as e:
        print(f"错误：读取65-NyuWa比较文件失败 - {e}")
    return nw_variants_by_type


def read_cmp_65_1019(file_path: str) -> Dict[str, Set[str]]:
    """读取65和1019比较文件，返回按类型分组的1019变异集合"""
    print("正在读取65-1019比较文件...")
    kg_variants_by_type = {}
    target_types = ['DEL', 'INS', 'DUP', 'INV', 'Alu', 'LINE', 'LTR', 'SVA']
    for sv_type in target_types:
        kg_variants_by_type[sv_type] = set()
    try:
        df = pd.read_csv(file_path, sep='\t', comment='#', header=None)
        if len(df.columns) < 17:
            print(f"错误：65-1019比较文件列数不足，只有 {len(df.columns)} 列")
            return kg_variants_by_type
        kg_id = df[2]
        kg_svtype = df[15]
        kg_svtype = kg_svtype.apply(normalize_sv_type)
        for i in range(len(df)):
            sv_type = kg_svtype.iloc[i] if i < len(kg_svtype) else None
            if pd.notna(sv_type) and sv_type in kg_variants_by_type and pd.notna(kg_id.iloc[i]):
                kg_variants_by_type[sv_type].add(str(kg_id.iloc[i]))
        print(f"65-1019比较文件读取完成，共 {len(df)} 行")
        for sv_type in target_types:
            print(f"  {sv_type}: 1019={len(kg_variants_by_type[sv_type])}")
    except Exception as e:
        print(f"错误：读取65-1019比较文件失败 - {e}")
    return kg_variants_by_type


def analyze_nyuwa_overlap(nw_all: Dict[str, int],
                          nw_1019_variants: Dict[str, Set[Tuple]],
                          nw_65_variants: Dict[str, Set[Tuple]]) -> pd.DataFrame:
    """分析NyuWa变异在不同比较文件中的重叠关系"""
    print("正在分析NyuWa变异重叠关系...")
    target_types = ['DEL', 'INS', 'DUP', 'INV', 'Alu', 'LINE', 'LTR', 'SVA']
    results = []
    for sv_type in target_types:
        total_count = nw_all.get(sv_type, 0)
        nw_1019_type = len(nw_1019_variants.get(sv_type, set()))
        nw_65_type = len(nw_65_variants.get(sv_type, set()))
        nw_1019_set = nw_1019_variants.get(sv_type, set())
        nw_65_set = nw_65_variants.get(sv_type, set())
        nw_both_type = len(nw_1019_set.intersection(nw_65_set))
        nw_only_1019 = nw_1019_type - nw_both_type
        nw_only_65 = nw_65_type - nw_both_type
        nw_unique = total_count - nw_1019_type - nw_only_65
        nw_unique = max(0, nw_unique)
        nw_only_1019 = max(0, nw_only_1019)
        nw_only_65 = max(0, nw_only_65)
        nw_both_type = max(0, nw_both_type)
        parts_sum = nw_unique + nw_both_type + nw_only_1019 + nw_only_65
        if parts_sum != total_count and total_count > 0:
            nw_unique = total_count - (nw_both_type + nw_only_1019 + nw_only_65)
            nw_unique = max(0, nw_unique)
        results.append({
            'SV_TYPE': sv_type,
            'DATASET': 'NyuWa',
            'TOTAL_COUNT': total_count,
            'UNIQUE': nw_unique,
            'BOTH': nw_both_type,
            'ONLY_1019': nw_only_1019,
            'ONLY_65': nw_only_65
        })
        print(f"  {sv_type}: 总数={total_count}, 唯一={nw_unique}, 两者都有={nw_both_type}, 仅1019={nw_only_1019}, 仅65={nw_only_65}")
    return pd.DataFrame(results)


def analyze_1019_overlap(kg_all: Dict[str, int],
                         kg_1019_nw_variants: Dict[str, Set[str]],
                         kg_65_1019_variants: Dict[str, Set[str]]) -> pd.DataFrame:
    """分析1019变异在不同比较文件中的重叠关系"""
    print("正在分析1019变异重叠关系...")
    target_types = ['DEL', 'INS', 'DUP', 'INV', 'Alu', 'LINE', 'LTR', 'SVA']
    results = []
    for sv_type in target_types:
        total_count = kg_all.get(sv_type, 0)
        kg_1019_nw_type = len(kg_1019_nw_variants.get(sv_type, set()))
        kg_65_1019_type = len(kg_65_1019_variants.get(sv_type, set()))
        kg_1019_nw_set = kg_1019_nw_variants.get(sv_type, set())
        kg_65_1019_set = kg_65_1019_variants.get(sv_type, set())
        kg_both_type = len(kg_1019_nw_set.intersection(kg_65_1019_set))
        kg_only_nyuwa = kg_1019_nw_type - kg_both_type
        kg_only_65 = kg_65_1019_type - kg_both_type
        kg_unique = total_count - kg_1019_nw_type - kg_only_65
        kg_unique = max(0, kg_unique)
        kg_only_nyuwa = max(0, kg_only_nyuwa)
        kg_only_65 = max(0, kg_only_65)
        kg_both_type = max(0, kg_both_type)
        parts_sum = kg_unique + kg_both_type + kg_only_nyuwa + kg_only_65
        if parts_sum != total_count and total_count > 0:
            kg_unique = total_count - (kg_both_type + kg_only_nyuwa + kg_only_65)
            kg_unique = max(0, kg_unique)
        results.append({
            'SV_TYPE': sv_type,
            'DATASET': '1019',
            'TOTAL_COUNT': total_count,
            'UNIQUE': kg_unique,
            'BOTH': kg_both_type,
            'ONLY_NYW': kg_only_nyuwa,
            'ONLY_65': kg_only_65
        })
        print(f"  {sv_type}: 总数={total_count}, 唯一={kg_unique}, 两者都有={kg_both_type}, 仅NyuWa={kg_only_nyuwa}, 仅65={kg_only_65}")
    return pd.DataFrame(results)


def save_detailed_statistics_with_precision(nyuwa_df: pd.DataFrame, kg_df: pd.DataFrame,
                                            output_file: str):
    """
    保存详细的统计结果，包含Precision列
    Precision = 1 - UNIQUE / TOTAL_COUNT  (若 TOTAL_COUNT=0 则为0)
    """
    print("正在保存详细统计结果（含Precision）...")
    try:
        with open(output_file, 'w') as f:
            f.write("SV_TYPE\tDATASET\tTOTAL_COUNT\tUNIQUE\tSHARE_BOTH\tSHARE_1019\tSHARE_65\tSHARE_NYW\tPRECISION\n")
            # NyuWa数据
            for _, row in nyuwa_df.iterrows():
                sv_type = row['SV_TYPE']
                total = row['TOTAL_COUNT']
                unique = row['UNIQUE']
                precision = 1 - (unique / total) if total > 0 else 0.0
                f.write(f"{sv_type}\t{row['DATASET']}\t{total}\t{unique}\t{row['BOTH']}\t{row['ONLY_1019']}\t{row['ONLY_65']}\t0\t{precision:.4f}\n")
            # 1019数据
            for _, row in kg_df.iterrows():
                sv_type = row['SV_TYPE']
                total = row['TOTAL_COUNT']
                unique = row['UNIQUE']
                precision = 1 - (unique / total) if total > 0 else 0.0
                f.write(f"{sv_type}\t{row['DATASET']}\t{total}\t{unique}\t{row['BOTH']}\t0\t{row['ONLY_65']}\t{row['ONLY_NYW']}\t{precision:.4f}\n")
        print(f"详细统计结果（含Precision）已保存到: {output_file}")
    except Exception as e:
        print(f"错误：保存详细统计结果失败 - {e}")


def create_precision_line_plot(nyuwa_df: pd.DataFrame, kg_df: pd.DataFrame,
                               output_file: str):
    """
    创建Precision曲线图（单图双线）
    NyuWa: 红色点线，1019: 蓝色点线
    """
    print("正在创建Precision曲线图...")
    try:
        # 设置SV类型顺序
        sv_types = ['DEL', 'INS', 'DUP', 'INV', 'Alu', 'LINE', 'LTR', 'SVA']
        nyuwa_df = nyuwa_df.set_index('SV_TYPE').reindex(sv_types).reset_index()
        kg_df = kg_df.set_index('SV_TYPE').reindex(sv_types).reset_index()

        # 计算Precision
        nw_precision = []
        for _, row in nyuwa_df.iterrows():
            total = row['TOTAL_COUNT']
            unique = row['UNIQUE']
            nw_precision.append(1 - unique / total if total > 0 else 0.0)

        kg_precision = []
        for _, row in kg_df.iterrows():
            total = row['TOTAL_COUNT']
            unique = row['UNIQUE']
            kg_precision.append(1 - unique / total if total > 0 else 0.0)

        # 绘图
        fig, ax = plt.subplots(figsize=(10, 7))
        x = np.arange(len(sv_types))

        # 绘制两条点线图
        ax.plot(x, nw_precision, 'r-o', linewidth=2, markersize=8, label='NyuWa')
        ax.plot(x, kg_precision, 'b-o', linewidth=2, markersize=8, label='1019')

        # 设置x轴标签
        ax.set_xticks(x)
        ax.set_xticklabels(sv_types, rotation=45, ha='right', fontsize=12)
        ax.set_xlabel('SV Type', fontsize=14, fontweight='bold')
        ax.set_ylabel('Precision', fontsize=14, fontweight='bold')
        ax.set_ylim(0, 1.05)
        ax.set_title('Precision of SV Calls across Different Datasets', fontsize=16, fontweight='bold')
        ax.legend(loc='best', fontsize=12)
        ax.grid(True, alpha=0.3, linestyle='--')

        plt.tight_layout()
        plt.savefig(output_file, bbox_inches='tight', facecolor='white', format='pdf')
        print(f"Precision曲线图已保存到: {output_file}")
    except Exception as e:
        print(f"错误：创建Precision曲线图失败 - {e}")


def main():
    # 文件路径配置
    kg_all_file = "/Parastor300s_G30S/haoxp/third/cmp/65/1019.svall.inv.addsvtype.rmsvtr.notduptoins_svlenge50.rm05mesv.rm05trsv_newsvtype.txt"
    nw_all_file = "/Parastor300s_G30S/haoxp/third/cmp/65/nyuwa.svall.addsvtype.svlenge50.rm05mesv.rm05trsv_newsvtype.txt"
    cmp_1019_nw_file = "/Parastor300s_G30S/haoxp/third/cmp/65/svcmp.1019.nw.sv.inv.over05.rmnan.addsvtype.rmdismatchINV.notduptoins_svlenge50.rm05mesv.rm05trsv.newsvtype.txt"
    cmp_65_nw_file = "/Parastor300s_G30S/haoxp/third/cmp/65/svcmp.65.nw.sv.inv.over05.rmnan.addsvtype.notduptoins_svlenge50.rm05mesv.rm05trsv.newsvtype.txt"
    cmp_65_1019_file = "/Parastor300s_G30S/haoxp/third/cmp/65/svcmp.65.1019.sv.inv.over05.rmnan.addsvtype.rmdismatchINV.notduptoins_svlenge50.rm05mesv.rm05trsv.newsvtype.txt"

    output_stat_file = "/Parastor300s_G30S/haoxp/third/cmp/65/65allsv_1019sv_1019inv_nwsv_nwme_allsvtype_detail_over05.notduptoins_svlenge50.rm05mesv.rm05trsv.precision.line.txt"
    output_plot_file = "/Parastor300s_G30S/haoxp/third/cmp/65/65allsv_1019sv_1019inv_nwsv_nwme_allsvtype_distribution_over05.notduptoins_svlenge50.rm05mesv.rm05trsv.precision.line.pdf"

    print("开始分析多数据集SV变异类型Precision曲线...")

    set_custom_style()

    # 1. 读取所有变异文件
    kg_all = read_all_variants(kg_all_file, "1019")
    nw_all = read_all_variants(nw_all_file, "NyuWa")

    # 2. 读取比较文件
    nw_1019_variants, kg_1019_nw_variants = read_cmp_1019_nw(cmp_1019_nw_file)
    nw_65_variants = read_cmp_65_nw(cmp_65_nw_file)
    kg_65_1019_variants = read_cmp_65_1019(cmp_65_1019_file)

    # 3. 分析重叠关系
    nyuwa_df = analyze_nyuwa_overlap(nw_all, nw_1019_variants, nw_65_variants)
    kg_df = analyze_1019_overlap(kg_all, kg_1019_nw_variants, kg_65_1019_variants)

    # 4. 保存详细统计结果（含Precision）
    save_detailed_statistics_with_precision(nyuwa_df, kg_df, output_stat_file)

    # 5. 创建Precision曲线图
    create_precision_line_plot(nyuwa_df, kg_df, output_plot_file)

    print("处理完成！")


if __name__ == "__main__":
    main()
