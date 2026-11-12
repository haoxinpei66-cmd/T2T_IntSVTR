#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from typing import Dict, Tuple, Set, List, Optional
import sys


# ==================== AI兼容性设置 ====================
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'SimHei']

plt.rcParams.update({
    'font.size': 24,
    'axes.titlesize': 20,
    'axes.labelsize': 28,
    'xtick.labelsize': 26,
    'ytick.labelsize': 26,
    'legend.fontsize': 24,
    'figure.titlesize': 22,
    'figure.dpi': 300,
})


def normalize_sv_type(sv_type: str) -> str:
    if pd.isna(sv_type):
        return sv_type
    sv_type = str(sv_type).strip()
    type_mapping = {'ALU': 'Alu', 'alu': 'Alu', 'L1': 'LINE', 'l1': 'LINE'}
    return type_mapping.get(sv_type, sv_type)


def read_all_variants(file_path: str, dataset: str) -> pd.DataFrame:
    print(f"正在读取{dataset}所有变异文件...")
    try:
        df = pd.read_csv(file_path, sep='\t', comment='#', header=None,
                        names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'SV_TYPE', 'SVLEN'])
        df['SV_TYPE'] = df['SV_TYPE'].apply(normalize_sv_type)
        df['SVLEN'] = pd.to_numeric(df['SVLEN'], errors='coerce')
        df = df[df['SVLEN'].notna()]
        df['SVLEN'] = df['SVLEN'].abs()
        print(f"{dataset}所有变异文件读取完成，共 {len(df)} 个变异")
        print(f"{dataset}变异类型分布: {df['SV_TYPE'].value_counts().to_dict()}")
        return df
    except Exception as e:
        print(f"错误：读取{dataset}所有变异文件失败 - {e}")
        return pd.DataFrame()


def read_cmp_1019_nw(file_path: str) -> Tuple[Dict[str, Set[Tuple]], Dict[str, Set[str]]]:
    print("正在读取1019-NyuWa比较文件...")
    nw_variants_by_type = {}
    kg_variants_by_type = {}
    target_types = ['Alu', 'DEL', 'INS', 'INV', 'LINE', 'LTR', 'SVA', 'DUP']
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
        nw_svlen = pd.to_numeric(nw_svlen, errors='coerce')
        for i in range(len(df)):
            nw_type = nw_svtype.iloc[i] if i < len(nw_svtype) else None
            kg_type = kg_svtype.iloc[i] if i < len(kg_svtype) else None
            if pd.notna(nw_type) and nw_type in nw_variants_by_type and pd.notna(nw_svlen.iloc[i]):
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
    print("正在读取65-NyuWa比较文件...")
    nw_variants_by_type = {}
    target_types = ['Alu', 'DEL', 'INS', 'INV', 'LINE', 'LTR', 'SVA', 'DUP']
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
        nw_svlen = pd.to_numeric(nw_svlen, errors='coerce')
        for i in range(len(df)):
            sv_type = nw_svtype.iloc[i] if i < len(nw_svtype) else None
            if pd.notna(sv_type) and sv_type in nw_variants_by_type and pd.notna(nw_svlen.iloc[i]):
                variant_id = (nw_chrom.iloc[i], nw_pos.iloc[i], sv_type, nw_svlen.iloc[i])
                nw_variants_by_type[sv_type].add(variant_id)
        print(f"65-NyuWa比较文件读取完成，共 {len(df)} 行")
        for sv_type in target_types:
            print(f"  {sv_type}: NyuWa={len(nw_variants_by_type[sv_type])}")
    except Exception as e:
        print(f"错误：读取65-NyuWa比较文件失败 - {e}")
    return nw_variants_by_type


def read_cmp_65_1019(file_path: str) -> Dict[str, Set[str]]:
    print("正在读取65-1019比较文件...")
    kg_variants_by_type = {}
    target_types = ['Alu', 'DEL', 'INS', 'INV', 'LINE', 'LTR', 'SVA', 'DUP']
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


def create_svlen_bins(sv_type: str, plot_type: str = "default") -> List[Tuple[str, Tuple]]:
    bins = []
    if sv_type == 'Alu':
        if plot_type == "detailed":
            for i in range(50, 100, 10):
                start = i; end = i + 9
                bins.append((f"{start}-{end}", (start, end)))
            for i in range(100, 400, 50):
                start = i + 1; end = i + 50
                bins.append((f"{start}-{end}", (start, end)))
            for i in range(400, 1000, 100):
                start = i + 1; end = i + 100
                bins.append((f"{start}-{end}", (start, end)))
            bins.append((">1000", (1001, float('inf'))))
        else:
            bins.append(("50-70", (50, 70)))
            bins.append(("71-100", (71, 100)))
            bins.append(("101-200", (101, 200)))
            bins.append(("201-400", (201, 400)))
            bins.append(("401-800", (401, 800)))
            bins.append((">800", (801, float('inf'))))
    elif sv_type == 'DEL':
        bins.append(("50-60", (50, 60))); bins.append(("61-70", (61, 70))); bins.append(("71-80", (71, 80)))
        bins.append(("81-100", (81, 100))); bins.append(("101-200", (101, 200))); bins.append(("201-300", (201, 300)))
        bins.append(("301-400", (301, 400))); bins.append(("401-600", (401, 600))); bins.append(("601-800", (601, 800)))
        bins.append(("801-1000", (801, 1000))); bins.append(("1001-2000", (1001, 2000))); bins.append(("2001-4000", (2001, 4000)))
        bins.append(("4001-10000", (4001, 10000))); bins.append(("10001-20000", (10001, 20000))); bins.append(("20001-40000", (20001, 40000)))
        bins.append(("40001-100000", (40001, 100000))); bins.append((">100000", (100001, float('inf'))))
    elif sv_type in ['INS', 'DUP', 'INSMERGE_DUP']:
        bins.append(("50-60", (50, 60))); bins.append(("61-70", (61, 70))); bins.append(("71-80", (71, 80)))
        bins.append(("81-100", (81, 100))); bins.append(("101-200", (101, 200))); bins.append(("201-300", (201, 300)))
        bins.append(("301-400", (301, 400))); bins.append(("401-600", (401, 600))); bins.append(("601-1000", (601, 1000)))
        bins.append(("1001-2000", (1001, 2000))); bins.append(("2001-3000", (2001, 3000))); bins.append(("3001-4000", (3001, 4000)))
        bins.append(("4001-6000", (4001, 6000))); bins.append(("6001-10000", (6001, 10000))); bins.append((">10000", (10001, float('inf'))))
    elif sv_type == 'LTR':
        bins.append(("50-100", (50, 100))); bins.append(("101-200", (101, 200))); bins.append(("201-300", (201, 300)))
        bins.append(("301-400", (301, 400))); bins.append(("401-600", (401, 600))); bins.append(("601-1000", (601, 1000)))
        bins.append(("1001-2000", (1001, 2000))); bins.append(("2001-4000", (2001, 4000))); bins.append(("4001-10000", (4001, 10000)))
        bins.append((">10000", (10001, float('inf'))))
    else:  # LINE, INV, SVA
        bins.append(("50-75", (50, 75))); bins.append(("76-100", (76, 100))); bins.append(("101-200", (101, 200)))
        bins.append(("201-300", (201, 300))); bins.append(("301-400", (301, 400))); bins.append(("401-600", (401, 600)))
        bins.append(("601-1000", (601, 1000))); bins.append(("1001-2000", (1001, 2000))); bins.append(("2001-3000", (2001, 3000)))
        bins.append(("3001-4000", (3001, 4000))); bins.append(("4001-6000", (4001, 6000))); bins.append(("6001-10000", (6001, 10000)))
        bins.append((">10000", (10001, float('inf'))))
    return bins


def assign_svlen_bin(svlen: float, bins: List[Tuple[str, Tuple]]) -> str:
    for bin_name, (start, end) in bins:
        if start <= svlen <= end:
            return bin_name
    return "Other"


def analyze_nyuwa_svlen_distribution(nw_all: pd.DataFrame,
                                   nw_1019_variants: Dict[str, Set[Tuple]],
                                   nw_65_variants: Dict[str, Set[Tuple]],
                                   sv_type: str, bins: List[Tuple[str, Tuple]]) -> pd.DataFrame:
    print(f"正在分析NyuWa {sv_type} SV长度分布...")
    nw_sv_type = nw_all[nw_all['SV_TYPE'] == sv_type].copy()
    if len(nw_sv_type) == 0:
        print(f"警告：NyuWa总集中没有找到 {sv_type} 类型的变异")
        return pd.DataFrame()
    nw_sv_type.loc[:, 'SVLEN_BIN'] = nw_sv_type['SVLEN'].apply(lambda x: assign_svlen_bin(x, bins))
    nw_sv_type.loc[:, 'composite_key'] = nw_sv_type['CHROM'] + '_' + nw_sv_type['POS'].astype(str) + '_' + nw_sv_type['SV_TYPE'] + '_' + nw_sv_type['SVLEN'].astype(str)
    nw_1019_keys = set()
    for variant in nw_1019_variants.get(sv_type, set()):
        chrom, pos, svtype, svlen = variant
        key = f"{chrom}_{pos}_{svtype}_{svlen}"
        nw_1019_keys.add(key)
    nw_65_keys = set()
    for variant in nw_65_variants.get(sv_type, set()):
        chrom, pos, svtype, svlen = variant
        key = f"{chrom}_{pos}_{svtype}_{svlen}"
        nw_65_keys.add(key)
    results = []
    bin_names = [bin_name for bin_name, _ in bins]
    for bin_name in bin_names:
        bin_variants = nw_sv_type[nw_sv_type['SVLEN_BIN'] == bin_name]
        total_count = len(bin_variants)
        if total_count == 0:
            results.append({'SVLEN_BIN': bin_name, 'TOTAL_COUNT': 0, 'UNIQUE': 0, 'BOTH': 0, 'ONLY_1019': 0, 'ONLY_65': 0})
            continue
        both_count = only_1019_count = only_65_count = unique_count = 0
        for _, variant in bin_variants.iterrows():
            key = variant['composite_key']
            in_1019 = key in nw_1019_keys
            in_65 = key in nw_65_keys
            if in_1019 and in_65:
                both_count += 1
            elif in_1019:
                only_1019_count += 1
            elif in_65:
                only_65_count += 1
            else:
                unique_count += 1
        calculated_total = both_count + only_1019_count + only_65_count + unique_count
        if calculated_total != total_count:
            unique_count = total_count - (both_count + only_1019_count + only_65_count)
            unique_count = max(0, unique_count)
        results.append({'SVLEN_BIN': bin_name, 'TOTAL_COUNT': total_count, 'UNIQUE': unique_count,
                        'BOTH': both_count, 'ONLY_1019': only_1019_count, 'ONLY_65': only_65_count})
    return pd.DataFrame(results)


def analyze_1019_svlen_distribution(kg_all: pd.DataFrame,
                                  kg_1019_nw_variants: Dict[str, Set[str]],
                                  kg_65_1019_variants: Dict[str, Set[str]],
                                  sv_type: str, bins: List[Tuple[str, Tuple]]) -> pd.DataFrame:
    print(f"正在分析1019 {sv_type} SV长度分布...")
    kg_sv_type = kg_all[kg_all['SV_TYPE'] == sv_type].copy()
    if len(kg_sv_type) == 0:
        print(f"警告：1019总集中没有找到 {sv_type} 类型的变异")
        return pd.DataFrame()
    kg_sv_type.loc[:, 'SVLEN_BIN'] = kg_sv_type['SVLEN'].apply(lambda x: assign_svlen_bin(x, bins))
    results = []
    bin_names = [bin_name for bin_name, _ in bins]
    for bin_name in bin_names:
        bin_variants = kg_sv_type[kg_sv_type['SVLEN_BIN'] == bin_name]
        total_count = len(bin_variants)
        if total_count == 0:
            results.append({'SVLEN_BIN': bin_name, 'TOTAL_COUNT': 0, 'UNIQUE': 0, 'BOTH': 0, 'ONLY_NYW': 0, 'ONLY_65': 0})
            continue
        both_count = only_nyw_count = only_65_count = unique_count = 0
        for _, variant in bin_variants.iterrows():
            variant_id = variant['ID']
            in_nyw = variant_id in kg_1019_nw_variants.get(sv_type, set())
            in_65 = variant_id in kg_65_1019_variants.get(sv_type, set())
            if in_nyw and in_65:
                both_count += 1
            elif in_nyw:
                only_nyw_count += 1
            elif in_65:
                only_65_count += 1
            else:
                unique_count += 1
        calculated_total = both_count + only_nyw_count + only_65_count + unique_count
        if calculated_total != total_count:
            unique_count = total_count - (both_count + only_nyw_count + only_65_count)
            unique_count = max(0, unique_count)
        results.append({'SVLEN_BIN': bin_name, 'TOTAL_COUNT': total_count, 'UNIQUE': unique_count,
                        'BOTH': both_count, 'ONLY_NYW': only_nyw_count, 'ONLY_65': only_65_count})
    return pd.DataFrame(results)


def create_precision_line_plot_svlen(nyuwa_df: pd.DataFrame, kg_df: pd.DataFrame,
                                     sv_type: str, plot_type: str, output_file: str,
                                     nw_precision_ratios: Dict[str, float],
                                     kg_precision_ratios: Dict[str, float]):
    """
    创建单张Precision曲线图，高度增加一倍（原2.67英寸 → 5.33英寸），宽度保持10英寸。
    移除图例、坐标轴标题、图标题。
    """
    print(f"正在为 {sv_type} 创建Precision曲线图...")
    try:
        plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Liberation Sans']
        plt.rcParams['axes.unicode_minus'] = False
        # 高度增加一倍：16/3 ≈ 5.33英寸，宽度10英寸
        fig, ax = plt.subplots(figsize=(10, 3))
        bin_names = nyuwa_df['SVLEN_BIN'].tolist()
        if not bin_names:
            bin_names = kg_df['SVLEN_BIN'].tolist()
        x = np.arange(len(bin_names))
        nw_precision = [nw_precision_ratios.get(bin_name, 0.0) for bin_name in bin_names]
        kg_precision = [kg_precision_ratios.get(bin_name, 0.0) for bin_name in bin_names]
        ax.plot(x, nw_precision, 'r-o', linewidth=1.5, markersize=4, label='NyuWa')
        ax.plot(x, kg_precision, 'b-o', linewidth=1.5, markersize=4, label='1019')
        ax.set_xticks(x)
        ax.set_xticklabels(bin_names, rotation=45, ha='right', fontsize=10)
        ax.set_xlabel('')  # 移除x轴标签
        ax.set_ylabel('')  # 移除y轴标签
        ax.set_ylim(0, 1.05)
        # ax.legend()  # 移除图例
        ax.grid(True, alpha=0.3, linestyle='--')
        plt.tight_layout()
        plt.savefig(output_file, bbox_inches='tight', facecolor='white', format='pdf')
        print(f"Precision曲线图已保存到: {output_file}")
    except Exception as e:
        print(f"错误：为 {sv_type} 创建Precision曲线图失败 - {e}")


def save_svlen_detailed_statistics(all_stats: Dict[str, Dict[str, pd.DataFrame]],
                                  nw_precision_stats: Dict[str, Dict[str, float]],
                                  kg_precision_stats: Dict[str, Dict[str, float]],
                                  output_file: str):
    print("正在保存详细的SV长度分布统计结果（含Precision）...")
    try:
        with open(output_file, 'w') as f:
            f.write("SV_TYPE\tSVLEN_BIN\tDATASET\tTOTAL_COUNT\tUNIQUE\tSHARE_BOTH\tSHARE_1019\tSHARE_65\tSHARE_NYW\tPRECISION\n")
            for sv_type_plot, datasets in all_stats.items():
                sv_type = sv_type_plot.split('_')[0]
                if 'nyuwa' in datasets:
                    for _, row in datasets['nyuwa'].iterrows():
                        bin_name = row['SVLEN_BIN']
                        precision = nw_precision_stats.get(sv_type_plot, {}).get(bin_name, 0.0)
                        f.write(f"{sv_type}\t{bin_name}\tNyuWa\t{row['TOTAL_COUNT']}\t{row['UNIQUE']}\t{row['BOTH']}\t{row['ONLY_1019']}\t{row['ONLY_65']}\t0\t{precision:.4f}\n")
                if '1019' in datasets:
                    for _, row in datasets['1019'].iterrows():
                        bin_name = row['SVLEN_BIN']
                        precision = kg_precision_stats.get(sv_type_plot, {}).get(bin_name, 0.0)
                        f.write(f"{sv_type}\t{bin_name}\t1019\t{row['TOTAL_COUNT']}\t{row['UNIQUE']}\t{row['BOTH']}\t0\t{row['ONLY_65']}\t{row['ONLY_NYW']}\t{precision:.4f}\n")
        print(f"详细的SV长度分布统计结果（含Precision）已保存到: {output_file}")
    except Exception as e:
        print(f"错误：保存详细的SV长度分布统计结果失败 - {e}")


def analyze_merged_ins_dup_svlen_distribution(nw_all: pd.DataFrame, kg_all: pd.DataFrame,
                                            nw_1019_variants: Dict[str, Set[Tuple]],
                                            nw_65_variants: Dict[str, Set[Tuple]],
                                            kg_1019_nw_variants: Dict[str, Set[str]],
                                            kg_65_1019_variants: Dict[str, Set[str]],
                                            bins: List[Tuple[str, Tuple]]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    print("正在分析INS和DUP合并的SV长度分布...")
    nw_ins_dup = nw_all[nw_all['SV_TYPE'].isin(['INS', 'DUP'])].copy()
    kg_ins_dup = kg_all[kg_all['SV_TYPE'].isin(['INS', 'DUP'])].copy()
    nw_ins_dup.loc[:, 'SVLEN_BIN'] = nw_ins_dup['SVLEN'].apply(lambda x: assign_svlen_bin(x, bins))
    kg_ins_dup.loc[:, 'SVLEN_BIN'] = kg_ins_dup['SVLEN'].apply(lambda x: assign_svlen_bin(x, bins))
    nw_ins_dup.loc[:, 'composite_key'] = nw_ins_dup['CHROM'] + '_' + nw_ins_dup['POS'].astype(str) + '_' + nw_ins_dup['SV_TYPE'] + '_' + nw_ins_dup['SVLEN'].astype(str)
    nw_1019_keys = set()
    for sv_type in ['INS', 'DUP']:
        for variant in nw_1019_variants.get(sv_type, set()):
            chrom, pos, svtype, svlen = variant
            key = f"{chrom}_{pos}_{svtype}_{svlen}"
            nw_1019_keys.add(key)
    nw_65_keys = set()
    for sv_type in ['INS', 'DUP']:
        for variant in nw_65_variants.get(sv_type, set()):
            chrom, pos, svtype, svlen = variant
            key = f"{chrom}_{pos}_{svtype}_{svlen}"
            nw_65_keys.add(key)
    nw_results = []
    bin_names = [bin_name for bin_name, _ in bins]
    for bin_name in bin_names:
        bin_variants = nw_ins_dup[nw_ins_dup['SVLEN_BIN'] == bin_name]
        total_count = len(bin_variants)
        if total_count == 0:
            nw_results.append({'SVLEN_BIN': bin_name, 'TOTAL_COUNT': 0, 'UNIQUE': 0, 'BOTH': 0, 'ONLY_1019': 0, 'ONLY_65': 0})
            continue
        both_count = only_1019_count = only_65_count = unique_count = 0
        for _, variant in bin_variants.iterrows():
            key = variant['composite_key']
            in_1019 = key in nw_1019_keys
            in_65 = key in nw_65_keys
            if in_1019 and in_65:
                both_count += 1
            elif in_1019:
                only_1019_count += 1
            elif in_65:
                only_65_count += 1
            else:
                unique_count += 1
        calculated_total = both_count + only_1019_count + only_65_count + unique_count
        if calculated_total != total_count:
            unique_count = total_count - (both_count + only_1019_count + only_65_count)
            unique_count = max(0, unique_count)
        nw_results.append({'SVLEN_BIN': bin_name, 'TOTAL_COUNT': total_count, 'UNIQUE': unique_count,
                           'BOTH': both_count, 'ONLY_1019': only_1019_count, 'ONLY_65': only_65_count})
    kg_results = []
    for bin_name in bin_names:
        bin_variants = kg_ins_dup[kg_ins_dup['SVLEN_BIN'] == bin_name]
        total_count = len(bin_variants)
        if total_count == 0:
            kg_results.append({'SVLEN_BIN': bin_name, 'TOTAL_COUNT': 0, 'UNIQUE': 0, 'BOTH': 0, 'ONLY_NYW': 0, 'ONLY_65': 0})
            continue
        both_count = only_nyw_count = only_65_count = unique_count = 0
        for _, variant in bin_variants.iterrows():
            variant_id = variant['ID']
            in_nyw = any(variant_id in kg_1019_nw_variants.get(st, set()) for st in ['INS', 'DUP'])
            in_65 = any(variant_id in kg_65_1019_variants.get(st, set()) for st in ['INS', 'DUP'])
            if in_nyw and in_65:
                both_count += 1
            elif in_nyw:
                only_nyw_count += 1
            elif in_65:
                only_65_count += 1
            else:
                unique_count += 1
        calculated_total = both_count + only_nyw_count + only_65_count + unique_count
        if calculated_total != total_count:
            unique_count = total_count - (both_count + only_nyw_count + only_65_count)
            unique_count = max(0, unique_count)
        kg_results.append({'SVLEN_BIN': bin_name, 'TOTAL_COUNT': total_count, 'UNIQUE': unique_count,
                           'BOTH': both_count, 'ONLY_NYW': only_nyw_count, 'ONLY_65': only_65_count})
    return pd.DataFrame(nw_results), pd.DataFrame(kg_results)


def parse_bin_range(bin_name: str) -> Tuple[float, float]:
    if '-' in bin_name:
        parts = bin_name.split('-')
        return float(parts[0]), float(parts[1])
    elif bin_name.startswith('>'):
        return float(bin_name[1:]) + 1, float('inf')
    else:
        return 0, float('inf')


def append_summary_statistics(output_file: str, all_stats: Dict[str, Dict[str, pd.DataFrame]]):
    range_rules = {
        'DEL': [('≤200 bp', (0, 200)), ('200 bp - 4 kb', (201, 4000)), ('>4 kb', (4001, float('inf')))],
        'INS': [('≤100 bp', (0, 100)), ('>100 bp', (101, float('inf')))],
        'DUP': [('≤1 kb', (0, 1000)), ('>1 kb', (1001, float('inf')))],
        'INV': [('≤300 bp', (0, 300)), ('>300 bp', (301, float('inf')))],
        'SVA': [('≤400 bp', (0, 400)), ('400–2,000 bp', (401, 2000)), ('>2,000 bp', (2001, float('inf'))), ('All', (0, float('inf')))],
        'Alu': [('All', (0, float('inf')))],
        'LINE': [('All', (0, float('inf')))],
        'LTR': [('All', (0, float('inf')))]
    }
    lines = ["\n\n===== Additional Summary Statistics for Specified Length Ranges (Precision) =====\n", "SV_TYPE\tDATASET\tLENGTH_RANGE\tTOTAL_COUNT\tPRECISION\n"]
    for sv_type, ranges in range_rules.items():
        key = f"{sv_type}_default"
        if key not in all_stats:
            continue
        datasets = all_stats[key]
        for dataset_name, df in datasets.items():
            if df.empty:
                continue
            for range_name, (start, end) in ranges:
                total_in_stats = 0
                unique_in_stats = 0
                for _, row in df.iterrows():
                    bin_name = row['SVLEN_BIN']
                    bin_start, bin_end = parse_bin_range(bin_name)
                    if bin_start <= end and bin_end >= start:
                        total_in_stats += row['TOTAL_COUNT']
                        unique_in_stats += row['UNIQUE']
                precision = 1 - (unique_in_stats / total_in_stats) if total_in_stats > 0 else 0.0
                lines.append(f"{sv_type}\t{dataset_name}\t{range_name}\t{total_in_stats}\t{precision:.4f}\n")
    try:
        with open(output_file, 'a') as f:
            f.writelines(lines)
        print(f"已追加指定区间统计（Precision）到: {output_file}")
    except Exception as e:
        print(f"错误：追加统计结果失败 - {e}")


def main():
    kg_all_file = "/Parastor300s_G30S/haoxp/third/cmp/65/1019.svall.inv.addsvtype.rmsvtr.notduptoins_svlenge50.rm05mesv.rm05trsv_newsvtype.txt"
    nw_all_file = "/Parastor300s_G30S/haoxp/third/cmp/65/nyuwa.svall.addsvtype.svlenge50.rm05mesv.rm05trsv_newsvtype.txt"
    cmp_1019_nw_file = "/Parastor300s_G30S/haoxp/third/cmp/65/svcmp.1019.nw.sv.inv.over05.rmnan.addsvtype.rmdismatchINV.notduptoins_svlenge50.rm05mesv.rm05trsv.newsvtype.txt"
    cmp_65_nw_file = "/Parastor300s_G30S/haoxp/third/cmp/65/svcmp.65.nw.sv.inv.over05.rmnan.addsvtype.notduptoins_svlenge50.rm05mesv.rm05trsv.newsvtype.txt"
    cmp_65_1019_file = "/Parastor300s_G30S/haoxp/third/cmp/65/svcmp.65.1019.sv.inv.over05.rmnan.addsvtype.rmdismatchINV.notduptoins_svlenge50.rm05mesv.rm05trsv.newsvtype.txt"
    output_dir = "/Parastor300s_G30S/haoxp/third/cmp/65/svlen_alltype"
    output_stat_file = os.path.join(output_dir, "1kallsv_nyuwaallsv_65_over05_svlen_distribution_svlenge50.rm05mesv.rm05trsv_splitdup.detail.precision.line.txt")
    output_plot_pattern = os.path.join(output_dir, "1kallsv_nyuwaallsv_65_over05_svlen_distribution_svlenge50.rm05mesv.rm05trsv_splitdup.precision.line.{}.{}.pdf")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    print("开始分析...")
    kg_all = read_all_variants(kg_all_file, "1019")
    nw_all = read_all_variants(nw_all_file, "NyuWa")
    nw_1019_variants, kg_1019_nw_variants = read_cmp_1019_nw(cmp_1019_nw_file)
    nw_65_variants = read_cmp_65_nw(cmp_65_nw_file)
    kg_65_1019_variants = read_cmp_65_1019(cmp_65_1019_file)
    target_types = [('Alu', 'default'), ('Alu', 'detailed'), ('DEL', 'default'), ('INS', 'default'),
                    ('INV', 'default'), ('LINE', 'default'), ('LTR', 'default'), ('SVA', 'default'), ('DUP', 'default')]
    all_stats = {}
    nw_precision_all = {}
    kg_precision_all = {}
    for sv_type, plot_type in target_types:
        print(f"\n正在处理 {sv_type} ({plot_type})...")
        bins = create_svlen_bins(sv_type, plot_type)
        bin_names = [bn for bn, _ in bins]
        nyuwa_svlen_df = analyze_nyuwa_svlen_distribution(nw_all, nw_1019_variants, nw_65_variants, sv_type, bins)
        kg_svlen_df = analyze_1019_svlen_distribution(kg_all, kg_1019_nw_variants, kg_65_1019_variants, sv_type, bins)
        nw_precision = {}
        kg_precision = {}
        for bin_name in bin_names:
            row_nw = nyuwa_svlen_df[nyuwa_svlen_df['SVLEN_BIN'] == bin_name]
            total = row_nw.iloc[0]['TOTAL_COUNT'] if not row_nw.empty else 0
            unique = row_nw.iloc[0]['UNIQUE'] if not row_nw.empty else 0
            nw_precision[bin_name] = 1 - unique/total if total > 0 else 0.0
            row_kg = kg_svlen_df[kg_svlen_df['SVLEN_BIN'] == bin_name]
            total = row_kg.iloc[0]['TOTAL_COUNT'] if not row_kg.empty else 0
            unique = row_kg.iloc[0]['UNIQUE'] if not row_kg.empty else 0
            kg_precision[bin_name] = 1 - unique/total if total > 0 else 0.0
        all_stats[f"{sv_type}_{plot_type}"] = {'nyuwa': nyuwa_svlen_df, '1019': kg_svlen_df}
        nw_precision_all[f"{sv_type}_{plot_type}"] = nw_precision
        kg_precision_all[f"{sv_type}_{plot_type}"] = kg_precision
        if not nyuwa_svlen_df.empty and not kg_svlen_df.empty:
            plot_file = output_plot_pattern.format(sv_type, plot_type)
            create_precision_line_plot_svlen(nyuwa_svlen_df, kg_svlen_df, sv_type, plot_type, plot_file, nw_precision, kg_precision)
    print("\n正在处理 INSMERGE_DUP...")
    ins_dup_bins = create_svlen_bins('INS', 'default')
    bin_names = [bn for bn, _ in ins_dup_bins]
    nw_ins_dup_df, kg_ins_dup_df = analyze_merged_ins_dup_svlen_distribution(
        nw_all, kg_all, nw_1019_variants, nw_65_variants, kg_1019_nw_variants, kg_65_1019_variants, ins_dup_bins)
    nw_precision_ins_dup = {}
    kg_precision_ins_dup = {}
    for bin_name in bin_names:
        row_nw = nw_ins_dup_df[nw_ins_dup_df['SVLEN_BIN'] == bin_name]
        total = row_nw.iloc[0]['TOTAL_COUNT'] if not row_nw.empty else 0
        unique = row_nw.iloc[0]['UNIQUE'] if not row_nw.empty else 0
        nw_precision_ins_dup[bin_name] = 1 - unique/total if total > 0 else 0.0
        row_kg = kg_ins_dup_df[kg_ins_dup_df['SVLEN_BIN'] == bin_name]
        total = row_kg.iloc[0]['TOTAL_COUNT'] if not row_kg.empty else 0
        unique = row_kg.iloc[0]['UNIQUE'] if not row_kg.empty else 0
        kg_precision_ins_dup[bin_name] = 1 - unique/total if total > 0 else 0.0
    all_stats["INSMERGE_DUP_default"] = {'nyuwa': nw_ins_dup_df, '1019': kg_ins_dup_df}
    nw_precision_all["INSMERGE_DUP_default"] = nw_precision_ins_dup
    kg_precision_all["INSMERGE_DUP_default"] = kg_precision_ins_dup
    if not nw_ins_dup_df.empty and not kg_ins_dup_df.empty:
        ins_dup_plot_file = output_plot_pattern.format("INSMERGE_DUP", "default")
        create_precision_line_plot_svlen(nw_ins_dup_df, kg_ins_dup_df, "INSMERGE_DUP", "default", ins_dup_plot_file, nw_precision_ins_dup, kg_precision_ins_dup)
    save_svlen_detailed_statistics(all_stats, nw_precision_all, kg_precision_all, output_stat_file)
    append_summary_statistics(output_stat_file, all_stats)
    print("\n处理完成！")
    print(f"详细统计结果: {output_stat_file}")
    print(f"图表文件: {output_plot_pattern.format('*', '*')}")


if __name__ == "__main__":
    main()
