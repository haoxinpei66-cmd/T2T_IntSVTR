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
    'font.size': 34,
    'axes.titlesize': 30,
    'axes.labelsize': 38,
    'xtick.labelsize': 36,
    'ytick.labelsize': 36,
    'legend.fontsize': 34,
    'figure.titlesize': 32,
    'figure.dpi': 300,
})


def normalize_sv_type(sv_type: str) -> str:
    """
    统一变异类型命名
    """
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


def read_all_variants(file_path: str, dataset: str) -> pd.DataFrame:
    """
    读取所有变异文件，返回包含SV类型和长度的DataFrame
    """
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
    """
    读取1019和NyuWa比较文件，返回按类型分组的变异集合
    """
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
    """
    读取65和NyuWa比较文件，返回按类型分组的NyuWa变异集合
    """
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
    """
    读取65和1019比较文件，返回按类型分组的1019变异集合
    """
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
    """
    根据变异类型创建SV长度区间
    """
    bins = []

    if sv_type == 'Alu':
        if plot_type == "detailed":
            for i in range(50, 100, 10):
                start = i
                end = i + 9
                bins.append((f"{start}-{end}", (start, end)))
            for i in range(100, 400, 50):
                start = i + 1
                end = i + 50
                bins.append((f"{start}-{end}", (start, end)))
            for i in range(400, 1000, 100):
                start = i + 1
                end = i + 100
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
        bins.append(("50-60", (50, 60)))
        bins.append(("61-70", (61, 70)))
        bins.append(("71-80", (71, 80)))
        bins.append(("81-100", (81, 100)))
        bins.append(("101-200", (101, 200)))
        bins.append(("201-300", (201, 300)))
        bins.append(("301-400", (301, 400)))
        bins.append(("401-600", (401, 600)))
        bins.append(("601-800", (601, 800)))
        bins.append(("801-1000", (801, 1000)))
        bins.append(("1001-2000", (1001, 2000)))
        bins.append(("2001-4000", (2001, 4000)))
        bins.append(("4001-10000", (4001, 10000)))
        bins.append(("10001-20000", (10001, 20000)))
        bins.append(("20001-40000", (20001, 40000)))
        bins.append(("40001-100000", (40001, 100000)))
        bins.append((">100000", (100001, float('inf'))))

    elif sv_type in ['INS', 'DUP', 'INSMERGE_DUP']:
        bins.append(("50-60", (50, 60)))
        bins.append(("61-70", (61, 70)))
        bins.append(("71-80", (71, 80)))
        bins.append(("81-100", (81, 100)))
        bins.append(("101-200", (101, 200)))
        bins.append(("201-300", (201, 300)))
        bins.append(("301-400", (301, 400)))
        bins.append(("401-600", (401, 600)))
        bins.append(("601-1000", (601, 1000)))
        bins.append(("1001-2000", (1001, 2000)))
        bins.append(("2001-3000", (2001, 3000)))
        bins.append(("3001-4000", (3001, 4000)))
        bins.append(("4001-6000", (4001, 6000)))
        bins.append(("6001-10000", (6001, 10000)))
        bins.append((">10000", (10001, float('inf'))))

    elif sv_type == 'LTR':
        bins.append(("50-100", (50, 100)))
        bins.append(("101-200", (101, 200)))
        bins.append(("201-300", (201, 300)))
        bins.append(("301-400", (301, 400)))
        bins.append(("401-600", (401, 600)))
        bins.append(("601-1000", (601, 1000)))
        bins.append(("1001-2000", (1001, 2000)))
        bins.append(("2001-4000", (2001, 4000)))
        bins.append(("4001-10000", (4001, 10000)))
        bins.append((">10000", (10001, float('inf'))))

    else:  # LINE, INV, SVA
        bins.append(("50-75", (50, 75)))
        bins.append(("76-100", (76, 100)))
        bins.append(("101-200", (101, 200)))
        bins.append(("201-300", (201, 300)))
        bins.append(("301-400", (301, 400)))
        bins.append(("401-600", (401, 600)))
        bins.append(("601-1000", (601, 1000)))
        bins.append(("1001-2000", (1001, 2000)))
        bins.append(("2001-3000", (2001, 3000)))
        bins.append(("3001-4000", (3001, 4000)))
        bins.append(("4001-6000", (4001, 6000)))
        bins.append(("6001-10000", (6001, 10000)))
        bins.append((">10000", (10001, float('inf'))))

    return bins


def assign_svlen_bin(svlen: float, bins: List[Tuple[str, Tuple]]) -> str:
    """
    将SV长度分配到对应的区间
    """
    for bin_name, (start, end) in bins:
        if start <= svlen <= end:
            return bin_name
    return "Other"


def analyze_nyuwa_svlen_distribution(nw_all: pd.DataFrame,
                                   nw_1019_variants: Dict[str, Set[Tuple]],
                                   nw_65_variants: Dict[str, Set[Tuple]],
                                   sv_type: str, bins: List[Tuple[str, Tuple]]) -> pd.DataFrame:
    """
    分析NyuWa变异在不同比较文件中的重叠关系，按SV长度区间分组
    """
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
            results.append({
                'SVLEN_BIN': bin_name,
                'TOTAL_COUNT': 0,
                'UNIQUE': 0,
                'BOTH': 0,
                'ONLY_1019': 0,
                'ONLY_65': 0
            })
            continue

        both_count = 0
        only_1019_count = 0
        only_65_count = 0
        unique_count = 0

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

        results.append({
            'SVLEN_BIN': bin_name,
            'TOTAL_COUNT': total_count,
            'UNIQUE': unique_count,
            'BOTH': both_count,
            'ONLY_1019': only_1019_count,
            'ONLY_65': only_65_count
        })

    return pd.DataFrame(results)


def analyze_1019_svlen_distribution(kg_all: pd.DataFrame,
                                  kg_1019_nw_variants: Dict[str, Set[str]],
                                  kg_65_1019_variants: Dict[str, Set[str]],
                                  sv_type: str, bins: List[Tuple[str, Tuple]]) -> pd.DataFrame:
    """
    分析1019变异在不同比较文件中的重叠关系，按SV长度区间分组
    """
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
            results.append({
                'SVLEN_BIN': bin_name,
                'TOTAL_COUNT': 0,
                'UNIQUE': 0,
                'BOTH': 0,
                'ONLY_NYW': 0,
                'ONLY_65': 0
            })
            continue

        both_count = 0
        only_nyw_count = 0
        only_65_count = 0
        unique_count = 0

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

        results.append({
            'SVLEN_BIN': bin_name,
            'TOTAL_COUNT': total_count,
            'UNIQUE': unique_count,
            'BOTH': both_count,
            'ONLY_NYW': only_nyw_count,
            'ONLY_65': only_65_count
        })

    return pd.DataFrame(results)


def create_svlen_stacked_plot(nyuwa_df: pd.DataFrame, kg_df: pd.DataFrame,
                             sv_type: str, plot_type: str, output_file: str,
                             nw_recall_ratios: Optional[Dict[str, float]] = None,
                             kg_recall_ratios: Optional[Dict[str, float]] = None):
    """
    创建按SV长度区间的堆叠柱状图，并添加Recall曲线
    """
    print(f"正在为 {sv_type} 创建SV长度分布堆叠柱状图（含Recall）...")

    try:
        plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Liberation Sans']
        plt.rcParams['axes.unicode_minus'] = False

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10))

        bin_names = nyuwa_df['SVLEN_BIN'].tolist()
        x = np.arange(len(bin_names))
        width = 0.7

        max_count = max(nyuwa_df['TOTAL_COUNT'].max(), kg_df['TOTAL_COUNT'].max())
        y_max = max_count * 1.15 if max_count > 0 else 100

        colors = {
            'UNIQUE': 'lightgray',
            'BOTH': 'gold',
            'ONLY_1019': 'lightcoral',
            'ONLY_65': 'lightblue',
            'ONLY_NYW': 'lightgreen'
        }

        # NyuWa图
        bottom = np.zeros(len(bin_names))
        parts_nw = ['BOTH', 'ONLY_1019', 'ONLY_65', 'UNIQUE']
        labels_nw = ['Share with both', 'Share with 1019', 'Share with 65', 'Unique']

        for i, part in enumerate(parts_nw):
            if part == 'BOTH':
                values = nyuwa_df['BOTH'].values
            elif part == 'ONLY_1019':
                values = nyuwa_df['ONLY_1019'].values
            elif part == 'ONLY_65':
                values = nyuwa_df['ONLY_65'].values
            else:
                values = nyuwa_df['UNIQUE'].values

            ax1.bar(x, values, width, bottom=bottom, color=colors[part],
                   label=labels_nw[i], edgecolor='white', linewidth=0.5)
            bottom += values

        ax1.set_ylabel('NyuWa Variant Count', fontsize=12, fontweight='bold')
        ax1.set_ylim(0, y_max)
        ax1.set_title(f'NyuWa {sv_type} SV Length Distribution - Overlap with 65 and 1019',
                     fontsize=14, fontweight='bold')
        ax1.set_xticks(x)
        ax1.set_xticklabels(bin_names, rotation=45, ha='right', fontsize=8)
        ax1.legend(loc='upper right')

        for i, total in enumerate(nyuwa_df['TOTAL_COUNT']):
            if total > 0:
                ax1.text(x[i], total + y_max * 0.01, f'{int(total):,}',
                        ha='center', va='bottom', fontsize=7, fontweight='bold')

        # 添加Recall曲线（右轴）
        if nw_recall_ratios is not None:
            ax1_right = ax1.twinx()
            recall_values = [nw_recall_ratios.get(bin_name, 0.0) for bin_name in bin_names]
            ax1_right.plot(x, recall_values, 'r-o', linewidth=2, markersize=6, label='Recall (65 dataset)')
            ax1_right.set_ylabel('Recall', fontsize=12, fontweight='bold', color='red')
            ax1_right.set_ylim(0, 1.05)
            ax1_right.tick_params(axis='y', labelcolor='red')
            ax1_right.legend(loc='upper left', fontsize=10)

        # 1019图
        bottom = np.zeros(len(bin_names))
        parts_kg = ['BOTH', 'ONLY_NYW', 'ONLY_65', 'UNIQUE']
        labels_kg = ['Share with both', 'Share with NyuWa', 'Share with 65', 'Unique']

        for i, part in enumerate(parts_kg):
            if part == 'BOTH':
                values = kg_df['BOTH'].values
            elif part == 'ONLY_NYW':
                values = kg_df['ONLY_NYW'].values
            elif part == 'ONLY_65':
                values = kg_df['ONLY_65'].values
            else:
                values = kg_df['UNIQUE'].values

            ax2.bar(x, values, width, bottom=bottom, color=colors[part],
                   label=labels_kg[i], edgecolor='white', linewidth=0.5)
            bottom += values

        ax2.set_ylabel('1019 Variant Count', fontsize=12, fontweight='bold')
        ax2.set_ylim(0, y_max)
        ax2.set_title(f'1019 {sv_type} SV Length Distribution - Overlap with 65 and NyuWa',
                     fontsize=14, fontweight='bold')
        ax2.set_xticks(x)
        ax2.set_xticklabels(bin_names, rotation=45, ha='right', fontsize=8)
        ax2.legend(loc='upper right')

        for i, total in enumerate(kg_df['TOTAL_COUNT']):
            if total > 0:
                ax2.text(x[i], total + y_max * 0.01, f'{int(total):,}',
                        ha='center', va='bottom', fontsize=7, fontweight='bold')

        # 添加Recall曲线（右轴）
        if kg_recall_ratios is not None:
            ax2_right = ax2.twinx()
            recall_values = [kg_recall_ratios.get(bin_name, 0.0) for bin_name in bin_names]
            ax2_right.plot(x, recall_values, 'r-o', linewidth=2, markersize=6, label='Recall (65 dataset)')
            ax2_right.set_ylabel('Recall', fontsize=12, fontweight='bold', color='red')
            ax2_right.set_ylim(0, 1.05)
            ax2_right.tick_params(axis='y', labelcolor='red')
            ax2_right.legend(loc='upper left', fontsize=10)

        ax1.grid(True, alpha=0.3, axis='y')
        ax2.grid(True, alpha=0.3, axis='y')
        plt.tight_layout()
        plt.savefig(output_file, bbox_inches='tight', facecolor='white', format='pdf')
        print(f"SV长度分布堆叠柱状图（含Recall）已保存到: {output_file}")

    except Exception as e:
        print(f"错误：为 {sv_type} 创建SV长度分布堆叠柱状图失败 - {e}")


def save_svlen_detailed_statistics(all_stats: Dict[str, Dict[str, pd.DataFrame]],
                                  nw_recall_stats: Dict[str, Dict[str, float]],
                                  kg_recall_stats: Dict[str, Dict[str, float]],
                                  output_file: str):
    """
    保存详细的SV长度分布统计结果，包含Recall列
    """
    print("正在保存详细的SV长度分布统计结果（含Recall）...")

    try:
        with open(output_file, 'w') as f:
            f.write("SV_TYPE\tSVLEN_BIN\tDATASET\tTOTAL_COUNT\tUNIQUE\tSHARE_BOTH\tSHARE_1019\tSHARE_65\tSHARE_NYW\tRECALL\n")

            for sv_type_plot, datasets in all_stats.items():
                sv_type = sv_type_plot.split('_')[0]  # 提取变异类型名称

                # 写入NyuWa数据
                if 'nyuwa' in datasets:
                    for _, row in datasets['nyuwa'].iterrows():
                        bin_name = row['SVLEN_BIN']
                        recall = nw_recall_stats.get(sv_type_plot, {}).get(bin_name, 0.0)
                        f.write(f"{sv_type}\t{bin_name}\tNyuWa\t{row['TOTAL_COUNT']}\t{row['UNIQUE']}\t{row['BOTH']}\t{row['ONLY_1019']}\t{row['ONLY_65']}\t0\t{recall:.4f}\n")

                # 写入1019数据
                if '1019' in datasets:
                    for _, row in datasets['1019'].iterrows():
                        bin_name = row['SVLEN_BIN']
                        recall = kg_recall_stats.get(sv_type_plot, {}).get(bin_name, 0.0)
                        f.write(f"{sv_type}\t{bin_name}\t1019\t{row['TOTAL_COUNT']}\t{row['UNIQUE']}\t{row['BOTH']}\t0\t{row['ONLY_65']}\t{row['ONLY_NYW']}\t{recall:.4f}\n")

        print(f"详细的SV长度分布统计结果（含Recall）已保存到: {output_file}")

    except Exception as e:
        print(f"错误：保存详细的SV长度分布统计结果失败 - {e}")


def analyze_merged_ins_dup_svlen_distribution(nw_all: pd.DataFrame, kg_all: pd.DataFrame,
                                            nw_1019_variants: Dict[str, Set[Tuple]],
                                            nw_65_variants: Dict[str, Set[Tuple]],
                                            kg_1019_nw_variants: Dict[str, Set[str]],
                                            kg_65_1019_variants: Dict[str, Set[str]],
                                            bins: List[Tuple[str, Tuple]]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    分析INS和DUP合并的SV长度分布
    """
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
            nw_results.append({
                'SVLEN_BIN': bin_name,
                'TOTAL_COUNT': 0,
                'UNIQUE': 0,
                'BOTH': 0,
                'ONLY_1019': 0,
                'ONLY_65': 0
            })
            continue

        both_count = 0
        only_1019_count = 0
        only_65_count = 0
        unique_count = 0

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

        nw_results.append({
            'SVLEN_BIN': bin_name,
            'TOTAL_COUNT': total_count,
            'UNIQUE': unique_count,
            'BOTH': both_count,
            'ONLY_1019': only_1019_count,
            'ONLY_65': only_65_count
        })

    kg_results = []
    for bin_name in bin_names:
        bin_variants = kg_ins_dup[kg_ins_dup['SVLEN_BIN'] == bin_name]
        total_count = len(bin_variants)

        if total_count == 0:
            kg_results.append({
                'SVLEN_BIN': bin_name,
                'TOTAL_COUNT': 0,
                'UNIQUE': 0,
                'BOTH': 0,
                'ONLY_NYW': 0,
                'ONLY_65': 0
            })
            continue

        both_count = 0
        only_nyw_count = 0
        only_65_count = 0
        unique_count = 0

        for _, variant in bin_variants.iterrows():
            variant_id = variant['ID']
            in_nyw = False
            in_65 = False

            for sv_type in ['INS', 'DUP']:
                if variant_id in kg_1019_nw_variants.get(sv_type, set()):
                    in_nyw = True
                if variant_id in kg_65_1019_variants.get(sv_type, set()):
                    in_65 = True

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

        kg_results.append({
            'SVLEN_BIN': bin_name,
            'TOTAL_COUNT': total_count,
            'UNIQUE': unique_count,
            'BOTH': both_count,
            'ONLY_NYW': only_nyw_count,
            'ONLY_65': only_65_count
        })

    return pd.DataFrame(nw_results), pd.DataFrame(kg_results)


# ==================== 新增：解析区间名称辅助函数 ====================
def parse_bin_range(bin_name: str) -> Tuple[float, float]:
    """
    解析区间名称，返回 (start, end)
    支持格式: "a-b" 或 ">X"
    """
    if '-' in bin_name:
        parts = bin_name.split('-')
        start = float(parts[0])
        end = float(parts[1])
        return start, end
    elif bin_name.startswith('>'):
        start = float(bin_name[1:]) + 1
        end = float('inf')
        return start, end
    else:
        # 其他情况（理论上不会出现）
        return 0, float('inf')


def append_summary_statistics(output_file: str, all_stats: Dict[str, Dict[str, pd.DataFrame]],
                              recall_df: pd.DataFrame, cmp_65_nw_ids: Set[str], cmp_65_1019_ids: Set[str]):
    """
    在原输出文件末尾追加指定长度区间的 Recall 和 Accuracy 统计
    recall_df: 16.sv.me.inv.info.txt 的 DataFrame
    cmp_65_nw_ids: 65-NyuWa 比较文件的ID集合
    cmp_65_1019_ids: 65-1019 比较文件的ID集合
    """
    # 定义要处理的变异类型和区间规则
    range_rules = {
        'DEL': [
            ('≤200 bp', (0, 200)),
            ('200 bp - 4 kb', (201, 4000)),
            ('>4 kb', (4001, float('inf')))
        ],
        'INS': [
            ('≤100 bp', (0, 100)),
            ('>100 bp', (101, float('inf')))
        ],
        'DUP': [
            ('≤1 kb', (0, 1000)),
            ('>1 kb', (1001, float('inf')))
        ],
        'INV': [
            ('≤300 bp', (0, 300)),
            ('>300 bp', (301, float('inf')))
        ],
        'SVA': [
            ('≤400 bp', (0, 400)),
            ('400–2,000 bp', (401, 2000)),
            ('>2,000 bp', (2001, float('inf'))),
            ('All', (0, float('inf')))   # 新增整体
        ],
        'Alu': [('All', (0, float('inf')))],
        'LINE': [('All', (0, float('inf')))],
        'LTR': [('All', (0, float('inf')))]
    }

    # 收集所有输出行
    lines = []
    lines.append("\n\n===== Additional Summary Statistics for Specified Length Ranges =====\n")
    lines.append("SV_TYPE\tDATASET\tLENGTH_RANGE\tTOTAL_COUNT\tRECALL\tACCURACY\n")

    for sv_type, ranges in range_rules.items():
        # 获取对应的统计键（default）
        key = f"{sv_type}_default"
        if key not in all_stats:
            continue

        datasets = all_stats[key]
        for dataset_name, df in datasets.items():
            if df.empty:
                continue

            # 根据数据集名称决定用哪个比较ID集合
            if dataset_name == 'nyuwa':
                cmp_ids = cmp_65_nw_ids
            elif dataset_name == '1019':
                cmp_ids = cmp_65_1019_ids
            else:
                continue

            # 从 recall_df 中提取该类型的所有变异
            recall_sub = recall_df[recall_df['type'] == sv_type]
            if recall_sub.empty:
                # 如果没有，则所有区间 TOTAL_COUNT = 0
                for range_name, (start, end) in ranges:
                    lines.append(f"{sv_type}\t{dataset_name}\t{range_name}\t0\t0.0\t0.0\n")
                continue

            for range_name, (start, end) in ranges:
                # 过滤落在该区间的 recall 变异
                bin_variants = recall_sub[(recall_sub['length'] >= start) & (recall_sub['length'] <= end)]
                total = len(bin_variants)
                if total == 0:
                    recall = 0.0
                    accuracy = 0.0
                else:
                    # 分子：ID 出现在比较文件中的数量
                    matched = sum(1 for vid in bin_variants['id'].astype(str).str.strip() if vid in cmp_ids)
                    recall = matched / total
                    # 计算 accuracy: 从 all_stats 中获取该区间内该数据集的 TOTAL 和 UNIQUE
                    # 根据已有的 df（按详细bin分组）计算该新区间内的统计
                    total_in_stats = 0
                    unique_in_stats = 0
                    for _, row in df.iterrows():
                        bin_name = row['SVLEN_BIN']
                        bin_start, bin_end = parse_bin_range(bin_name)
                        # 检查该bin是否与新区间重叠
                        if bin_start <= end and bin_end >= start:
                            total_in_stats += row['TOTAL_COUNT']
                            unique_in_stats += row['UNIQUE']
                    if total_in_stats == 0:
                        accuracy = 0.0
                    else:
                        accuracy = 1 - (unique_in_stats / total_in_stats)
                lines.append(f"{sv_type}\t{dataset_name}\t{range_name}\t{total}\t{recall:.4f}\t{accuracy:.4f}\n")

    # 追加到文件末尾
    try:
        with open(output_file, 'a') as f:
            f.writelines(lines)
        print(f"已追加指定区间统计到: {output_file}")
    except Exception as e:
        print(f"错误：追加统计结果失败 - {e}")


def main():
    # 文件路径配置
    kg_all_file = "/Parastor300s_G30S/haoxp/third/cmp/65/1019.svall.inv.addsvtype.rmsvtr.notduptoins_svlenge50.rm05mesv.rm05trsv_newsvtype.txt"
    nw_all_file = "/Parastor300s_G30S/haoxp/third/cmp/65/nyuwa.svall.addsvtype.svlenge50.rm05mesv.rm05trsv_newsvtype.txt"

    cmp_1019_nw_file = "/Parastor300s_G30S/haoxp/third/cmp/65/svcmp.1019.nw.sv.inv.over05.rmnan.addsvtype.rmdismatchINV.notduptoins_svlenge50.rm05mesv.rm05trsv.newsvtype.txt"
    cmp_65_nw_file = "/Parastor300s_G30S/haoxp/third/cmp/65/svcmp.65.nw.sv.inv.over05.rmnan.addsvtype.notduptoins_svlenge50.rm05mesv.rm05trsv.newsvtype.txt"
    cmp_65_1019_file = "/Parastor300s_G30S/haoxp/third/cmp/65/svcmp.65.1019.sv.inv.over05.rmnan.addsvtype.rmdismatchINV.notduptoins_svlenge50.rm05mesv.rm05trsv.newsvtype.txt"
    recall_file = "/Parastor300s_G30S/haoxp/third/cmp/65/16.sv.me.inv.info.txt"

    output_dir = "/Parastor300s_G30S/haoxp/third/cmp/65/svlen_alltype"
    output_stat_file = os.path.join(output_dir, "1kallsv_nyuwaallsv_65_over05_svlen_distribution_svlenge50.rm05mesv.rm05trsv_splitdup.detail.recall.txt")
    output_plot_pattern = os.path.join(output_dir, "1kallsv_nyuwaallsv_65_over05_svlen_distribution_svlenge50.rm05mesv.rm05trsv_splitdup.recall.{}.{}.pdf")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"创建输出目录: {output_dir}")

    print("开始分析多数据集SV变异类型SV长度分布重叠关系（含Recall）...")

    # 1. 读取所有变异文件
    kg_all = read_all_variants(kg_all_file, "1019")
    nw_all = read_all_variants(nw_all_file, "NyuWa")

    # 2. 读取比较文件
    nw_1019_variants, kg_1019_nw_variants = read_cmp_1019_nw(cmp_1019_nw_file)
    nw_65_variants = read_cmp_65_nw(cmp_65_nw_file)
    kg_65_1019_variants = read_cmp_65_1019(cmp_65_1019_file)

    # 3. 加载 recall 文件并提取65比较文件的ID集合
    print("加载 recall 文件和65比较文件ID集合...")
    if not os.path.exists(recall_file):
        print(f"错误：recall文件不存在 {recall_file}")
        sys.exit(1)

    recall_df = pd.read_csv(recall_file, sep='\t', comment='#')
    # 检查列名，可能有无头情况
    if 'type' not in recall_df.columns or 'id' not in recall_df.columns:
        print("警告：recall文件列名不符合预期，尝试重新读取（无头）")
        recall_df = pd.read_csv(recall_file, sep='\t', header=None,
                                names=['chr','pos','end','id','length','type','direction','region'])
    recall_df['type'] = recall_df['type'].apply(normalize_sv_type)
    recall_df['length'] = pd.to_numeric(recall_df['length'], errors='coerce')
    recall_df = recall_df.dropna(subset=['length', 'type'])

    # 读取65-NyuWa比较文件（第三列ID）
    cmp_65_nw_ids = set()
    if os.path.exists(cmp_65_nw_file):
        df_nw = pd.read_csv(cmp_65_nw_file, sep='\t', comment='#', header=None)
        if len(df_nw.columns) >= 3:
            cmp_65_nw_ids = set(df_nw[2].astype(str).str.strip())
    else:
        print(f"警告：文件不存在 {cmp_65_nw_file}")

    # 读取65-1019比较文件（第十一列ID，索引10）
    cmp_65_1019_ids = set()
    if os.path.exists(cmp_65_1019_file):
        df_1019 = pd.read_csv(cmp_65_1019_file, sep='\t', comment='#', header=None)
        if len(df_1019.columns) >= 11:
            cmp_65_1019_ids = set(df_1019[10].astype(str).str.strip())
    else:
        print(f"警告：文件不存在 {cmp_65_1019_file}")

    # 4. 定义需要统计的变异类型和对应的图表类型
    target_types = [
        ('Alu', 'default'),
        ('Alu', 'detailed'),
        ('DEL', 'default'),
        ('INS', 'default'),
        ('INV', 'default'),
        ('LINE', 'default'),
        ('LTR', 'default'),
        ('SVA', 'default'),
        ('DUP', 'default')
    ]

    # 5. 存储所有统计结果和Recall数据
    all_stats = {}
    nw_recall_all = {}
    kg_recall_all = {}

    for sv_type, plot_type in target_types:
        print(f"\n正在处理 {sv_type} 变异类型 ({plot_type})...")

        # 创建对应的区间划分
        bins = create_svlen_bins(sv_type, plot_type)
        bin_names = [bin_name for bin_name, _ in bins]

        # 计算 recall 比例
        recall_sub = recall_df[recall_df['type'] == sv_type].copy()
        if recall_sub.empty:
            print(f"警告：recall文件中没有 {sv_type} 变异")
            nw_recall = {bin_name: 0.0 for bin_name in bin_names}
            kg_recall = {bin_name: 0.0 for bin_name in bin_names}
        else:
            recall_sub['SVLEN_BIN'] = recall_sub['length'].apply(lambda x: assign_svlen_bin(x, bins))
            nw_recall = {}
            kg_recall = {}
            for bin_name in bin_names:
                bin_variants = recall_sub[recall_sub['SVLEN_BIN'] == bin_name]
                total = len(bin_variants)
                if total == 0:
                    nw_recall[bin_name] = 0.0
                    kg_recall[bin_name] = 0.0
                else:
                    matched_nw = sum(1 for vid in bin_variants['id'].astype(str).str.strip() if vid in cmp_65_nw_ids)
                    nw_recall[bin_name] = matched_nw / total
                    matched_1019 = sum(1 for vid in bin_variants['id'].astype(str).str.strip() if vid in cmp_65_1019_ids)
                    kg_recall[bin_name] = matched_1019 / total

        # 计算NyuWa SV长度分布
        nyuwa_svlen_df = analyze_nyuwa_svlen_distribution(
            nw_all, nw_1019_variants, nw_65_variants, sv_type, bins
        )

        # 计算1019 SV长度分布
        kg_svlen_df = analyze_1019_svlen_distribution(
            kg_all, kg_1019_nw_variants, kg_65_1019_variants, sv_type, bins
        )

        # 保存统计结果
        all_stats[f"{sv_type}_{plot_type}"] = {
            'nyuwa': nyuwa_svlen_df,
            '1019': kg_svlen_df
        }
        nw_recall_all[f"{sv_type}_{plot_type}"] = nw_recall
        kg_recall_all[f"{sv_type}_{plot_type}"] = kg_recall

        # 创建SV长度分布堆叠柱状图
        if not nyuwa_svlen_df.empty and not kg_svlen_df.empty:
            plot_file = output_plot_pattern.format(sv_type, plot_type)
            create_svlen_stacked_plot(nyuwa_svlen_df, kg_svlen_df, sv_type, plot_type, plot_file,
                                      nw_recall, kg_recall)
        else:
            print(f"警告：{sv_type} 没有数据，跳过绘图")

    # 6. 计算INS和DUP合并的SV长度分布
    print(f"\n正在处理 INSMERGE_DUP 变异类型...")

    ins_dup_bins = create_svlen_bins('INS', 'default')
    bin_names = [bin_name for bin_name, _ in ins_dup_bins]

    # 计算 recall 比例（INS和DUP合并）
    recall_ins_dup = recall_df[recall_df['type'].isin(['INS', 'DUP'])].copy()
    recall_ins_dup['SVLEN_BIN'] = recall_ins_dup['length'].apply(lambda x: assign_svlen_bin(x, ins_dup_bins))

    nw_recall_ins_dup = {bin_name: 0.0 for bin_name in bin_names}
    kg_recall_ins_dup = {bin_name: 0.0 for bin_name in bin_names}
    for bin_name in bin_names:
        bin_variants = recall_ins_dup[recall_ins_dup['SVLEN_BIN'] == bin_name]
        total = len(bin_variants)
        if total == 0:
            nw_recall_ins_dup[bin_name] = 0.0
            kg_recall_ins_dup[bin_name] = 0.0
        else:
            matched_nw = sum(1 for vid in bin_variants['id'].astype(str).str.strip() if vid in cmp_65_nw_ids)
            nw_recall_ins_dup[bin_name] = matched_nw / total
            matched_1019 = sum(1 for vid in bin_variants['id'].astype(str).str.strip() if vid in cmp_65_1019_ids)
            kg_recall_ins_dup[bin_name] = matched_1019 / total

    # 计算INS和DUP合并的SV长度分布
    nw_ins_dup_df, kg_ins_dup_df = analyze_merged_ins_dup_svlen_distribution(
        nw_all, kg_all, nw_1019_variants, nw_65_variants,
        kg_1019_nw_variants, kg_65_1019_variants, ins_dup_bins
    )

    # 保存统计结果
    all_stats["INSMERGE_DUP_default"] = {
        'nyuwa': nw_ins_dup_df,
        '1019': kg_ins_dup_df
    }
    nw_recall_all["INSMERGE_DUP_default"] = nw_recall_ins_dup
    kg_recall_all["INSMERGE_DUP_default"] = kg_recall_ins_dup

    # 创建INS+DUP合并图表
    if not nw_ins_dup_df.empty and not kg_ins_dup_df.empty:
        ins_dup_plot_file = output_plot_pattern.format("INSMERGE_DUP", "default")
        create_svlen_stacked_plot(nw_ins_dup_df, kg_ins_dup_df, "INSMERGE_DUP", "default", ins_dup_plot_file,
                                  nw_recall_ins_dup, kg_recall_ins_dup)
    else:
        print("警告：INSMERGE_DUP 没有数据，跳过绘图")

    # 7. 保存详细的SV长度分布统计结果（含Recall）
    save_svlen_detailed_statistics(all_stats, nw_recall_all, kg_recall_all, output_stat_file)

    # 8. 追加指定区间统计（新增）
    append_summary_statistics(output_stat_file, all_stats, recall_df, cmp_65_nw_ids, cmp_65_1019_ids)

    print("\n处理完成！")
    print(f"详细统计结果: {output_stat_file}")
    print(f"图表文件: {output_plot_pattern.format('*', '*')}")


if __name__ == "__main__":
    main()
