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


def compute_recall_ratios(recall_file: str,
                          cmp_65_nw_file: str,
                          cmp_65_1019_file: str,
                          sv_types: List[str]) -> Tuple[Dict[str, float], Dict[str, float]]:
    """
    计算每个变异类型的Recall比例
    返回两个字典：nw_recall (NyuWa图用) 和 kg_recall (1019图用)
    """
    print("正在计算Recall比例...")
    # 读取 recall_file (16.sv.me.inv.info.txt)
    if not os.path.exists(recall_file):
        print(f"错误：文件不存在 {recall_file}")
        return {sv: 0.0 for sv in sv_types}, {sv: 0.0 for sv in sv_types}
    recall_df = pd.read_csv(recall_file, sep='\t', comment='#')
    # 检查列名，可能有无头情况
    if 'type' not in recall_df.columns or 'id' not in recall_df.columns:
        recall_df = pd.read_csv(recall_file, sep='\t', header=None,
                                names=['chr','pos','end','id','length','type','direction','region'])
    recall_df['type'] = recall_df['type'].apply(normalize_sv_type)
    # 读取65-NyuWa比较文件（第三列ID）
    cmp_65_nw_ids = set()
    if os.path.exists(cmp_65_nw_file):
        df_nw = pd.read_csv(cmp_65_nw_file, sep='\t', comment='#', header=None)
        if len(df_nw.columns) >= 3:
            cmp_65_nw_ids = set(df_nw[2].astype(str).str.strip())
    # 读取65-1019比较文件（第十一列ID）
    cmp_65_1019_ids = set()
    if os.path.exists(cmp_65_1019_file):
        df_1019 = pd.read_csv(cmp_65_1019_file, sep='\t', comment='#', header=None)
        if len(df_1019.columns) >= 11:
            cmp_65_1019_ids = set(df_1019[10].astype(str).str.strip())
    # 计算每个变异类型的Recall
    nw_recall = {}
    kg_recall = {}
    for sv_type in sv_types:
        total = len(recall_df[recall_df['type'] == sv_type])
        if total == 0:
            nw_recall[sv_type] = 0.0
            kg_recall[sv_type] = 0.0
            continue
        ids_in_group = recall_df[recall_df['type'] == sv_type]['id'].astype(str).str.strip()
        matched_nw = sum(1 for vid in ids_in_group if vid in cmp_65_nw_ids)
        matched_1019 = sum(1 for vid in ids_in_group if vid in cmp_65_1019_ids)
        nw_recall[sv_type] = matched_nw / total
        kg_recall[sv_type] = matched_1019 / total
        print(f"  {sv_type}: total={total}, recall_nw={matched_nw}/{total}={nw_recall[sv_type]:.4f}, recall_kg={matched_1019}/{total}={kg_recall[sv_type]:.4f}")
    return nw_recall, kg_recall


def save_detailed_statistics_with_recall(nyuwa_df: pd.DataFrame, kg_df: pd.DataFrame,
                                         nw_recall: Dict[str, float], kg_recall: Dict[str, float],
                                         output_file: str):
    """保存详细的统计结果，包含Recall列"""
    print("正在保存详细统计结果（含Recall）...")
    try:
        with open(output_file, 'w') as f:
            f.write("SV_TYPE\tDATASET\tTOTAL_COUNT\tUNIQUE\tSHARE_BOTH\tSHARE_1019\tSHARE_65\tSHARE_NYW\tRECALL\n")
            # NyuWa数据
            for _, row in nyuwa_df.iterrows():
                sv_type = row['SV_TYPE']
                recall = nw_recall.get(sv_type, 0.0)
                f.write(f"{sv_type}\t{row['DATASET']}\t{row['TOTAL_COUNT']}\t{row['UNIQUE']}\t{row['BOTH']}\t{row['ONLY_1019']}\t{row['ONLY_65']}\t0\t{recall:.4f}\n")
            # 1019数据
            for _, row in kg_df.iterrows():
                sv_type = row['SV_TYPE']
                recall = kg_recall.get(sv_type, 0.0)
                f.write(f"{sv_type}\t{row['DATASET']}\t{row['TOTAL_COUNT']}\t{row['UNIQUE']}\t{row['BOTH']}\t0\t{row['ONLY_65']}\t{row['ONLY_NYW']}\t{recall:.4f}\n")
        print(f"详细统计结果（含Recall）已保存到: {output_file}")
    except Exception as e:
        print(f"错误：保存详细统计结果失败 - {e}")


def create_stacked_bar_plot_with_recall(nyuwa_df: pd.DataFrame, kg_df: pd.DataFrame,
                                        nw_recall: Dict[str, float], kg_recall: Dict[str, float],
                                        output_file: str):
    """创建堆叠柱状图，并在右轴添加Recall曲线"""
    print("正在创建堆叠柱状图（含Recall曲线）...")
    try:
        plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Liberation Sans']
        plt.rcParams['axes.unicode_minus'] = False
        sv_types = ['DEL', 'INS', 'DUP', 'INV', 'Alu', 'LINE', 'LTR', 'SVA']
        nyuwa_df = nyuwa_df.set_index('SV_TYPE').reindex(sv_types).reset_index()
        kg_df = kg_df.set_index('SV_TYPE').reindex(sv_types).reset_index()
        fig_width = 10
        fig_height = 10
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(fig_width, fig_height))
        x = np.arange(len(sv_types))
        width = 0.5
        x_compression = 0.85
        x_adjusted = x * x_compression
        max_count = max(nyuwa_df['TOTAL_COUNT'].max(), kg_df['TOTAL_COUNT'].max())
        y_max = max_count * 1.15 if max_count > 0 else 100
        colors = {
            'UNIQUE': 'lightgray',
            'SHARE_BOTH': 'gold',
            'SHARE_1019': 'lightcoral',
            'SHARE_65': 'lightblue',
            'SHARE_NYW': 'lightgreen'
        }
        # 绘制NyuWa图（上部）
        bottom = np.zeros(len(sv_types))
        parts_nw = ['SHARE_BOTH', 'SHARE_1019', 'SHARE_65', 'UNIQUE']
        labels_nw = ['Share with both', 'Share with 1019', 'Share with 65', 'Unique']
        for i, part in enumerate(parts_nw):
            if part == 'SHARE_BOTH':
                values = nyuwa_df['BOTH'].values
            elif part == 'SHARE_1019':
                values = nyuwa_df['ONLY_1019'].values
            elif part == 'SHARE_65':
                values = nyuwa_df['ONLY_65'].values
            else:
                values = nyuwa_df['UNIQUE'].values
            ax1.bar(x_adjusted, values, width, bottom=bottom, color=colors[part],
                    label=labels_nw[i], edgecolor='white', linewidth=0.5)
            bottom += values
        ax1.set_ylabel('NyuWa Variant Count', fontsize=12, fontweight='bold')
        ax1.set_ylim(0, y_max)
        ax1.set_title('NyuWa SV Type Distribution - Overlap with 65 and 1019', fontsize=14, fontweight='bold')
        ax1.set_xticks(x_adjusted)
        ax1.set_xticklabels(sv_types, rotation=45, ha='right', fontsize=10)
        left_margin = width * 0.8
        right_margin = width * 0.8
        if len(x_adjusted) > 1:
            ax1.set_xlim(x_adjusted[0] - left_margin, x_adjusted[-1] + right_margin)
        else:
            ax1.set_xlim(x_adjusted[0] - width, x_adjusted[0] + width)
        ax1.legend(loc='upper right', fontsize=10)
        for i, total in enumerate(nyuwa_df['TOTAL_COUNT']):
            if total > 0:
                ax1.text(x_adjusted[i], total + y_max * 0.01, f'{int(total):,}',
                         ha='center', va='bottom', fontsize=9, fontweight='bold')
        # 添加Recall曲线（右轴）
        ax1_right = ax1.twinx()
        recall_values_nw = [nw_recall.get(sv, 0.0) for sv in sv_types]
        ax1_right.plot(x_adjusted, recall_values_nw, 'r-o', linewidth=2, markersize=6, label='Recall (65 dataset)')
        ax1_right.set_ylabel('Recall', fontsize=12, fontweight='bold', color='red')
        ax1_right.set_ylim(0, 1.05)
        ax1_right.tick_params(axis='y', labelcolor='red')
        ax1_right.legend(loc='upper left', fontsize=10)

        # 绘制1019图（下部）
        bottom = np.zeros(len(sv_types))
        parts_kg = ['SHARE_BOTH', 'SHARE_NYW', 'SHARE_65', 'UNIQUE']
        labels_kg = ['Share with both', 'Share with NyuWa', 'Share with 65', 'Unique']
        for i, part in enumerate(parts_kg):
            if part == 'SHARE_BOTH':
                values = kg_df['BOTH'].values
            elif part == 'SHARE_NYW':
                values = kg_df['ONLY_NYW'].values
            elif part == 'SHARE_65':
                values = kg_df['ONLY_65'].values
            else:
                values = kg_df['UNIQUE'].values
            ax2.bar(x_adjusted, values, width, bottom=bottom, color=colors[part],
                    label=labels_kg[i], edgecolor='white', linewidth=0.5)
            bottom += values
        ax2.set_ylabel('1019 Variant Count', fontsize=12, fontweight='bold')
        ax2.set_ylim(0, y_max)
        ax2.set_title('1019 SV Type Distribution - Overlap with 65 and NyuWa', fontsize=14, fontweight='bold')
        ax2.set_xticks(x_adjusted)
        ax2.set_xticklabels(sv_types, rotation=45, ha='right', fontsize=10)
        if len(x_adjusted) > 1:
            ax2.set_xlim(x_adjusted[0] - left_margin, x_adjusted[-1] + right_margin)
        else:
            ax2.set_xlim(x_adjusted[0] - width, x_adjusted[0] + width)
        ax2.legend(loc='upper right', fontsize=10)
        for i, total in enumerate(kg_df['TOTAL_COUNT']):
            if total > 0:
                ax2.text(x_adjusted[i], total + y_max * 0.01, f'{int(total):,}',
                         ha='center', va='bottom', fontsize=9, fontweight='bold')
        # 添加Recall曲线（右轴）
        ax2_right = ax2.twinx()
        recall_values_kg = [kg_recall.get(sv, 0.0) for sv in sv_types]
        ax2_right.plot(x_adjusted, recall_values_kg, 'r-o', linewidth=2, markersize=6, label='Recall (65 dataset)')
        ax2_right.set_ylabel('Recall', fontsize=12, fontweight='bold', color='red')
        ax2_right.set_ylim(0, 1.05)
        ax2_right.tick_params(axis='y', labelcolor='red')
        ax2_right.legend(loc='upper left', fontsize=10)

        ax1.grid(True, alpha=0.3, axis='y')
        ax2.grid(True, alpha=0.3, axis='y')
        plt.tight_layout(pad=1.5, h_pad=2.5, w_pad=1.5)
        plt.savefig(output_file, bbox_inches='tight', facecolor='white', format='pdf')
        print(f"堆叠柱状图（含Recall）已保存到: {output_file}")
    except Exception as e:
        print(f"错误：创建堆叠柱状图失败 - {e}")


def main():
    # 文件路径配置
    kg_all_file = "/Parastor300s_G30S/haoxp/third/cmp/65/1019.svall.inv.addsvtype.rmsvtr.notduptoins_svlenge50.rm05mesv.rm05trsv_newsvtype.txt"
    nw_all_file = "/Parastor300s_G30S/haoxp/third/cmp/65/nyuwa.svall.addsvtype.svlenge50.rm05mesv.rm05trsv_newsvtype.txt"
    cmp_1019_nw_file = "/Parastor300s_G30S/haoxp/third/cmp/65/svcmp.1019.nw.sv.inv.over05.rmnan.addsvtype.rmdismatchINV.notduptoins_svlenge50.rm05mesv.rm05trsv.newsvtype.txt"
    cmp_65_nw_file = "/Parastor300s_G30S/haoxp/third/cmp/65/svcmp.65.nw.sv.inv.over05.rmnan.addsvtype.notduptoins_svlenge50.rm05mesv.rm05trsv.newsvtype.txt"
    cmp_65_1019_file = "/Parastor300s_G30S/haoxp/third/cmp/65/svcmp.65.1019.sv.inv.over05.rmnan.addsvtype.rmdismatchINV.notduptoins_svlenge50.rm05mesv.rm05trsv.newsvtype.txt"
    recall_file = "/Parastor300s_G30S/haoxp/third/cmp/65/16.sv.me.inv.info.txt"

    output_stat_file = "/Parastor300s_G30S/haoxp/third/cmp/65/65allsv_1019sv_1019inv_nwsv_nwme_allsvtype_detail_over05.notduptoins_svlenge50.rm05mesv.rm05trsv.recall.txt"
    output_plot_file = "/Parastor300s_G30S/haoxp/third/cmp/65/65allsv_1019sv_1019inv_nwsv_nwme_allsvtype_distribution_over05.notduptoins_svlenge50.rm05mesv.rm05trsv.recall.pdf"

    print("开始分析多数据集SV变异类型重叠关系（含Recall）...")

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
    # 4. 计算Recall比例
    sv_types = ['DEL', 'INS', 'DUP', 'INV', 'Alu', 'LINE', 'LTR', 'SVA']
    nw_recall, kg_recall = compute_recall_ratios(recall_file, cmp_65_nw_file, cmp_65_1019_file, sv_types)
    # 5. 保存详细统计结果（含Recall）
    save_detailed_statistics_with_recall(nyuwa_df, kg_df, nw_recall, kg_recall, output_stat_file)
    # 6. 创建堆叠柱状图（含Recall曲线）
    create_stacked_bar_plot_with_recall(nyuwa_df, kg_df, nw_recall, kg_recall, output_plot_file)
    print("处理完成！")


if __name__ == "__main__":
    main()
