#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, Tuple, List, Set
import time
from tqdm import tqdm

# ==================== 绘图全局设置 ====================
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

def get_motif_category(motif_len: int) -> str:
    if motif_len == 1:
        return "1"
    elif motif_len == 2:
        return "2"
    elif motif_len == 3:
        return "3"
    elif motif_len == 4:
        return "4"
    elif motif_len == 5:
        return "5"
    elif motif_len == 6:
        return "6"
    elif 7 <= motif_len <= 25:
        return "7-25"
    elif motif_len > 25:
        return ">25"
    else:
        return "Unknown"

# ==================== 加载注释文件：返回两套字典 ====================
def load_annotated_region_file(file_path: str, dataset_name: str):
    print(f"正在加载 {dataset_name} 区域注释文件: {file_path}")
    strict_key_dict = {}    # 子集用：带VARID
    loose_key_dict = {}     # 总集用：忽略VARID
    
    df = pd.read_csv(
        file_path, sep='\t', comment='#', header=0,
        dtype={'CHROM': str, 'POS': int, 'VARID': str, 'REFLEN': int, 'MOTIF_LEN': int, 'REGION': str},
        low_memory=False
    )
    df.columns = ['CHROM', 'POS', 'VARID', 'REFLEN', 'MOTIF_LEN', 'REGION']
    
    for _, row in df.iterrows():
        reg = row['REGION'] if row['REGION'] != 'NR' else 'Unique'
        # 严格key：子集匹配
        strict_key = (str(row['CHROM']), int(row['POS']), str(row['VARID']), 
                     int(row['REFLEN']), int(row['MOTIF_LEN']))
        strict_key_dict[strict_key] = reg
        
        # 宽松key：总集匹配，忽略VARID
        loose_key = (str(row['CHROM']), int(row['POS']), 
                     int(row['REFLEN']), int(row['MOTIF_LEN']))
        loose_key_dict[loose_key] = reg
    
    print(f"加载完成：{len(strict_key_dict)} 条 {dataset_name} 注释")
    return strict_key_dict, loose_key_dict

# ==================== 读取总集文件 ====================
def read_tr_file_total(file_path: str, dataset_name: str):
    print(f"正在读取 {dataset_name} TR总集文件: {file_path}")
    try:
        df = pd.read_csv(file_path, sep='\t', header=None,
                        names=['CHROM', 'START', 'REFLEN', 'MOTIF_LEN', 'MOTIF_SCORE'],
                        low_memory=False)
        chroms = [f'chr{i}' for i in range(1, 23)]
        chr_df = df[df['CHROM'].isin(chroms)].copy()
        print(f"{dataset_name} 总集过滤后：{len(chr_df)} 条（1-22号染色体）")
        return chr_df
    except Exception as e:
        print(f"读取失败：{e}")
        return pd.DataFrame()

# ==================== 读取子集文件 ====================
def read_tr_file_subset_new(file_path: str, dataset_name: str):
    print(f"正在读取 {dataset_name} 子集文件: {file_path}")
    try:
        df = pd.read_csv(file_path, sep='\t', header=None,
                        names=['CHROM', 'POS', 'VARID', 'REFLEN', 'MOTIF_LEN'],
                        dtype={'CHROM': str, 'POS': int, 'VARID': str, 'REFLEN': int, 'MOTIF_LEN': int},
                        low_memory=False)
        chroms = [f'chr{i}' for i in range(1, 23)]
        chr_df = df[df['CHROM'].isin(chroms)].copy()
        print(f"{dataset_name} 子集过滤后：{len(chr_df)} 条（1-22号染色体）")
        return chr_df
    except Exception as e:
        print(f"读取失败：{e}")
        return pd.DataFrame()

# ==================== 总集注释：用宽松key，忽略VARID ====================
def annotate_total_with_loose_key(df: pd.DataFrame, loose_key_dict: dict, desc: str):
    print(f"为 {desc} 匹配区域（忽略VARID宽松匹配）...")
    regions = []
    for _, row in tqdm(df.iterrows(), total=len(df), desc=desc):
        lkey = (str(row['CHROM']), int(row['START']),
                int(row['REFLEN']), int(row['MOTIF_LEN']))
        region = loose_key_dict.get(lkey, 'Unique')
        regions.append(region)
    df['REGION'] = regions
    df['MOTIF_CATEGORY'] = df['MOTIF_LEN'].apply(get_motif_category)
    return df

# ==================== 子集注释：用严格key，带VARID ====================
def annotate_subset_with_strict_key(df: pd.DataFrame, strict_key_dict: dict, desc: str):
    print(f"为 {desc} 匹配区域（带VARID严格匹配）...")
    regions = []
    for _, row in tqdm(df.iterrows(), total=len(df), desc=desc):
        skey = (str(row['CHROM']), int(row['POS']), str(row['VARID']),
               int(row['REFLEN']), int(row['MOTIF_LEN']))
        region = strict_key_dict.get(skey, 'Unique')
        regions.append(region)
    df['REGION'] = regions
    df['MOTIF_CATEGORY'] = df['MOTIF_LEN'].apply(get_motif_category)
    return df

# ==================== 统计：按MOTIF_CATEGORY + REGION 独立计数 ====================
def create_statistics(kg_total_df: pd.DataFrame, kg_subset_df: pd.DataFrame,
                     nw_total_df: pd.DataFrame, nw_subset_df: pd.DataFrame):
    print("开始按motif+区域独立统计数量...")
    motif_categories = ["1", "2", "3", "4", "5", "6", "7-25", ">25"]
    region_types = ["SD", "SR", "RM", "Unique"]
    results = []

    for motif_cat in motif_categories:
        for region in region_types:
            # 每个组合单独统计
            kg_total_cnt = len(kg_total_df[
                (kg_total_df['MOTIF_CATEGORY'] == motif_cat) &
                (kg_total_df['REGION'] == region)
            ])
            kg_subset_cnt = len(kg_subset_df[
                (kg_subset_df['MOTIF_CATEGORY'] == motif_cat) &
                (kg_subset_df['REGION'] == region)
            ])

            nw_total_cnt = len(nw_total_df[
                (nw_total_df['MOTIF_CATEGORY'] == motif_cat) &
                (nw_total_df['REGION'] == region)
            ])
            nw_subset_cnt = len(nw_subset_df[
                (nw_subset_df['MOTIF_CATEGORY'] == motif_cat) &
                (nw_subset_df['REGION'] == region)
            ])

            kg_ratio = kg_subset_cnt / kg_total_cnt if kg_total_cnt > 0 else 0.0
            nw_ratio = nw_subset_cnt / nw_total_cnt if nw_total_cnt > 0 else 0.0

            row = {
                'MOTIF_CATEGORY': motif_cat,
                'REGION': region,
                '1KG_TOTAL': kg_total_cnt,
                '1KG_SUBSET': kg_subset_cnt,
                '1KG_RATIO': kg_ratio,
                'NyuWa_TOTAL': nw_total_cnt,
                'NyuWa_SUBSET': nw_subset_cnt,
                'NyuWa_RATIO': nw_ratio
            }
            results.append(row)
    return pd.DataFrame(results)

def save_detailed_stats(df: pd.DataFrame, output_file: str):
    with open(output_file, 'w') as f:
        f.write("MOTIF_CATEGORY\tREGION\t1KG_TOTAL\t1KG_SUBSET\t1KG_RATIO\tNyuWa_TOTAL\tNyuWa_SUBSET\tNyuWa_RATIO\n")
        for _, row in df.iterrows():
            f.write(f"{row['MOTIF_CATEGORY']}\t{row['REGION']}\t{row['1KG_TOTAL']}\t{row['1KG_SUBSET']}\t{row['1KG_RATIO']:.6f}\t{row['NyuWa_TOTAL']}\t{row['NyuWa_SUBSET']}\t{row['NyuWa_RATIO']:.6f}\n")
    print(f"统计文件已保存：{output_file}")

# ==================== 绘图函数 ====================
def plot_ratio_figure(stats_df, output_plot_file):
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Liberation Sans']
    plt.rcParams['axes.unicode_minus'] = False
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    region_types = ['SD', 'SR', 'RM', 'Unique']
    motif_categories = ["1", "2", "3", "4", "5", "6", "7-25", ">25"]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(18, 12))
    x = np.arange(len(motif_categories))
    width = 0.18
    spacing = 0.02

    nw_data = {}
    kg_data = {}
    for region in region_types:
        rd = stats_df[stats_df['REGION'] == region]
        nw_data[region] = [rd[rd['MOTIF_CATEGORY']==m]['NyuWa_RATIO'].values[0] for m in motif_categories]
        kg_data[region] = [rd[rd['MOTIF_CATEGORY']==m]['1KG_RATIO'].values[0] for m in motif_categories]

    # NyuWa
    for i, region in enumerate(region_types):
        pos = x + i*(width+spacing)
        ax1.bar(pos, [1]*len(pos), width, color=colors[i], alpha=0.6)
        ax1.bar(pos, nw_data[region], width, color=colors[i], alpha=1, edgecolor='black', linewidth=0.5)
    ax1.set_ylim(0,1.1)
    ax1.set_ylabel('NyuWa Ratio', fontsize=12, weight='bold')
    ax1.set_title('NyuWa TR Ratio Distribution by Region and Motif Length (1-22)', fontsize=14, weight='bold')
    ax1.set_xticks(x + 1.5*(width+spacing))
    ax1.set_xticklabels(motif_categories, rotation=45, ha='right')
    ax1.grid(True, alpha=0.3, axis='y')

    # 1KG
    for i, region in enumerate(region_types):
        pos = x + i*(width+spacing)
        ax2.bar(pos, [1]*len(pos), width, color=colors[i], alpha=0.6)
        ax2.bar(pos, kg_data[region], width, color=colors[i], alpha=1, edgecolor='black', linewidth=0.5)
    ax2.set_ylim(0,1.1)
    ax2.set_ylabel('1KG Ratio', fontsize=12, weight='bold')
    ax2.set_title('1KG TR Ratio Distribution by Region and Motif Length (1-22)', fontsize=14, weight='bold')
    ax2.set_xticks(x + 1.5*(width+spacing))
    ax2.set_xticklabels(motif_categories, rotation=45, ha='right')
    ax2.grid(True, alpha=0.3, axis='y')

    plt.tight_layout(rect=[0,0,0.9,0.95])
    plt.savefig(output_plot_file, bbox_inches='tight', facecolor='white', format='pdf')
    print(f"绘图已保存：{output_plot_file}")

# ==================== 主入口 ====================
def main():
    # 路径保持你原来的不变
    ktr_region_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/svcmp.nwtr_1ktr/svinfo.1ktr.addregion.txt"
    nwtr_region_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/svcmp.nwtr_1ktr/svinfo.nwtr.addregion.txt"

    kg_total_files = [
        "/Parastor300s_G30S/haoxp/third/cmp/1019/extract/reflen_vntr_rmXY_rmdup_addmotif.txt",
        "/Parastor300s_G30S/haoxp/third/cmp/1019/extract/reflen_svtr_rmover05vntr_realreflen_addmotif.txt"
    ]

    kg_subset_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/reflen_1ktr_rmXY_rmdup_addmotif.txt"
    nw_total_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/extract/reflen_nwtr_filt_SD_addmotif.txt"
    nw_subset_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/reflen_nwtr_rmXY_rmdup_addmotif.txt"

    output_stat_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/RU_motif_ALTANNO/1019svtr_1019vntr_nwtr_motif_region_ratio_detail_over05_adjust_new.txt"
    output_plot_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/RU_motif_ALTANNO/1019svtr_1019vntr_nwtr_motif_region_ratio_distribution_over05_adjust_new.pdf"

    print("===== 开始统计：总集忽略VARID匹配，子集严格带VARID匹配 =====")
    start_time = time.time()

    # 加载注释：两套字典
    kg_strict_dict, kg_loose_dict = load_annotated_region_file(ktr_region_file, "1KG")
    nw_strict_dict, nw_loose_dict = load_annotated_region_file(nwtr_region_file, "NyuWa")

    # 处理1KG总集：宽松匹配
    kg_total_dfs = []
    for f in kg_total_files:
        df = read_tr_file_total(f, "1KG")
        if not df.empty:
            df = annotate_total_with_loose_key(df, kg_loose_dict, "1KG总集")
            kg_total_dfs.append(df)
    kg_total_combined = pd.concat(kg_total_dfs, ignore_index=True)

    # 处理1KG子集：严格匹配
    kg_subset_df = read_tr_file_subset_new(kg_subset_file, "1KG子集")
    kg_subset_classified = annotate_subset_with_strict_key(kg_subset_df, kg_strict_dict, "1KG子集")

    # 处理NyuWa总集：宽松匹配
    nw_total_df = read_tr_file_total(nw_total_file, "NyuWa")
    nw_total_classified = annotate_total_with_loose_key(nw_total_df, nw_loose_dict, "NyuWa总集")

    # 处理NyuWa子集：严格匹配
    nw_subset_df = read_tr_file_subset_new(nw_subset_file, "NyuWa子集")
    nw_subset_classified = annotate_subset_with_strict_key(nw_subset_df, nw_strict_dict, "NyuWa子集")

    # 统计+输出
    stats_df = create_statistics(kg_total_combined, kg_subset_classified,
                                nw_total_classified, nw_subset_classified)
    save_detailed_stats(stats_df, output_stat_file)

    # 绘图
    plot_ratio_figure(stats_df, output_plot_file)

    print(f"全部完成，耗时：{(time.time()-start_time)/60:.2f} 分钟")

if __name__ == "__main__":
    main()
