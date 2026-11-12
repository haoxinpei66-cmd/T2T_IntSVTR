#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import time

# ==================== 绘图兼容设置 ====================
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']

plt.rcParams.update({
    'font.size': 12,
    'axes.titlesize': 13,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
    'figure.dpi': 300,
})

# ==================== 核心功能 ====================
def get_reflen_category(reflen: int) -> str:
    reflen_int = int(reflen)
    if reflen_int <= 10:
        return "1-10"
    elif 11 <= reflen_int <= 20:
        return "11-20"
    elif 21 <= reflen_int <= 30:
        return "21-30"
    elif 31 <= reflen_int <= 40:
        return "31-40"
    elif 41 <= reflen_int <= 60:
        return "41-60"
    elif 61 <= reflen_int <= 120:
        return "61-120"
    elif 121 <= reflen_int <= 200:
        return "121-200"
    else:
        return ">200"

# -----------------------------------------------------------------------------
# 加载注释文件：生成两套 key
# 严格key（子集）：CHROM, POS, VARID, REFLEN, MOTIF_LEN
# 宽松key（总集）：CHROM, POS, REFLEN, MOTIF_LEN
# -----------------------------------------------------------------------------
def load_anno_dict(anno_file):
    print(f"Loading annotation: {anno_file}")
    df = pd.read_csv(anno_file, sep='\t', comment='#', header=0, low_memory=False)
    df.columns = ['CHROM', 'POS', 'VARID', 'REFLEN', 'MOTIF_LEN', 'REGION']
    df['REGION'] = df['REGION'].replace('NR', 'Unique')

    strict_dict = {}
    loose_dict = {}
    for _, r in df.iterrows():
        strict_key = (str(r['CHROM']), int(r['POS']), str(r['VARID']), int(r['REFLEN']), int(r['MOTIF_LEN']))
        loose_key = (str(r['CHROM']), int(r['POS']), int(r['REFLEN']), int(r['MOTIF_LEN']))
        strict_dict[strict_key] = r['REGION']
        loose_dict[loose_key] = r['REGION']
    return strict_dict, loose_dict

# -----------------------------------------------------------------------------
# 读取【标准格式】文件（总集：CHROM, START, REFLEN, MOTIF_LEN, ...）
# -----------------------------------------------------------------------------
def read_standard_file(f):
    df = pd.read_csv(f, sep='\t', header=None, names=['CHROM', 'START', 'REFLEN', 'MOTIF_LEN', 'MOTIF_SCORE'], low_memory=False)
    chroms = [f'chr{i}' for i in range(1,23)]
    df = df[df['CHROM'].isin(chroms)].copy()
    df['POS'] = df['START']
    df['VARID'] = 'N/A'
    return df

# -----------------------------------------------------------------------------
# 读取【新子集格式】文件（CHROM, POS, VARID, REFLEN, MOTIF_LEN）
# -----------------------------------------------------------------------------
def read_new_subset_file(f):
    df = pd.read_csv(f, sep='\t', header=None, names=['CHROM', 'POS', 'VARID', 'REFLEN', 'MOTIF_LEN'], low_memory=False)
    chroms = [f'chr{i}' for i in range(1,23)]
    df = df[df['CHROM'].isin(chroms)].copy()
    return df

# -----------------------------------------------------------------------------
# 总集注释：宽松匹配（忽略 VARID）
# -----------------------------------------------------------------------------
def annotate_total_loose(df, loose_dict):
    regs = []
    for _, r in df.iterrows():
        key = (str(r['CHROM']), int(r['POS']), int(r['REFLEN']), int(r['MOTIF_LEN']))
        regs.append(loose_dict.get(key, 'Unique'))
    df['REGION'] = regs
    df['REFLEN_CATEGORY'] = df['REFLEN'].apply(get_reflen_category)
    return df

# -----------------------------------------------------------------------------
# 子集注释：严格匹配（带 VARID）
# -----------------------------------------------------------------------------
def annotate_subset_strict(df, strict_dict):
    regs = []
    for _, r in df.iterrows():
        key = (str(r['CHROM']), int(r['POS']), str(r['VARID']), int(r['REFLEN']), int(r['MOTIF_LEN']))
        regs.append(strict_dict.get(key, 'Unique'))
    df['REGION'] = regs
    df['REFLEN_CATEGORY'] = df['REFLEN'].apply(get_reflen_category)
    return df

# -----------------------------------------------------------------------------
# 统计表格
# -----------------------------------------------------------------------------
def create_statistics(kg_total_df, kg_subset_df, nw_total_df, nw_subset_df):
    reflen_cats = ["1-10","11-20","21-30","31-40","41-60","61-120","121-200",">200"]
    regions = ["SD","SR","RM","Unique"]
    rows = []
    for rl in reflen_cats:
        for reg in regions:
            kt = len(kg_total_df[(kg_total_df['REFLEN_CATEGORY']==rl) & (kg_total_df['REGION']==reg)])
            ks = len(kg_subset_df[(kg_subset_df['REFLEN_CATEGORY']==rl) & (kg_subset_df['REGION']==reg)])
            nt = len(nw_total_df[(nw_total_df['REFLEN_CATEGORY']==rl) & (nw_total_df['REGION']==reg)])
            ns = len(nw_subset_df[(nw_subset_df['REFLEN_CATEGORY']==rl) & (nw_subset_df['REGION']==reg)])
            kr = ks/kt if kt>0 else 0
            nr = ns/nt if nt>0 else 0
            rows.append([rl, reg, kt, ks, kr, nt, ns, nr])
    return pd.DataFrame(rows, columns=["REFLEN_CATEGORY","REGION","1KG_TOTAL","1KG_SUBSET","1KG_RATIO","NyuWa_TOTAL","NyuWa_SUBSET","NyuWa_RATIO"])

# -----------------------------------------------------------------------------
# 保存结果
# -----------------------------------------------------------------------------
def save_stats(df, out):
    df.to_csv(out, sep='\t', index=False)
    print(f"Stats saved: {out}")

# -----------------------------------------------------------------------------
# 绘图
# -----------------------------------------------------------------------------
def plot_figure(stats, out_plot):
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    regions = ["SD","SR","RM","Unique"]
    reflens = ["1-10","11-20","21-30","31-40","41-60","61-120","121-200",">200"]
    fig, (ax1, ax2) = plt.subplots(2,1,figsize=(16,10))
    x = np.arange(len(reflens))
    w = 0.18
    sp = 0.03

    for i, reg in enumerate(regions):
        d = stats[stats['REGION']==reg].set_index('REFLEN_CATEGORY')
        ax1.bar(x+i*(w+sp), [1]*len(x), w, color=colors[i], alpha=0.6)
        ax1.bar(x+i*(w+sp), d.loc[reflens]['NyuWa_RATIO'].values, w, color=colors[i], edgecolor='k', linewidth=0.5)
        ax2.bar(x+i*(w+sp), [1]*len(x), w, color=colors[i], alpha=0.6)
        ax2.bar(x+i*(w+sp), d.loc[reflens]['1KG_RATIO'].values, w, color=colors[i], edgecolor='k', linewidth=0.5)

    ax1.set_title('NyuWa')
    ax2.set_title('1KG')
    for ax in [ax1, ax2]:
        ax.set_xticks(x + 1.5*(w+sp))
        ax.set_xticklabels(reflens, rotation=45, ha='right')
        ax.set_ylim(0,1.1)
        ax.grid(alpha=0.3)
    ax1.legend(regions, loc='upper left', bbox_to_anchor=(1,1))
    plt.tight_layout()
    plt.savefig(out_plot, bbox_inches='tight', format='pdf')
    print(f"Plot saved: {out_plot}")

# ==================== 主函数 ====================
def main():
    # ========== 路径配置（完全按你要求） ==========
    kg_anno = "/Parastor300s_G30S/haoxp/third/cmp/1019/svcmp.nwtr_1ktr/svinfo.1ktr.addregion.txt"
    nw_anno = "/Parastor300s_G30S/haoxp/third/cmp/1019/svcmp.nwtr_1ktr/svinfo.nwtr.addregion.txt"

    kg_total_files = [
        "/Parastor300s_G30S/haoxp/third/cmp/1019/extract/reflen_vntr_rmXY_rmdup_addmotif.txt",
        "/Parastor300s_G30S/haoxp/third/cmp/1019/extract/reflen_svtr_rmover05vntr_realreflen_addmotif.txt"
    ]
    kg_subset_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/reflen_1ktr_rmXY_rmdup_addmotif.txt"
    nw_total_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/extract/reflen_nwtr_filt_SD_addmotif.txt"
    nw_subset_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/reflen_nwtr_rmXY_rmdup_addmotif.txt"

    output_stat_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/RU_motif_ALTANNO/1019svtr_1019vntr_nwtr_reflen_region_ratio_detail_over05_adjust_new.txt"
    output_plot_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/RU_motif_ALTANNO/1019svtr_1019vntr_nwtr_reflen_region_ratio_distribution_over05_adjust_new.pdf"

    print("===== 开始运行（总集宽松，子集严格）=====")
    start = time.time()

    # 1. 加载注释
    kg_strict, kg_loose = load_anno_dict(kg_anno)
    nw_strict, nw_loose = load_anno_dict(nw_anno)

    # 2. 1KG 总集
    kg_total = pd.concat([read_standard_file(f) for f in kg_total_files], ignore_index=True)
    kg_total = annotate_total_loose(kg_total, kg_loose)

    # 3. 1KG 子集（新格式）
    kg_sub = read_new_subset_file(kg_subset_file)
    kg_sub = annotate_subset_strict(kg_sub, kg_strict)

    # 4. NyuWa 总集
    nw_total = read_standard_file(nw_total_file)
    nw_total = annotate_total_loose(nw_total, nw_loose)

    # 5. NyuWa 子集（新格式）
    nw_sub = read_new_subset_file(nw_subset_file)
    nw_sub = annotate_subset_strict(nw_sub, nw_strict)

    # 6. 统计 + 输出
    stats = create_statistics(kg_total, kg_sub, nw_total, nw_sub)
    save_stats(stats, output_stat_file)
    plot_figure(stats, output_plot_file)

    print(f"全部完成！耗时：{(time.time()-start)/60:.1f} min")

if __name__ == "__main__":
    main()
