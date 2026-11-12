#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from pathlib import Path
import gzip
import re

# ==================== 配置路径 ====================
BASE_DIR = Path("/Parastor300s_G30S/haoxp/third/cmp/25")
SVINFO_1KTR = BASE_DIR / "svinfo.1ktr.nomac.addregion.txt"
SVINFO_NWTR = BASE_DIR / "svinfo.nwtr.nomac.addregion.txt"

VCF_DIR = BASE_DIR / "extract/len_vcf"

# 输入文件（只保留与25tr比较的）
CMP_1KTR_DENOM_25TR = BASE_DIR / "svcmp.1ktr.25tr.nomac.rmnoreflen.addspear.over05.txt"
CMP_1KTR_NUM_25TR   = BASE_DIR / "svcmp.1ktr.25tr.nomac.rmnoreflen.addspear.over05.or05.txt"

CMP_NWTR_DENOM_25TR = BASE_DIR / "svcmp.25tr.nwtr.nomac.rmnoreflen.addspear.over05.txt"
CMP_NWTR_NUM_25TR   = BASE_DIR / "svcmp.25tr.nwtr.nomac.rmnoreflen.addspear.over05.or05.txt"

# 输出过滤后的文件（加入 .rmonedata）
OUT_1KTR_DENOM_25TR = BASE_DIR / "svcmp.1ktr.25tr.nomac.rmnoreflen.addspear.over05.1kmacmotif.rmonedata.txt"
OUT_1KTR_NUM_25TR   = BASE_DIR / "svcmp.1ktr.25tr.nomac.rmnoreflen.addspear.over05.or05.1kmacmotif.rmonedata.txt"

OUT_NWTR_DENOM_25TR = BASE_DIR / "svcmp.25tr.nwtr.nomac.rmnoreflen.addspear.over05.nwmacmotif.rmonedata.txt"
OUT_NWTR_NUM_25TR   = BASE_DIR / "svcmp.25tr.nwtr.nomac.rmnoreflen.addspear.over05.or05.nwmacmotif.rmonedata.txt"

# 最终绘图和统计文件
OUTPUT_PLOT = BASE_DIR / "1ktr_nwtr.motif.precision.base05.or05.macmotif.rmnoreflen.rmonedata.pdf"
OUTPUT_STAT = BASE_DIR / "1ktr_nwtr.motif.precision.base05.or05.macmotif.rmnoreflen.rmonedata.txt"

REGIONS = ["SD", "SR", "RM", "NR"]
MOTIF_BINS = [
    (2, 2, "2"),
    (3, 3, "3"),
    (4, 4, "4"),
    (5, 5, "5"),
    (6, 6, "6"),
    (7, 25, "7-25"),
    (26, float('inf'), ">25")
]
MOTIF_GROUPS = [label for (low, high, label) in MOTIF_BINS]
CHROMOSOMES = [f"{i}" for i in range(1, 23)]

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False

def get_motif_group(motif_len):
    for low, high, label in MOTIF_BINS:
        if low <= motif_len <= high:
            return label
    return ">25"

def load_svinfo_map(filepath, dataset_name):
    """从 svinfo 提取 (五列键) -> (motif_group, region, motif_len)，排除 motif=1，去重"""
    print(f"读取 svinfo ({dataset_name}): {filepath}")
    mapping = {}
    seen = set()
    with open(filepath, 'r') as f:
        header = f.readline().strip()
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 6:
                continue
            chrom, pos, vid, reflen_str, motif_str, region = parts[:6]
            try:
                motif = int(motif_str)
            except:
                continue
            if motif < 2:
                continue
            key = (chrom, str(pos), vid, str(reflen_str), motif_str)
            if key in seen:
                continue
            seen.add(key)
            motif_group = get_motif_group(motif)
            mapping[key] = (motif_group, region, motif)
    print(f"  加载映射数: {len(mapping):,}")
    return mapping

def build_valid_map_from_vcf(svinfo_map, vcf_dir, file_pattern):
    """构建有效标志：在原始VCF中存在样本绝对值 ≥ motif_len"""
    key3_to_full_and_motif = {}
    for full_key, (_, _, motif_len) in svinfo_map.items():
        chrom, pos, vid, reflen, motif_str = full_key
        key3 = (chrom, pos, vid)
        if key3 not in key3_to_full_and_motif:
            key3_to_full_and_motif[key3] = (full_key, motif_len)

    valid = {}
    for chrom in CHROMOSOMES:
        vcf_path = vcf_dir / file_pattern.format(chrom)
        if not vcf_path.exists():
            print(f"  跳过: {vcf_path}")
            continue
        print(f"  扫描 {vcf_path.name} ...")
        with gzip.open(vcf_path, 'rt') as f:
            for line in f:
                if line.startswith('#CHROM'):
                    break
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 10:
                    continue
                chrom_vcf, pos, vid = fields[0], fields[1], fields[2]
                key3 = (chrom_vcf, pos, vid)
                if key3 not in key3_to_full_and_motif:
                    continue
                full_key, motif_len = key3_to_full_and_motif[key3]

                sample_vals = []
                for val in fields[9:]:
                    val_str = val.split(':')[0]
                    if val_str in ('.', '-', './.', '.|.'):
                        sample_vals.append(None)
                    else:
                        try:
                            sample_vals.append(float(val_str))
                        except ValueError:
                            sample_vals.append(None)
                is_valid = any(v is not None and abs(v) >= motif_len for v in sample_vals)
                valid[full_key] = is_valid
    print(f"  满足条件变异数: {sum(valid.values())} / {len(valid)}")
    return valid

def load_cmp_and_write_filtered(filepath, side, valid_set, sv_map, outpath):
    """
    读取文件，提取key、最大相关性，写入过滤后的文件
    side: 'left' 或 'right' 指定哪5列作为键
    valid_set: 有效的五列key集合
    返回：key集合, key->最大相关性
    """
    key_to_max_corr = defaultdict(float)
    valid_keys = set()
    valid_lines = []

    with open(filepath, 'r') as f:
        # 保留第一行如果是注释行
        first_line = f.readline().rstrip('\n')
        if first_line.startswith('#'):
            valid_lines.append(first_line)
            # 继续读剩余行
        else:
            # 否则将第一行作为数据行处理，需要重置文件指针
            f.seek(0)

        # 重新遍历所有行（包括可能的表头行）
        f.seek(0)
        for line in f:
            raw_line = line.rstrip('\n')
            if not raw_line:
                continue
            # 保留注释行
            if raw_line.startswith('#'):
                valid_lines.append(raw_line)
                continue
            parts = raw_line.split('\t')
            if len(parts) < 13:
                continue
            # 跳过可能的列名行（包含chrom/pos等）
            if any(c in parts[0].lower() for c in ['chrom', 'bed']):
                continue

            if side == 'left':
                key_cols = parts[:5]
            else:
                key_cols = parts[5:10]
            key = tuple(str(x) for x in key_cols)

            if key not in valid_set or key not in sv_map:
                continue

            # 提取相关性最大值（列索引11和12）
            try:
                c1 = float(parts[11])
            except:
                c1 = -float('inf')
            try:
                c2 = float(parts[12])
            except:
                c2 = -float('inf')
            curr = max(c1, c2)
            if curr == -float('inf'):
                curr = 0.0   # 无有效相关性值，贡献0

            # 更新该键的最大相关性
            if curr > key_to_max_corr[key]:
                key_to_max_corr[key] = curr
            valid_keys.add(key)
            valid_lines.append(raw_line)

    # 写入输出文件
    with open(outpath, 'w') as f:
        for line in valid_lines:
            f.write(line + '\n')
    print(f"  输出 {len([l for l in valid_lines if l and not l.startswith('#')])} 行 → {outpath.name}")
    return valid_keys, key_to_max_corr

def main():
    # 1. 加载 svinfo
    map_1ktr = load_svinfo_map(SVINFO_1KTR, "1ktr")
    map_nwtr = load_svinfo_map(SVINFO_NWTR, "nwtr")

    # 2. 构建有效标志
    print("\n===== 构建1ktr有效标志 =====")
    valid_1ktr = build_valid_map_from_vcf(map_1ktr, VCF_DIR, "chr{}.1ktr.len.vcf.gz")
    print("\n===== 构建nwtr有效标志 =====")
    valid_nwtr = build_valid_map_from_vcf(map_nwtr, VCF_DIR, "chr{}.nwtr.nomac.len.vcf.gz")
    valid1_set = set(valid_1ktr.keys())
    validn_set = set(valid_nwtr.keys())

    # 3. 读取并输出4个过滤文件（只保留与25tr比较的）
    print("\n===== 处理1ktr输入文件（仅25tr比较）=====")
    denom1_keys, corr1_denom = load_cmp_and_write_filtered(
        CMP_1KTR_DENOM_25TR, 'left', valid1_set, map_1ktr, OUT_1KTR_DENOM_25TR)
    num1_keys,   corr1_num   = load_cmp_and_write_filtered(
        CMP_1KTR_NUM_25TR,   'left', valid1_set, map_1ktr, OUT_1KTR_NUM_25TR)

    print("\n===== 处理nwtr输入文件（仅25tr比较）=====")
    denom2_keys, corr2_denom = load_cmp_and_write_filtered(
        CMP_NWTR_DENOM_25TR, 'right', validn_set, map_nwtr, OUT_NWTR_DENOM_25TR)
    num2_keys,   corr2_num   = load_cmp_and_write_filtered(
        CMP_NWTR_NUM_25TR,   'right', validn_set, map_nwtr, OUT_NWTR_NUM_25TR)

    # 4. 分母：去重后的变异集合（只使用 denom 文件中的键）
    denom_1ktr = denom1_keys   # 直接使用 denom 文件提供的变异集合
    denom_nwtr = denom2_keys

    # 5. 分子：每个变异取最大相关性（来自 num 文件）
    #    注意：num 文件和 denom 文件可能包含不同键，但实际计算时分母和分子应一致
    #    这里分子只从 num 文件中获取；如果某个变异只出现在 denom 中但不在 num 中，则其相关性视为0
    all_corr_1ktr = {}
    for k, v in corr1_num.items():
        all_corr_1ktr[k] = v
    # 对于在 denom 中但不在 num 中的键，相关性默认为0，已在 compute_weighted_stats 中处理
    all_corr_nwtr = {}
    for k, v in corr2_num.items():
        all_corr_nwtr[k] = v

    # 6. 计算加权统计
    def compute_weighted_stats(sv_map, denom_keys, corr_dict, valid_map):
        denom_count = defaultdict(int)
        sum_corr = defaultdict(float)
        for key in denom_keys:
            if key not in sv_map or key not in valid_map or not valid_map[key]:
                continue
            g, r = sv_map[key][:2]
            denom_count[(g, r)] += 1
        for key, corr in corr_dict.items():
            if key not in sv_map or key not in valid_map or not valid_map[key]:
                continue
            g, r = sv_map[key][:2]
            sum_corr[(g, r)] += corr
        result = {}
        all_keys = set(denom_count.keys()) | set(sum_corr.keys())
        for k in all_keys:
            d = denom_count.get(k, 0)
            s = sum_corr.get(k, 0.0)
            val = s / d if d > 0 else 0.0
            result[k] = (d, s, val)
        return result

    stats_1ktr = compute_weighted_stats(map_1ktr, denom_1ktr, all_corr_1ktr, valid_1ktr)
    stats_nwtr = compute_weighted_stats(map_nwtr, denom_nwtr, all_corr_nwtr, valid_nwtr)

    total_denom1 = sum(d for d,_,_ in stats_1ktr.values())
    total_sum1   = sum(s for _,s,_ in stats_1ktr.values())
    total_denom2 = sum(d for d,_,_ in stats_nwtr.values())
    total_sum2   = sum(s for _,s,_ in stats_nwtr.values())
    print(f"\n1ktr 分母总数（去重）: {total_denom1}, 分子总和: {total_sum1:.6f}")
    print(f"nwtr 分母总数（去重）: {total_denom2}, 分子总和: {total_sum2:.6f}")

    # 7. 准备绘图数据
    x_pos = np.arange(len(MOTIF_GROUPS) * len(REGIONS))
    y1 = []
    y2 = []
    stat_out1 = []
    stat_out2 = []

    for mg in MOTIF_GROUPS:
        for r in REGIONS:
            k = (mg, r)
            d1, s1, v1 = stats_1ktr.get(k, (0,0,0))
            d2, s2, v2 = stats_nwtr.get(k, (0,0,0))
            y1.append(v1)
            if r == "SD" or d2 < 5:
                y2.append(np.nan)
            else:
                y2.append(v2)
            stat_out1.append((mg, r, d1, s1, v1))
            stat_out2.append((mg, r, d2, s2, v2, (r=="SD" or d2<5)))

    # 8. 绘图
    plt.figure(figsize=(14,6))
    plt.plot(x_pos, y1, 'b-o', lw=2, ms=6, label='1ktr')
    plt.plot(x_pos, y2, 'r-o', lw=2, ms=6, label='nwtr')
    plt.xticks(x_pos, [f"{mg}\n{r}" for mg in MOTIF_GROUPS for r in REGIONS], rotation=45, ha='right', fontsize=9)
    plt.ylim(0, 1.05)
    plt.grid(alpha=0.3, ls='--')
    plt.ylabel('Weighted Correlation Mean', fontsize=12)
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, bbox_inches='tight', facecolor='w')
    print(f"绘图保存: {OUTPUT_PLOT}")

    # 9. 输出统计文件
    with open(OUTPUT_STAT, 'w') as f:
        f.write("# motif_group\tregion\tdenominator\tmax_corr_sum\tweighted_mean\n")
        f.write("# 1ktr results (only vs 25tr)\n")
        for mg, r, d, s, v in stat_out1:
            f.write(f"{mg}\t{r}\t{d}\t{s:.6f}\t{v:.6f}\n")
        f.write("\n# nwtr results (only vs 25tr, kept_in_plot=no 表示该点在图中被隐藏)\n")
        for mg, r, d, s, v, kept in stat_out2:
            kept_str = "no" if kept else "yes"
            f.write(f"{mg}\t{r}\t{d}\t{s:.6f}\t{v:.6f}\t{kept_str}\n")
        f.write(f"\n# 1ktr 分母总键数（去重）: {total_denom1}\n")
        f.write(f"# 1ktr 分子加权总和: {total_sum1:.6f}\n")
        f.write(f"# nwtr 分母总键数（去重）: {total_denom2}\n")
        f.write(f"# nwtr 分子加权总和: {total_sum2:.6f}\n")
    print(f"统计文件保存: {OUTPUT_STAT}")

    # 10. 控制台输出分组分布
    print("\n各 motif 分组 1ktr 变异总数（有效统计）：")
    for mg in MOTIF_GROUPS:
        for r in REGIONS:
            d = stats_1ktr.get((mg, r), (0,0,0))[0]
            if d > 0:
                print(f"  {mg} {r}: {d}")

if __name__ == "__main__":
    main()
