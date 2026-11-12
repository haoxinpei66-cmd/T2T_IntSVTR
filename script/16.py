#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from pathlib import Path
import gzip

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
OUT_1KTR_DENOM_25TR = BASE_DIR / "svcmp.1ktr.25tr.nomac.rmnoreflen.addspear.over05.1k2macmotif.rmonedata.txt"
OUT_1KTR_NUM_25TR   = BASE_DIR / "svcmp.1ktr.25tr.nomac.rmnoreflen.addspear.over05.or05.1k2macmotif.rmonedata.txt"

OUT_NWTR_DENOM_25TR = BASE_DIR / "svcmp.25tr.nwtr.nomac.rmnoreflen.addspear.over05.nw2macmotif.rmonedata.txt"
OUT_NWTR_NUM_25TR   = BASE_DIR / "svcmp.25tr.nwtr.nomac.rmnoreflen.addspear.over05.or05.nw2macmotif.rmonedata.txt"

OUTPUT_PLOT = BASE_DIR / "1ktr_nwtr.reflen.precision.base05.or05.macmotif.rmnoreflen.rmonedata.pdf"
OUTPUT_STAT = BASE_DIR / "1ktr_nwtr.reflen.precision.base05.or05.macmotif.rmnoreflen.rmonedata.txt"

REGIONS = ["SD", "SR", "RM", "NR"]
REFLEN_BINS = [
    (1, 10, "1-10"),
    (11, 20, "11-20"),
    (21, 30, "21-30"),
    (31, 40, "31-40"),
    (41, 60, "41-60"),
    (61, 120, "61-120"),
    (121, 200, "121-200"),
    (201, float('inf'), ">200")
]
REFLEN_GROUPS = [label for (low, high, label) in REFLEN_BINS]

CHROMOSOMES = [f"{i}" for i in range(1, 23)]

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False

def get_reflen_group(reflen):
    if reflen <= 0:
        return None
    for low, high, label in REFLEN_BINS:
        if low <= reflen <= high:
            return label
    return ">200"

def load_svinfo_map(filepath, dataset_name):
    """加载 svinfo，返回 (五列键) -> (reflen_group, region, motif_len) ，过滤 motif<2, reflen<=0"""
    print(f"读取 {dataset_name}: {filepath}")
    mapping = {}
    seen = set()
    with open(filepath, 'r') as f:
        header = f.readline().strip()
        for line in f:
            line = line.rstrip('\n')
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
            if region not in REGIONS:
                continue
            try:
                reflen = int(reflen_str)
            except:
                reflen = 0
            if reflen <= 0:
                continue
            reflen_group = get_reflen_group(reflen)
            if reflen_group is None:
                continue
            key = (chrom, str(pos), vid, str(reflen_str), motif_str)
            if key in seen:
                continue
            seen.add(key)
            mapping[key] = (reflen_group, region, motif)  # 保存 motif_len 用于后续校验
    print(f"  加载映射数: {len(mapping):,}")
    return mapping

def load_cmp_keys(filepath, position, label):
    """从比较文件提取变异键集合（五列），跳过注释和表头（仅供检查，不用）"""
    # 此函数保留但不使用，因为我们在 process_cmp_file 中同时做了输出和键提取
    pass

def build_valid_map(svinfo_map, file_suffix):
    """
    预加载VCF文件，检查每个变异是否存在样本值 绝对值 >= motif_len。
    返回字典: key -> bool
    """
    key3_to_info = {}
    for full_key, (_, _, motif_len) in svinfo_map.items():
        chrom, pos, vid, reflen, motif_str = full_key
        key3 = (chrom, pos, vid)
        if key3 not in key3_to_info:
            key3_to_info[key3] = (full_key, motif_len)
    valid = {}
    total_checked = 0
    total_valid = 0
    for chrom in CHROMOSOMES:
        vcf_path = VCF_DIR / file_suffix.format(chrom)
        if not vcf_path.exists():
            print(f"  警告: {vcf_path} 不存在，跳过")
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
                chrom_vcf = fields[0]
                pos = fields[1]
                vid = fields[2]
                key3 = (chrom_vcf, pos, vid)
                if key3 not in key3_to_info:
                    continue
                full_key, motif_len = key3_to_info[key3]
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
                # 使用绝对值判断
                is_valid = any(v is not None and abs(v) >= motif_len for v in sample_vals)
                valid[full_key] = is_valid
                total_checked += 1
                if is_valid:
                    total_valid += 1
    print(f"  已评估 {total_checked} 个变异，满足条件（绝对值 ≥ motif_len）: {total_valid}")
    return valid

def process_cmp_file(filepath, position, valid_dict, sv_map, outpath, is_num_file):
    """
    处理一个比较文件：
      - 输出所有通过过滤的行（保留原始行）到 outpath
      - 若 is_num_file=False: 返回该文件中所有有效变异键的集合（用于分母）
      - 若 is_num_file=True: 返回字典 {key: 该变异的全局最大相关性}（用于分子）
    valid_dict: {key: bool}
    """
    print(f"处理 {'分子' if is_num_file else '分母'} 文件: {filepath}")
    out_f = open(outpath, 'w')
    keys_set = set()
    key_to_max = defaultdict(float)
    skipped_no_info = 0
    skipped_invalid = 0
    kept_lines = 0

    with open(filepath, 'r') as f:
        for line in f:
            original_line = line.rstrip('\n')
            # 保留注释行（以#开头）但通常不处理数据
            if original_line.startswith('#'):
                out_f.write(original_line + '\n')
                continue
            if not original_line.strip():
                out_f.write(original_line + '\n')
                continue
            parts = original_line.split('\t')
            if len(parts) < 13:
                continue
            # 跳过可能的表头行
            if any(c in parts[0].lower() for c in ['chrom', 'bed']):
                out_f.write(original_line + '\n')
                continue

            # 提取 key
            if position == 'left':
                key_cols = parts[:5]
            else:
                key_cols = parts[5:10]
            key = tuple(str(x) for x in key_cols)

            # 检查有效性
            if key not in sv_map:
                skipped_no_info += 1
                continue
            if key not in valid_dict or not valid_dict[key]:
                skipped_invalid += 1
                continue

            # 写入输出文件
            out_f.write(original_line + '\n')
            kept_lines += 1

            if not is_num_file:
                keys_set.add(key)
            else:
                # 分子文件：提取相关性最大值
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
                    curr = 0.0
                if curr > key_to_max[key]:
                    key_to_max[key] = curr

    out_f.close()
    if is_num_file:
        print(f"  输出 {kept_lines} 行，记录 {len(key_to_max)} 个变异的有效相关性")
        return key_to_max
    else:
        print(f"  输出 {kept_lines} 行，包含 {len(keys_set)} 个唯一变异")
        return keys_set

def main():
    # 1. 加载 svinfo 映射
    map_1ktr = load_svinfo_map(SVINFO_1KTR, "1ktr")
    map_nwtr = load_svinfo_map(SVINFO_NWTR, "nwtr")

    # 2. 构建有效标志（绝对值条件）
    print("\n===== 构建 1ktr 变异有效性 =====")
    valid_1ktr = build_valid_map(map_1ktr, "chr{}.1ktr.len.vcf.gz")
    print("\n===== 构建 nwtr 变异有效性 =====")
    valid_nwtr = build_valid_map(map_nwtr, "chr{}.nwtr.nomac.len.vcf.gz")

    # 3. 处理所有比较文件（只保留与25tr比较的）
    print("\n===== 处理 1ktr 分母文件（vs 25tr）=====")
    denom1_25 = process_cmp_file(CMP_1KTR_DENOM_25TR, 'left', valid_1ktr, map_1ktr, OUT_1KTR_DENOM_25TR, is_num_file=False)
    print("\n===== 处理 1ktr 分子文件（vs 25tr）=====")
    num1_25   = process_cmp_file(CMP_1KTR_NUM_25TR,   'left', valid_1ktr, map_1ktr, OUT_1KTR_NUM_25TR,   is_num_file=True)

    print("\n===== 处理 nwtr 分母文件（vs 25tr）=====")
    denomn_25 = process_cmp_file(CMP_NWTR_DENOM_25TR, 'right', valid_nwtr, map_nwtr, OUT_NWTR_DENOM_25TR, is_num_file=False)
    print("\n===== 处理 nwtr 分子文件（vs 25tr）=====")
    numn_25   = process_cmp_file(CMP_NWTR_NUM_25TR,   'right', valid_nwtr, map_nwtr, OUT_NWTR_NUM_25TR,   is_num_file=True)

    # 4. 合并分母键（并集）
    denom_1ktr = denom1_25   # 只有与25tr比较的
    denom_nwtr = denomn_25
    print(f"1ktr 分母去重总键数: {len(denom_1ktr)}")
    print(f"nwtr 分母去重总键数: {len(denom_nwtr)}")

    # 5. 合并分子相关性（每个变异取最大值，此处只有一个来源，但保持结构）
    all_corr_1 = num1_25
    all_corr_nw = numn_25
    print(f"1ktr 分子有效变异数: {len(all_corr_1)}")
    print(f"nwtr 分子有效变异数: {len(all_corr_nw)}")

    # 6. 统计加权均值
    def compute_weighted_stats(sv_map, denom_keys, corr_dict, valid_dict):
        denom_count = defaultdict(int)
        sum_corr = defaultdict(float)
        for key in denom_keys:
            if key not in sv_map or key not in valid_dict or not valid_dict[key]:
                continue
            rg, reg = sv_map[key][:2]
            denom_count[(rg, reg)] += 1
        for key, corr in corr_dict.items():
            if key not in sv_map or key not in valid_dict or not valid_dict[key]:
                continue
            rg, reg = sv_map[key][:2]
            sum_corr[(rg, reg)] += corr
        result = {}
        all_keys = set(denom_count.keys()) | set(sum_corr.keys())
        for k in all_keys:
            d = denom_count.get(k, 0)
            s = sum_corr.get(k, 0.0)
            val = s / d if d > 0 else 0.0
            result[k] = (d, s, val)
        return result

    stats_1ktr = compute_weighted_stats(map_1ktr, denom_1ktr, all_corr_1, valid_1ktr)
    stats_nwtr = compute_weighted_stats(map_nwtr, denom_nwtr, all_corr_nw, valid_nwtr)

    total_denom1 = sum(d for d,_,_ in stats_1ktr.values())
    total_sum1   = sum(s for _,s,_ in stats_1ktr.values())
    total_denom2 = sum(d for d,_,_ in stats_nwtr.values())
    total_sum2   = sum(s for _,s,_ in stats_nwtr.values())
    print(f"\n1ktr 分母总数（去重）: {total_denom1}, 分子加权和: {total_sum1:.6f}, 加权均值: {total_sum1/total_denom1 if total_denom1>0 else 0:.6f}")
    print(f"nwtr 分母总数（去重）: {total_denom2}, 分子加权和: {total_sum2:.6f}, 加权均值: {total_sum2/total_denom2 if total_denom2>0 else 0:.6f}")

    # 7. 准备绘图数据
    x_positions = np.arange(len(REFLEN_GROUPS) * len(REGIONS))
    y_1ktr = []
    y_nwtr = []
    stats_1ktr_out = []
    stats_nwtr_out = []

    for rg in REFLEN_GROUPS:
        for reg in REGIONS:
            key = (rg, reg)
            d1, s1, v1 = stats_1ktr.get(key, (0, 0.0, 0.0))
            d2, s2, v2 = stats_nwtr.get(key, (0, 0.0, 0.0))
            y_1ktr.append(v1)
            discard = (reg == "SD") or (d2 < 5)
            if discard:
                y_nwtr.append(np.nan)
            else:
                y_nwtr.append(v2)
            stats_1ktr_out.append((rg, reg, d1, s1, v1))
            stats_nwtr_out.append((rg, reg, d2, s2, v2, discard))

    # 8. 绘图
    fig, ax = plt.subplots(figsize=(16, 6))
    ax.plot(x_positions, y_1ktr, 'b-o', linewidth=2, markersize=6, label='_nolegend_')
    ax.plot(x_positions, y_nwtr, 'r-o', linewidth=2, markersize=6, label='_nolegend_')

    x_labels = [f"{rg}\n{reg}" for rg in REFLEN_GROUPS for reg in REGIONS]
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=8)
    ax.set_xlim(-0.5, len(x_positions)-0.5)
    ax.set_ylim(0, 1.05)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylabel('')
    legend = ax.get_legend()
    if legend:
        legend.remove()
    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, bbox_inches='tight', facecolor='white', format='pdf')
    print(f"\n图片已保存: {OUTPUT_PLOT}")

    # 9. 输出统计文件
    with open(OUTPUT_STAT, 'w') as f:
        f.write("# 加权均值统计（分母=去重变异数，分子=每个变异最大相关性之和）\n")
        f.write("# 过滤条件: motif>=2, reflen>=1, 原始VCF中存在绝对值≥motif_len的样本值\n")
        f.write("# 仅使用与25tr比较的结果\n")
        f.write("# 1ktr 结果\n")
        f.write("# reflen_group\tregion\tdenom\tweighted_sum\tweighted_mean\n")
        for (rg, reg, d, s, v) in stats_1ktr_out:
            f.write(f"{rg}\t{reg}\t{d}\t{s:.6f}\t{v:.6f}\n")
        f.write("\n# nwtr 结果（剔除SD区域和分母<5的点）\n")
        f.write("# reflen_group\tregion\tdenom\tweighted_sum\tweighted_mean\tkept_in_plot\n")
        for (rg, reg, d, s, v, discard) in stats_nwtr_out:
            kept = "no" if discard else "yes"
            f.write(f"{rg}\t{reg}\t{d}\t{s:.6f}\t{v:.6f}\t{kept}\n")
        f.write(f"\n# 1ktr 分母总键数: {total_denom1}\n")
        f.write(f"# 1ktr 分子加权总和: {total_sum1:.6f}\n")
        f.write(f"# nwtr 分母总键数: {total_denom2}\n")
        f.write(f"# nwtr 分子加权总和: {total_sum2:.6f}\n")
    print(f"统计文件已保存: {OUTPUT_STAT}")

    # 打印分组统计
    print("\n各 REFLEN 分组 1ktr 变异总数（有效统计）：")
    for rg in REFLEN_GROUPS:
        for reg in REGIONS:
            d = stats_1ktr.get((rg, reg), (0,0,0))[0]
            if d > 0:
                print(f"  {rg} {reg}: {d}")

if __name__ == "__main__":
    main()
