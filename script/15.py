#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from pathlib import Path
import gzip

# ==================== 路径配置 ====================
BASE_DIR = Path("/Parastor300s_G30S/haoxp/third/cmp/25")
SVINFO_25TR = BASE_DIR / "svinfo.25tr.nomac.rmnoreflen.addregion.txt"
BED_FILE = BASE_DIR / "extract/len_vcf/25tr.len.nomac.rmnoreflen.bed.gz"

# 输入文件
CMP_1KTR_DENOM = BASE_DIR / "svcmp.1ktr.25tr.nomac.rmnoreflen.addspear.over05.txt"
CMP_1KTR_NUM   = BASE_DIR / "svcmp.1ktr.25tr.nomac.rmnoreflen.addspear.over05.or05.txt"
CMP_NWTR_DENOM = BASE_DIR / "svcmp.25tr.nwtr.nomac.rmnoreflen.addspear.over05.txt"
CMP_NWTR_NUM   = BASE_DIR / "svcmp.25tr.nwtr.nomac.rmnoreflen.addspear.over05.or05.txt"

# 输出过滤后的文件（添加 .25macmotif 后缀）
OUT_1KTR_DENOM = BASE_DIR / "svcmp.1ktr.25tr.nomac.rmnoreflen.addspear.over05.25macmotif.txt"
OUT_1KTR_NUM   = BASE_DIR / "svcmp.1ktr.25tr.nomac.rmnoreflen.addspear.over05.or05.25macmotif.txt"
OUT_NWTR_DENOM = BASE_DIR / "svcmp.25tr.nwtr.nomac.rmnoreflen.addspear.over05.25macmotif.txt"
OUT_NWTR_NUM   = BASE_DIR / "svcmp.25tr.nwtr.nomac.rmnoreflen.addspear.over05.or05.25macmotif.txt"

OUTPUT_PLOT = BASE_DIR / "1ktr_nwtr.reflen.recall.base05.or05.macmotif.rmnoreflen.pdf"
OUTPUT_STAT = BASE_DIR / "1ktr_nwtr.reflen.recall.base05.or05.macmotif.rmnoreflen.txt"

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

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False

def get_reflen_group(reflen):
    for low, high, label in REFLEN_BINS:
        if low <= reflen <= high:
            return label
    return ">200"

def load_svinfo_map(filepath):
    print(f"读取 svinfo: {filepath}")
    mapping = {}
    seen = set()
    with open(filepath) as f:
        f.readline()
        for line in f:
            line = line.rstrip('\n')
            if not line: continue
            parts = line.split('\t')
            if len(parts) < 6: continue
            chrom, pos, vid, reflen_str, motif_str, region = parts[:6]
            try:
                motif = int(motif_str)
            except:
                continue
            if motif < 2: continue
            key = (chrom, pos, vid, reflen_str, motif_str)
            if key in seen: continue
            seen.add(key)
            try:
                reflen = int(reflen_str)
            except:
                reflen = 0
            group = get_reflen_group(reflen)
            mapping[key] = (group, region, motif)
    print(f"  加载映射数: {len(mapping):,}")
    return mapping

def load_bed_samples(bed_path):
    print(f"读取 BED 文件: {bed_path}")
    bed_map = {}
    with gzip.open(bed_path, 'rt') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            parts = line.split('\t')
            if len(parts) < 5: continue
            chrom = parts[0]
            if chrom not in [f"chr{i}" for i in range(1, 23)]: continue
            start = int(parts[1])
            sample_vals = []
            for val_str in parts[4:]:
                if val_str in ('.', '-', './.', '.|.'):
                    sample_vals.append(None)
                else:
                    try:
                        sample_vals.append(float(val_str))
                    except ValueError:
                        sample_vals.append(None)
            if len(sample_vals) != 25:
                print(f"  警告: {chrom}:{start} 样本数 {len(sample_vals)}，期望25")
            bed_map[(chrom, start)] = sample_vals
    print(f"  加载 BED 条目数: {len(bed_map):,}")
    return bed_map

def build_valid_map(svinfo_map, bed_map):
    valid = {}
    total = len(svinfo_map)
    valid_count = 0
    missing_bed = 0
    print(f"  开始评估 {total:,} 个变异...")
    for key, (_, _, motif_len) in svinfo_map.items():
        chrom, pos, vid, reflen, motif_str = key
        bed_key = (chrom, int(pos))
        if bed_key not in bed_map:
            valid[key] = False
            missing_bed += 1
            continue
        sample_vals = bed_map[bed_key]
        ok = any(v is not None and abs(v) >= motif_len for v in sample_vals)
        valid[key] = ok
        if ok:
            valid_count += 1
    print(f"  有效变异数（满足样本条件）: {valid_count}，缺失 BED 条目: {missing_bed}")
    return valid

def process_cmp_file(input_path, output_path, position, svinfo_map, valid_map, label, is_num_file=False):
    print(f"处理 {label}: {input_path}")
    if is_num_file:
        # 分子文件：记录每个变异的最大 max_corr
        per_key_max = {}
    else:
        # 分母文件：计数唯一键
        denom_counts = defaultdict(int)
        seen_keys = set()
    out_f = open(output_path, 'w')
    output_line_count = 0
    skip_no_svinfo = 0
    skip_invalid = 0

    with open(input_path) as f:
        for line in f:
            original_line = line.rstrip('\n')
            if not original_line or original_line.startswith('#') or original_line.startswith('bed'):
                continue
            parts = original_line.split('\t')
            if len(parts) < 13:
                continue
            cols = parts[:5] if position == 'left' else parts[5:10]
            if any(('chrom' in c.lower() or 'pos' in c.lower()) for c in cols):
                continue
            key = tuple(cols)
            if key not in svinfo_map:
                skip_no_svinfo += 1
                continue
            if not valid_map.get(key, False):
                skip_invalid += 1
                continue

            # 提取相关性最大值
            corr1 = None
            corr2 = None
            try:
                if parts[11] not in ('.','-','NA',''):
                    corr1 = float(parts[11])
            except:
                pass
            try:
                if parts[12] not in ('.','-','NA',''):
                    corr2 = float(parts[12])
            except:
                pass
            max_corr = 0.0
            if corr1 is not None and corr2 is not None:
                max_corr = max(corr1, corr2)
            elif corr1 is not None:
                max_corr = corr1
            elif corr2 is not None:
                max_corr = corr2

            out_f.write(original_line + '\n')
            output_line_count += 1

            if is_num_file:
                # 分子：记录该变异出现的最大 max_corr
                current = per_key_max.get(key, -float('inf'))
                if max_corr > current:
                    per_key_max[key] = max_corr
            else:
                # 分母：唯一键计数
                if key not in seen_keys:
                    seen_keys.add(key)
                    rg, reg = svinfo_map[key][:2]
                    denom_counts[(rg, reg)] += 1

    out_f.close()

    if is_num_file:
        # 汇聚分子：按组累加每个变异的最终 max_corr
        num_sums = defaultdict(float)
        for key, max_val in per_key_max.items():
            rg, reg = svinfo_map[key][:2]
            num_sums[(rg, reg)] += max_val
        total_valid = len(per_key_max)
        print(f"  唯一有效键数（分子）: {total_valid}, 输出行数: {output_line_count}, 跳过(无svinfo): {skip_no_svinfo}, 跳过(样本条件不满足): {skip_invalid}")
        return num_sums, total_valid
    else:
        total_valid = sum(denom_counts.values())
        print(f"  唯一有效键数（分母）: {total_valid}, 输出行数: {output_line_count}, 跳过(无svinfo): {skip_no_svinfo}, 跳过(样本条件不满足): {skip_invalid}")
        return denom_counts, total_valid

def compute_weighted_recall(denom_counts, num_sums):
    all_keys = set(denom_counts.keys()) | set(num_sums.keys())
    result = {}
    for key in all_keys:
        d = denom_counts.get(key, 0)
        n = num_sums.get(key, 0.0)
        r = n / d if d > 0 else 0.0
        result[key] = (d, n, r)
    return result

def main():
    svinfo_map = load_svinfo_map(SVINFO_25TR)
    bed_map = load_bed_samples(BED_FILE)
    print("\n===== 构建 25tr 变异有效性 =====")
    valid_25tr = build_valid_map(svinfo_map, bed_map)

    print("\n===== 处理 1ktr 分母文件 =====")
    denom_1ktr, _ = process_cmp_file(CMP_1KTR_DENOM, OUT_1KTR_DENOM, 'right', svinfo_map, valid_25tr, "1ktr分母", is_num_file=False)
    print("\n===== 处理 1ktr 分子文件 =====")
    num_1ktr, _ = process_cmp_file(CMP_1KTR_NUM, OUT_1KTR_NUM, 'right', svinfo_map, valid_25tr, "1ktr分子", is_num_file=True)

    print("\n===== 处理 nwtr 分母文件 =====")
    denom_nwtr, _ = process_cmp_file(CMP_NWTR_DENOM, OUT_NWTR_DENOM, 'left', svinfo_map, valid_25tr, "nwtr分母", is_num_file=False)
    print("\n===== 处理 nwtr 分子文件 =====")
    num_nwtr, _ = process_cmp_file(CMP_NWTR_NUM, OUT_NWTR_NUM, 'left', svinfo_map, valid_25tr, "nwtr分子", is_num_file=True)

    recall_1ktr = compute_weighted_recall(denom_1ktr, num_1ktr)
    recall_nwtr = compute_weighted_recall(denom_nwtr, num_nwtr)

    total_denom1 = sum(d for d,_,_ in recall_1ktr.values())
    total_num1   = sum(n for _,n,_ in recall_1ktr.values())
    total_denom2 = sum(d for d,_,_ in recall_nwtr.values())
    total_num2   = sum(n for _,n,_ in recall_nwtr.values())
    print(f"\n1ktr 分母总键数（去重）: {total_denom1}, 分子总和（每个变异最大相关性之和）: {total_num1:.6f}")
    print(f"nwtr 分母总键数（去重）: {total_denom2}, 分子总和（每个变异最大相关性之和）: {total_num2:.6f}")

    x_positions = np.arange(len(REFLEN_GROUPS) * len(REGIONS))
    y_1ktr = []
    y_nwtr = []
    stats = []
    for rg in REFLEN_GROUPS:
        for reg in REGIONS:
            key = (rg, reg)
            d1, n1, r1 = recall_1ktr.get(key, (0, 0.0, 0.0))
            d2, n2, r2 = recall_nwtr.get(key, (0, 0.0, 0.0))
            y_1ktr.append(r1)
            y_nwtr.append(np.nan if reg == "SD" else r2)
            stats.append((rg, reg, d1, n1, r1, d2, n2, r2))

    fig, ax = plt.subplots(figsize=(15, 6))
    ax.plot(x_positions, y_1ktr, 'b-o', linewidth=2, markersize=6, label='_nolegend_')
    ax.plot(x_positions, y_nwtr, 'r-o', linewidth=2, markersize=6, label='_nolegend_')
    x_labels = [f"{rg}\n{reg}" for rg in REFLEN_GROUPS for reg in REGIONS]
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=9)
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
    print(f"图片已保存: {OUTPUT_PLOT}")

    with open(OUTPUT_STAT, 'w') as f:
        f.write("# 加权 Recall 统计（分母去重，分子每个变异只取最大相关性并求和）\n")
        f.write("# reflen_group\tregion\t1ktr_denom\t1ktr_sum_corr\t1ktr_weighted_recall\tnwtr_denom\tnwtr_sum_corr\tnwtr_weighted_recall\n")
        for (rg, reg, d1, n1, r1, d2, n2, r2) in stats:
            f.write(f"{rg}\t{reg}\t{d1}\t{n1:.6f}\t{r1:.6f}\t{d2}\t{n2:.6f}\t{r2:.6f}\n")
        f.write(f"\n# 1ktr 分母总键数（去重）: {total_denom1}\n")
        f.write(f"# 1ktr 分子总和（每个变异最大相关性之和）: {total_num1:.6f}\n")
        f.write(f"# nwtr 分母总键数（去重）: {total_denom2}\n")
        f.write(f"# nwtr 分子总和（每个变异最大相关性之和）: {total_num2:.6f}\n")
    print(f"统计文件已保存: {OUTPUT_STAT}")

if __name__ == "__main__":
    main()
