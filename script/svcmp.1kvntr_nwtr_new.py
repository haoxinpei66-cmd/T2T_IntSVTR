#!/usr/bin/env python3
"""
VNTR变异相关性分析脚本 - 只输出有overlap的变异对
输入输出文件都在脚本内部定义
改进版：基于实际序列长度计算相关性，正确计算AC计数
"""

import gzip
import sys
import numpy as np
from zbio import io

def ints(s):
    """安全转换字符串为整数，无法转换时返回原字符串"""
    try: 
        return int(s)
    except: 
        return s

def get_end_from_info(info_str):
    """从INFO字段中提取END位置"""
    for info in info_str.split(';'):
        if '=' in info:
            key, value = info.split('=', 1)
            if key == 'END':
                try:
                    return int(value)
                except:
                    return None
    return None

def calculate_overlap_ratio(pos1, end1, pos2, end2):
    """计算两个变异之间的overlap比例"""
    if end1 is None or end2 is None:
        return 0
    
    # 计算overlap长度
    start_max = max(pos1, pos2)
    end_min = min(end1, end2)
    
    if end_min < start_max:
        overlap_length = 0
    else:
        overlap_length = end_min - start_max
    
    # 计算两个变异的长度
    length1 = end1 - pos1 + 1 if end1 >= pos1 else 0
    length2 = end2 - pos2 + 1 if end2 >= pos2 else 0
    
    # 取较大的长度作为分母
    max_length = max(length1, length2)
    
    if max_length == 0:
        return 0
    
    return overlap_length / max_length

def calculate_length_correlation(ref_gts, input_gts, ref_alleles, input_alleles):
    """计算两个基因型列表的长度相关性（基于实际序列长度）"""
    # 计算每个等位基因的实际长度
    ref_lengths = [len(allele) for allele in ref_alleles]  # 包括ref和所有alt
    input_lengths = [len(allele) for allele in input_alleles]  # 包括ref和所有alt
    
    l1, l2 = [], []
    valid_samples = 0
    
    for i in range(min(len(ref_gts), len(input_gts))):
        ref_gt = ref_gts[i].replace('|', '/').split('/')
        input_gt = input_gts[i].replace('|', '/').split('/')
        
        # 跳过缺失基因型
        if (ref_gt[0] == '.' or input_gt[0] == '.' or 
            not ref_gt[0].isdigit() or not input_gt[0].isdigit()):
            continue
        
        try:
            # 计算参考变异的样本总长度
            ref_total = 0
            for allele_str in ref_gt:
                if allele_str.isdigit():
                    allele_idx = int(allele_str)
                    if allele_idx < len(ref_lengths):
                        ref_total += ref_lengths[allele_idx]
            
            # 计算输入变异的样本总长度
            input_total = 0
            for allele_str in input_gt:
                if allele_str.isdigit():
                    allele_idx = int(allele_str)
                    if allele_idx < len(input_lengths):
                        input_total += input_lengths[allele_idx]
            
            l1.append(ref_total)
            l2.append(input_total)
            valid_samples += 1
            
        except (ValueError, IndexError, TypeError):
            continue
    
    # 计算相关性
    if len(l1) > 1:
        try:
            correlation = np.corrcoef(l1, l2)[0][1]
            #lensumr2 = correlation ** 2
            lensumr2 = correlation
        except:
            lensumr2 = float('nan')
    else:
        lensumr2 = float('nan')
    
    return lensumr2, valid_samples

def calculate_ac_counts(gts, alt_count):
    """计算AC计数，每个alt对应一个AC值（包括0值）"""
    # 初始化AC计数数组，长度为alt_count，全部设为0
    ac_counts = [0] * alt_count
    
    for gt_str in gts:
        gt = gt_str.replace('|', '/').split('/')
        if gt[0] == '.' or not gt[0].isdigit():
            continue
        
        try:
            # 统计每个alt等位基因的出现次数
            # alt等位基因的索引从1开始（0是ref）
            for allele in gt:
                if allele.isdigit():
                    allele_idx = int(allele)
                    if 1 <= allele_idx <= alt_count:  # alt等位基因索引从1开始
                        ac_counts[allele_idx - 1] += 1
        except (ValueError, IndexError):
            continue
    
    return ac_counts

def calculate_an(gts):
    """计算等位基因数（AN）"""
    an = 0
    for gt_str in gts:
        gt = gt_str.replace('|', '/').split('/')
        if gt[0] != '.' and gt[0].isdigit():
            an += 2  # 二倍体基因组
    return an

def read_vcf_file(vcf_file, target_chr):
    """读取VCF文件并返回变异字典"""
    variants = {}
    sample_indices = None
    samples_list = []
    
    print(f"正在读取VCF文件: {vcf_file}")
    
    try:
        with gzip.open(vcf_file, 'rt') as f:
            for line in f:
                line = line.strip()
                
                if line.startswith('##'):
                    continue
                    
                if line.startswith('#CHROM'):
                    # 处理头行，获取所有样本
                    headers = line.split('\t')
                    if len(headers) > 9:
                        samples_list = headers[9:]
                        sample_indices = list(range(9, len(headers)))
                    continue
                
                if sample_indices is None:
                    continue
                    
                fields = line.split('\t')
                if len(fields) < 10:  # 至少要有样本列
                    continue
                    
                chrom = fields[0]
                # 统一染色体格式比较（去掉chr前缀进行比较）
                chrom_compare = chrom.replace('chr', '')
                if chrom_compare != target_chr:
                    continue
                    
                try:
                    pos = int(fields[1])
                    var_id = fields[2]
                    ref = fields[3]
                    alt = fields[4]
                    info = fields[7]
                    
                    # 获取END位置
                    end = get_end_from_info(info)
                    if end is None:
                        # 如果没有END信息，使用REF长度估计
                        end = pos + len(ref) - 1
                    
                    # 提取所有样本的基因型
                    gts = []
                    for idx in sample_indices:
                        if idx < len(fields):
                            gt_field = fields[idx].split(':')[0]
                            gts.append(gt_field)
                        else:
                            gts.append('./.')
                    
                    # 计算alt数量
                    alt_count = len(alt.split(',')) if alt != '.' else 0
                    
                    # 使用位置和等位基因作为唯一键
                    key = (chrom, pos, ref, alt)
                    variants[key] = {
                        'chrom': chrom,  # 保留原始染色体名称（可能带chr）
                        'pos': pos,
                        'id': var_id,
                        'ref': ref,
                        'alt': alt,
                        'end': end,
                        'gts': gts,
                        'info': info,
                        'alt_count': alt_count
                    }
                    
                except (ValueError, IndexError) as e:
                    continue
                    
    except Exception as e:
        print(f"读取文件 {vcf_file} 时出错: {e}")
        return {}
    
    print(f"从 {vcf_file} 读取到 {len(variants)} 个变异")
    return variants, samples_list

def find_overlapping_variants(ref_variants, input_variants, min_overlap_ratio=0):
    """查找有overlap的变异对"""
    overlapping_pairs = []
    
    print("正在查找有overlap的变异对...")
    
    # 将输入变异按位置排序以便快速查找
    input_list = sorted(input_variants.values(), key=lambda x: x['pos'])
    
    for ref_key, ref_var in ref_variants.items():
        ref_pos = ref_var['pos']
        ref_end = ref_var['end']
        
        # 在输入变异中查找重叠的变异
        for input_var in input_list:
            input_pos = input_var['pos']
            input_end = input_var['end']
            
            # 如果输入变异的位置已经超过参考变异的结束位置，可以提前终止
            if input_pos > ref_end:
                break
                
            # 计算overlap比例
            overlap_ratio = calculate_overlap_ratio(ref_pos, ref_end, input_pos, input_end)
            
            if overlap_ratio >= min_overlap_ratio:
                overlapping_pairs.append((ref_var, input_var, overlap_ratio))
    
    print(f"找到 {len(overlapping_pairs)} 对有overlap的变异")
    return overlapping_pairs

def format_ac_counts(ac_counts):
    """格式化AC计数，保留所有值（包括0）"""
    return ','.join(map(str, ac_counts))

def main():
    # ==================== 输入文件配置 ====================
    # 第一个VCF文件（参考文件）
    REF_VCF_FILE = "/Parastor300s_G30S/haoxp/third/cmp/1019/extract/modif_vntr_967_vcf/chrsurban.processed.vcf.gz"
    
    # 第二个VCF文件（输入文件）
    INPUT_VCF_FILE = "/Parastor300s_G30S/haoxp/third/cmp/1019/extract/tr.mac1.chrsurban.vcf.gz"
    
    # 输出文件
    OUTPUT_FILE = "/Parastor300s_G30S/haoxp/third/cmp/1019/svcmp_1kvntr_nwtr/svcmp.1kvntr_nwtr_new_chrsurban.txt"
    
    # 目标染色体（比较时使用，不包含chr前缀）
    TARGET_CHR = "surban"
    
    # 最小overlap比例阈值（>0.1才输出）
    MIN_OVERLAP_RATIO = 0.1
    # ==================== 配置结束 ====================
    
    # 读取两个VCF文件
    ref_variants, ref_samples = read_vcf_file(REF_VCF_FILE, TARGET_CHR)
    input_variants, input_samples = read_vcf_file(INPUT_VCF_FILE, TARGET_CHR)
    
    if not ref_variants or not input_variants:
        print("错误: 无法读取VCF文件")
        return
    
    print(f"参考文件样本数: {len(ref_samples)}")
    print(f"输入文件样本数: {len(input_samples)}")
    
    # 查找有overlap的变异对
    overlapping_pairs = find_overlapping_variants(ref_variants, input_variants, MIN_OVERLAP_RATIO)
    
    if not overlapping_pairs:
        print("没有找到有overlap的变异对")
        return
    
    # 打开输出文件
    print(f"正在写入结果到: {OUTPUT_FILE}")
    with open(OUTPUT_FILE, 'w') as outfile:
        # 写入头信息（按新顺序）
        header = [
            'CHROM1', 'POS1', 'ID1', 'END1', 'REF1', 'ALT1', 'AN1', 'AC1', 
            'VALID_SAMPLES', 'LEN_R2', 'CHROM2', 'POS2', 'ID2', 'END2', 
            'REF2', 'ALT2', 'AN2', 'AC2', 'OVERLAP_RATIO'
        ]
        outfile.write('\t'.join(header) + '\n')
        
        processed_count = 0
        valid_count = 0
        
        for ref_var, input_var, overlap_ratio in overlapping_pairs:
            processed_count += 1
            
            # 构建等位基因列表（包括ref和所有alt）
            ref_alleles = [ref_var['ref']] + ref_var['alt'].split(',')
            input_alleles = [input_var['ref']] + input_var['alt'].split(',')
            
            # 计算相关性（只使用两个文件共有的样本）
            min_samples = min(len(ref_var['gts']), len(input_var['gts']))
            ref_gts = ref_var['gts'][:min_samples]
            input_gts = input_var['gts'][:min_samples]
            
            # 计算长度相关性（基于实际序列长度）
            lensumr2, valid_samples = calculate_length_correlation(
                ref_gts, input_gts, ref_alleles, input_alleles)
            
            # 计算各自的AC计数
            ac1 = calculate_ac_counts(ref_gts, ref_var['alt_count'])
            ac2 = calculate_ac_counts(input_gts, input_var['alt_count'])
            
            # 计算AN
            an1 = calculate_an(ref_gts)
            an2 = calculate_an(input_gts)
            
            # 只输出有有效相关性的结果，且overlap比例>0.1
            if not np.isnan(lensumr2) and valid_samples > 0 and overlap_ratio > 0.1:
                valid_count += 1
                
                # 格式化AC计数（保留所有值，包括0）
                ac1_str = format_ac_counts(ac1)
                ac2_str = format_ac_counts(ac2)
                
                # 验证AC计数数量是否正确
                expected_ac1_count = ref_var['alt_count']
                actual_ac1_count = len(ac1)
                expected_ac2_count = input_var['alt_count']
                actual_ac2_count = len(ac2)
                
                if expected_ac1_count != actual_ac1_count:
                    print(f"警告: 变异1 AC计数数量不匹配: 期望{expected_ac1_count}, 实际{actual_ac1_count}")
                if expected_ac2_count != actual_ac2_count:
                    print(f"警告: 变异2 AC计数数量不匹配: 期望{expected_ac2_count}, 实际{actual_ac2_count}")
                
                # 准备输出行（按新顺序）
                output_fields = [
                    ref_var['chrom'],                    # CHROM1（带chr）
                    str(ref_var['pos']),                # POS1
                    ref_var['id'],                      # ID1
                    str(ref_var['end']),                # END1
                    ref_var['ref'],                     # REF1
                    ref_var['alt'],                     # ALT1
                    str(an1),                           # AN1
                    ac1_str,                            # AC1
                    str(valid_samples),                 # VALID_SAMPLES
                    f"{lensumr2:.6f}",                  # LEN_R2
                    input_var['chrom'],                 # CHROM2（带chr）
                    str(input_var['pos']),              # POS2
                    input_var['id'],                    # ID2
                    str(input_var['end']),              # END2
                    input_var['ref'],                   # REF2
                    input_var['alt'],                   # ALT2
                    str(an2),                           # AN2
                    ac2_str,                            # AC2
                    f"{overlap_ratio:.6f}"              # OVERLAP_RATIO
                ]
                
                outfile.write('\t'.join(output_fields) + '\n')
            
            if processed_count % 100 == 0:
                print(f"已处理 {processed_count} 对变异，有效 {valid_count} 对")
        
        print(f"处理完成: 总共 {len(overlapping_pairs)} 对有overlap的变异，有效 {valid_count} 对")
        if len(overlapping_pairs) > 0:
            print(f"有效率: {valid_count/len(overlapping_pairs)*100:.2f}%")

if __name__ == "__main__":
    main()
