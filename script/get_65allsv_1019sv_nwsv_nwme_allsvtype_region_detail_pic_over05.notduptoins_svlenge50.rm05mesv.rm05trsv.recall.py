#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, Tuple, List, Set, Optional
import sys
import time
from tqdm import tqdm
import os

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

try:
    import intervaltree
    INTERVALTREE_AVAILABLE = True
except ImportError:
    print("警告：intervaltree模块未安装，将使用较慢的原始方法")
    INTERVALTREE_AVAILABLE = False

# ==================== 原有辅助函数（保持不变） ====================
def normalize_sv_type(sv_type: str) -> str:
    """统一变异类型命名"""
    if pd.isna(sv_type):
        return sv_type
    sv_type = str(sv_type).strip()
    type_mapping = {'ALU': 'Alu', 'alu': 'Alu', 'L1': 'LINE', 'l1': 'LINE'}
    return type_mapping.get(sv_type, sv_type)

def build_interval_tree(regions: List[Tuple[str, int, int]]) -> Dict[str, intervaltree.IntervalTree]:
    """构建区间树用于快速区间查询"""
    tree_dict = {}
    for chrom, start, end in regions:
        if chrom not in tree_dict:
            tree_dict[chrom] = intervaltree.IntervalTree()
        tree_dict[chrom][start:end+1] = True
    return tree_dict

def load_sr_t2t_file(sr_t2t_file: str) -> Tuple[List[Tuple[str, int, int]], Dict]:
    """加载单独的SR_t2t区域文件并构建区间树"""
    print("正在加载SR_t2t区域文件...")
    sr_t2t_regions = []
    try:
        df = pd.read_csv(sr_t2t_file, sep='\t', compression='gzip', comment='#')
        required_cols = ['#chrom', 'chromStart', 'chromEnd']
        if all(col in df.columns for col in required_cols):
            for _, row in df.iterrows():
                sr_t2t_regions.append((row['#chrom'], row['chromStart'], row['chromEnd']))
        else:
            print("警告：SR_t2t文件列名不匹配，尝试使用位置索引")
            for _, row in df.iterrows():
                sr_t2t_regions.append((row.iloc[0], row.iloc[1], row.iloc[2]))
        if INTERVALTREE_AVAILABLE:
            sr_t2t_tree = build_interval_tree(sr_t2t_regions)
            print(f"SR_t2t区域文件加载完成，共 {len(sr_t2t_regions):,} 个区域（使用intervaltree优化）")
            return sr_t2t_regions, sr_t2t_tree
        else:
            print(f"SR_t2t区域文件加载完成，共 {len(sr_t2t_regions):,} 个区域（使用原始方法）")
            return sr_t2t_regions, {}
    except Exception as e:
        print(f"错误：加载SR_t2t区域文件失败 - {e}")
        return [], {}

def load_sd_regions(sd_file: str) -> Tuple[List[Tuple[str, int, int]], Dict]:
    """加载SD区域文件并构建区间树"""
    print("正在加载SD区域文件...")
    sd_regions = []
    try:
        df = pd.read_csv(sd_file, sep='\t', header=None, names=['CHROM', 'START', 'END'])
        for _, row in df.iterrows():
            sd_regions.append((row['CHROM'], row['START'], row['END']))
        if INTERVALTREE_AVAILABLE:
            sd_tree = build_interval_tree(sd_regions)
            print(f"SD区域文件加载完成，共 {len(sd_regions)} 个区域（使用intervaltree优化）")
            return sd_regions, sd_tree
        else:
            print(f"SD区域文件加载完成，共 {len(sd_regions)} 个区域（使用原始方法）")
            return sd_regions, {}
    except Exception as e:
        print(f"错误：加载SD区域文件失败 - {e}")
        return [], {}

def load_rm_regions(rm_file: str) -> Tuple:
    """加载RepeatMasker区域文件，返回SR区域和RM区域及对应的区间树"""
    print("正在加载RepeatMasker区域文件...")
    sr_regions = []
    rm_regions = []
    try:
        df = pd.read_csv(rm_file, sep='\t', header=None,
                        names=['CHROM', 'START', 'END', 'NAME', 'SCORE', 'STRAND', 'TYPE', 'FAMILY', 'DIV', 'ID'])
        total_rows = len(df)
        print(f"RepeatMasker文件共有 {total_rows:,} 行数据")
        for i, (_, row) in enumerate(tqdm(df.iterrows(), total=total_rows, desc="加载RepeatMasker",
                                         mininterval=5.0, maxinterval=10.0,
                                         bar_format='{l_bar}{bar:20}{r_bar}')):
            if row['TYPE'] == "Simple_repeat":
                sr_regions.append((row['CHROM'], row['START'], row['END']))
            else:
                rm_regions.append((row['CHROM'], row['START'], row['END']))
        if INTERVALTREE_AVAILABLE:
            sr_tree = build_interval_tree(sr_regions)
            rm_tree = build_interval_tree(rm_regions)
            print(f"RepeatMasker区域文件加载完成，共 {len(sr_regions):,} 个SR区域，{len(rm_regions):,} 个RM区域（使用intervaltree优化）")
            return sr_regions, sr_tree, rm_regions, rm_tree
        else:
            print(f"RepeatMasker区域文件加载完成，共 {len(sr_regions):,} 个SR区域，{len(rm_regions):,} 个RM区域（使用原始方法）")
            return sr_regions, {}, rm_regions, {}
    except Exception as e:
        print(f"错误：加载RepeatMasker区域文件失败 - {e}")
        return [], {}, [], {}

def point_in_region_fixed(chrom: str, pos: int, tree_dict: Dict) -> bool:
    """检查点是否在区域中（使用intervaltree优化）"""
    if chrom in tree_dict:
        intervals = tree_dict[chrom][pos]
        return len(intervals) > 0
    return False

def point_in_region_list_fixed(chrom: str, pos: int, regions_list: List[Tuple[str, int, int]]) -> bool:
    """检查点是否在区域中（使用原始列表方法）"""
    for region_chrom, region_start, region_end in regions_list:
        if chrom == region_chrom and region_start <= pos <= region_end:
            return True
    return False

def is_insertion_from_ref_alt(ref: str, alt: str) -> bool:
    """根据REF和ALT序列判断是否为insertion"""
    return len(alt) > len(ref)

def classify_variant_fixed(row: pd.Series,
                       sr_t2t_tree: Dict, sr_tree: Dict, sd_tree: Dict, rm_tree: Dict,
                       sr_t2t_regions_list: List, sr_regions_list: List, sd_regions_list: List, rm_regions_list: List) -> str:
    """分类变异 - 区域判断优先级：SR_t2t -> SR -> SD -> RM -> Unique"""
    chrom = row['CHROM']
    start = row['POS']
    svlen = abs(row['SVLEN'])
    variant_id = str(row['ID'])
    sv_type = row['SV_TYPE']
    sv_type = normalize_sv_type(sv_type)
    end = start + svlen - 1
    use_fast_method = INTERVALTREE_AVAILABLE and sr_t2t_tree and sr_tree and sd_tree and rm_tree
    # 判断是否为INS类型或特定类型且为insertion
    is_ins_type = sv_type == "INS"
    is_te_ins = sv_type in ["Alu", "LINE", "LTR", "SVA"] and is_insertion_from_ref_alt(row['REF'], row['ALT'])
    is_phase2_ins = "phase2_ins" in variant_id
    if is_ins_type or is_te_ins or is_phase2_ins:
        # 只检查起始点
        if use_fast_method:
            if point_in_region_fixed(chrom, start, sr_t2t_tree):
                return "SR"
            elif point_in_region_fixed(chrom, start, sr_tree):
                return "SR"
            elif point_in_region_fixed(chrom, start, sd_tree):
                return "SD"
            elif point_in_region_fixed(chrom, start, rm_tree):
                return "RM"
            else:
                return "Unique"
        else:
            if point_in_region_list_fixed(chrom, start, sr_t2t_regions_list):
                return "SR"
            elif point_in_region_list_fixed(chrom, start, sr_regions_list):
                return "SR"
            elif point_in_region_list_fixed(chrom, start, sd_regions_list):
                return "SD"
            elif point_in_region_list_fixed(chrom, start, rm_regions_list):
                return "RM"
            else:
                return "Unique"
    else:
        # 检查起始点和结束点
        start_in_sr_t2t = False
        start_in_sr = False
        start_in_sd = False
        start_in_rm = False
        end_in_sr_t2t = False
        end_in_sr = False
        end_in_sd = False
        end_in_rm = False
        if use_fast_method:
            start_in_sr_t2t = point_in_region_fixed(chrom, start, sr_t2t_tree)
            start_in_sr = point_in_region_fixed(chrom, start, sr_tree)
            start_in_sd = point_in_region_fixed(chrom, start, sd_tree)
            start_in_rm = point_in_region_fixed(chrom, start, rm_tree)
            end_in_sr_t2t = point_in_region_fixed(chrom, end, sr_t2t_tree)
            end_in_sr = point_in_region_fixed(chrom, end, sr_tree)
            end_in_sd = point_in_region_fixed(chrom, end, sd_tree)
            end_in_rm = point_in_region_fixed(chrom, end, rm_tree)
        else:
            start_in_sr_t2t = point_in_region_list_fixed(chrom, start, sr_t2t_regions_list)
            start_in_sr = point_in_region_list_fixed(chrom, start, sr_regions_list)
            start_in_sd = point_in_region_list_fixed(chrom, start, sd_regions_list)
            start_in_rm = point_in_region_list_fixed(chrom, start, rm_regions_list)
            end_in_sr_t2t = point_in_region_list_fixed(chrom, end, sr_t2t_regions_list)
            end_in_sr = point_in_region_list_fixed(chrom, end, sr_regions_list)
            end_in_sd = point_in_region_list_fixed(chrom, end, sd_regions_list)
            end_in_rm = point_in_region_list_fixed(chrom, end, rm_regions_list)
        if start_in_sr_t2t or end_in_sr_t2t:
            return "SR"
        elif start_in_sr or end_in_sr:
            return "SR"
        elif start_in_sd or end_in_sd:
            return "SD"
        elif start_in_rm or end_in_rm:
            return "RM"
        else:
            return "Unique"

def classify_variants_batch_fixed(df: pd.DataFrame,
                          sr_t2t_tree: Dict, sr_tree: Dict, sd_tree: Dict, rm_tree: Dict,
                          sr_t2t_regions_list: List, sr_regions_list: List, sd_regions_list: List, rm_regions_list: List,
                          desc: str = "分类变异") -> pd.Series:
    """批量分类变异"""
    print(f"正在{desc}...")
    classifications = []
    for _, row in tqdm(df.iterrows(), total=len(df), desc=desc,
                      mininterval=3.0, maxinterval=10.0,
                      bar_format='{l_bar}{bar:20}{r_bar}'):
        region_type = classify_variant_fixed(row, sr_t2t_tree, sr_tree, sd_tree, rm_tree,
                                            sr_t2t_regions_list, sr_regions_list, sd_regions_list, rm_regions_list)
        classifications.append(region_type)
    return pd.Series(classifications, index=df.index)

def read_and_classify_variants(file_path: str, dataset_name: str,
                              sr_t2t_tree: Dict, sr_tree: Dict, sd_tree: Dict, rm_tree: Dict,
                              sr_t2t_regions_list: List, sr_regions_list: List, sd_regions_list: List, rm_regions_list: List) -> Tuple[Dict, pd.DataFrame]:
    """读取变异文件并分类统计"""
    print(f"正在读取并分类 {dataset_name} 变异文件...")
    target_types = ['DEL', 'INS', 'DUP', 'INV', 'Alu', 'LINE', 'LTR', 'SVA']
    region_types = ['SD', 'SR', 'RM', 'Unique']
    stats = {sv_type: {region: 0 for region in region_types} for sv_type in target_types}
    try:
        start_time = time.time()
        print(f"正在读取 {dataset_name} 文件...")
        df = pd.read_csv(file_path, sep='\t', comment='#', header=None,
                        names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'SV_TYPE', 'SVLEN'])
        df['SV_TYPE'] = df['SV_TYPE'].apply(normalize_sv_type)
        original_count = len(df)
        df = df[df['SV_TYPE'].isin(target_types)].copy()
        df['SVLEN'] = df['SVLEN'].abs()
        print(f"{dataset_name} 文件读取完成，共 {original_count:,} 个变异，其中 {len(df):,} 个目标变异，耗时 {time.time() - start_time:.2f} 秒")
        start_time = time.time()
        df['REGION'] = classify_variants_batch_fixed(df, sr_t2t_tree, sr_tree, sd_tree, rm_tree,
                                             sr_t2t_regions_list, sr_regions_list, sd_regions_list, rm_regions_list,
                                             f"分类{dataset_name}变异")
        classification_time = time.time() - start_time
        print(f"{dataset_name} 分类完成，耗时 {classification_time:.2f} 秒，平均 {classification_time/len(df)*1000:.2f} 毫秒/变异")
        print("正在统计分类结果...")
        for sv_type in target_types:
            sv_type_df = df[df['SV_TYPE'] == sv_type]
            for region in region_types:
                stats[sv_type][region] = len(sv_type_df[sv_type_df['REGION'] == region])
        print(f"{dataset_name} 变异类型区域分布:")
        for sv_type in target_types:
            total = sum(stats[sv_type].values())
            if total > 0:
                print(f"  {sv_type}: SD={stats[sv_type]['SD']:,}, SR={stats[sv_type]['SR']:,}, RM={stats[sv_type]['RM']:,}, Unique={stats[sv_type]['Unique']:,}")
    except Exception as e:
        print(f"错误：读取并分类 {dataset_name} 变异文件失败 - {e}")
    return stats, df

def extract_variant_ids_from_cmp_files(cmp_1019_nw_file: str, cmp_65_nw_file: str, cmp_65_1019_file: str) -> Tuple[Set, Set, Set, Set]:
    """从比较文件中提取变异ID集合"""
    print("正在从比较文件中提取变异ID...")
    nw_1019_ids = set()
    kg_1019_nw_ids = set()
    nw_65_ids = set()
    kg_65_1019_ids = set()
    try:
        print("读取1019-NyuWa比较文件...")
        df_1019_nw = pd.read_csv(cmp_1019_nw_file, sep='\t', comment='#', header=None)
        if len(df_1019_nw.columns) >= 19:
            for _, row in df_1019_nw.iterrows():
                nw_chrom = row[8]
                nw_pos = row[9]
                nw_svtype = normalize_sv_type(row[18])
                if pd.notna(nw_chrom) and pd.notna(nw_pos) and pd.notna(nw_svtype):
                    nw_id = f"{nw_chrom}_{nw_pos}_{nw_svtype}"
                    nw_1019_ids.add(nw_id)
                kg_id = row[2]
                if pd.notna(kg_id):
                    kg_1019_nw_ids.add(str(kg_id))
        print(f"1019-NyuWa: 找到 {len(nw_1019_ids)} 个NyuWa变异, {len(kg_1019_nw_ids)} 个1019变异")
        print("读取65-NyuWa比较文件...")
        df_65_nw = pd.read_csv(cmp_65_nw_file, sep='\t', comment='#', header=None)
        if len(df_65_nw.columns) >= 19:
            for _, row in df_65_nw.iterrows():
                nw_chrom = row[8]
                nw_pos = row[9]
                nw_svtype = normalize_sv_type(row[18])
                if pd.notna(nw_chrom) and pd.notna(nw_pos) and pd.notna(nw_svtype):
                    nw_id = f"{nw_chrom}_{nw_pos}_{nw_svtype}"
                    nw_65_ids.add(nw_id)
        print(f"65-NyuWa: 找到 {len(nw_65_ids)} 个NyuWa变异")
        print("读取65-1019比较文件...")
        df_65_1019 = pd.read_csv(cmp_65_1019_file, sep='\t', comment='#', header=None)
        if len(df_65_1019.columns) >= 17:
            for _, row in df_65_1019.iterrows():
                kg_id = row[2]
                if pd.notna(kg_id):
                    kg_65_1019_ids.add(str(kg_id))
        print(f"65-1019: 找到 {len(kg_65_1019_ids)} 个1019变异")
    except Exception as e:
        print(f"错误：提取变异ID失败 - {e}")
    return nw_1019_ids, kg_1019_nw_ids, nw_65_ids, kg_65_1019_ids

def analyze_overlap_by_region_with_classified_data(nw_df: pd.DataFrame, kg_df: pd.DataFrame,
                                                  nw_1019_ids: Set, kg_1019_nw_ids: Set,
                                                  nw_65_ids: Set, kg_65_1019_ids: Set) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """使用已分类的区域数据来分析重叠关系"""
    print("正在分析区域重叠关系...")
    target_types = ['DEL', 'INS', 'DUP', 'INV', 'Alu', 'LINE', 'LTR', 'SVA']
    region_types = ['SD', 'SR', 'RM', 'Unique']
    nyuwa_results = []
    kg_results = []
    nw_df['UNIQUE_ID'] = nw_df['CHROM'].astype(str) + '_' + nw_df['POS'].astype(str) + '_' + nw_df['SV_TYPE']
    print("分析NyuWa数据重叠关系...")
    for sv_type in target_types:
        for region in region_types:
            region_variants = nw_df[(nw_df['SV_TYPE'] == sv_type) & (nw_df['REGION'] == region)]
            total_count = len(region_variants)
            if total_count > 0:
                in_1019 = 0
                in_65 = 0
                in_both = 0
                for _, variant in region_variants.iterrows():
                    variant_id = variant['UNIQUE_ID']
                    in_1019_flag = variant_id in nw_1019_ids
                    in_65_flag = variant_id in nw_65_ids
                    if in_1019_flag and in_65_flag:
                        in_both += 1
                    elif in_1019_flag:
                        in_1019 += 1
                    elif in_65_flag:
                        in_65 += 1
                unique_count = total_count - (in_both + in_1019 + in_65)
                nyuwa_results.append({
                    'SV_TYPE': sv_type,
                    'REGION': region,
                    'DATASET': 'NyuWa',
                    'TOTAL_COUNT': total_count,
                    'UNIQUE': unique_count,
                    'BOTH': in_both,
                    'ONLY_1019': in_1019,
                    'ONLY_65': in_65
                })
                if total_count > 0:
                    print(f"  {sv_type}-{region}: 总数={total_count}, 唯一={unique_count}, 两者都有={in_both}, 仅1019={in_1019}, 仅65={in_65}")
    print("分析1019数据重叠关系...")
    for sv_type in target_types:
        for region in region_types:
            region_variants = kg_df[(kg_df['SV_TYPE'] == sv_type) & (kg_df['REGION'] == region)]
            total_count = len(region_variants)
            if total_count > 0:
                in_nyuwa = 0
                in_65 = 0
                in_both = 0
                for _, variant in region_variants.iterrows():
                    variant_id = variant['ID']
                    in_nyuwa_flag = str(variant_id) in kg_1019_nw_ids
                    in_65_flag = str(variant_id) in kg_65_1019_ids
                    if in_nyuwa_flag and in_65_flag:
                        in_both += 1
                    elif in_nyuwa_flag:
                        in_nyuwa += 1
                    elif in_65_flag:
                        in_65 += 1
                unique_count = total_count - (in_both + in_nyuwa + in_65)
                kg_results.append({
                    'SV_TYPE': sv_type,
                    'REGION': region,
                    'DATASET': '1019',
                    'TOTAL_COUNT': total_count,
                    'UNIQUE': unique_count,
                    'BOTH': in_both,
                    'ONLY_1019': in_nyuwa,
                    'ONLY_65': in_65
                })
                if total_count > 0:
                    print(f"  {sv_type}-{region}: 总数={total_count}, 唯一={unique_count}, 两者都有={in_both}, 仅NyuWa={in_nyuwa}, 仅65={in_65}")
    return pd.DataFrame(nyuwa_results), pd.DataFrame(kg_results)

# ==================== 新增：计算Recall比例的函数（修正版） ====================
def compute_recall_ratios(recall_file: str, cmp_nw_file: str, cmp_1019_file: str, 
                          sv_types: List[str], region_types: List[str]) -> Tuple[Dict, Dict]:
    """
    计算每个 (sv_type, region) 分组的 recall 比例
    返回两个字典：nw_ratios 和 kg_ratios，键为 (sv_type, region)，值为比例
    """
    print("正在计算与65数据集的重叠比例...")
    # 读取 recall_file (16.sv.me.inv.info.txt)
    if not os.path.exists(recall_file):
        print(f"错误：文件不存在 {recall_file}")
        return {}, {}
    recall_df = pd.read_csv(recall_file, sep='\t', comment='#')
    # 确保列名正确
    if 'type' not in recall_df.columns or 'region' not in recall_df.columns or 'id' not in recall_df.columns:
        print("警告：recall文件列名不符合预期，尝试重新读取（无头）")
        recall_df = pd.read_csv(recall_file, sep='\t', header=None,
                                names=['chr','pos','end','id','length','type','direction','region'])
    # 读取 cmp_nw_file (65-NyuWa比较文件)，提取第三列ID
    cmp_nw_ids = set()
    if os.path.exists(cmp_nw_file):
        df_nw = pd.read_csv(cmp_nw_file, sep='\t', comment='#', header=None)
        if len(df_nw.columns) >= 3:
            cmp_nw_ids = set(df_nw[2].astype(str).str.strip())
    else:
        print(f"警告：文件不存在 {cmp_nw_file}")
    # 读取 cmp_1019_file (65-1019比较文件)，提取第十一列ID（索引10）
    cmp_1019_ids = set()
    if os.path.exists(cmp_1019_file):
        df_1019 = pd.read_csv(cmp_1019_file, sep='\t', comment='#', header=None)
        if len(df_1019.columns) >= 11:  # 确保有第11列
            cmp_1019_ids = set(df_1019[10].astype(str).str.strip())  # 第十一列索引10
        else:
            print(f"警告：文件 {cmp_1019_file} 列数不足 {len(df_1019.columns)}，期望至少11列")
    else:
        print(f"警告：文件不存在 {cmp_1019_file}")
    # 计算每个分组的比例
    nw_ratios = {}
    kg_ratios = {}
    for sv_type in sv_types:
        for region in region_types:
            # 分母：recall文件中 type==sv_type 且 region==region 的变异数
            total = len(recall_df[(recall_df['type'] == sv_type) & (recall_df['region'] == region)])
            if total == 0:
                nw_ratios[(sv_type, region)] = 0.0
                kg_ratios[(sv_type, region)] = 0.0
                continue
            # 分子：这些变异中ID出现在cmp_nw_ids中的数目
            ids_in_group = recall_df[(recall_df['type'] == sv_type) & (recall_df['region'] == region)]['id'].astype(str).str.strip()
            matched_nw = sum(1 for vid in ids_in_group if vid in cmp_nw_ids)
            matched_1019 = sum(1 for vid in ids_in_group if vid in cmp_1019_ids)
            nw_ratios[(sv_type, region)] = matched_nw / total
            kg_ratios[(sv_type, region)] = matched_1019 / total
    return nw_ratios, kg_ratios

# ==================== 修改绘图函数以添加红线 ====================
def create_combined_region_plot_with_recall(nyuwa_df: pd.DataFrame, kg_df: pd.DataFrame,
                                           nw_ratios: Dict, kg_ratios: Dict,
                                           output_file: str):
    """创建区域分布和数据集重叠的综合图表，并添加Recall曲线"""
    print("正在创建综合区域分布对比图（含Recall）...")
    try:
        plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Liberation Sans']
        plt.rcParams['axes.unicode_minus'] = False
        overlap_colors = {
            'UNIQUE': 'lightgray',
            'BOTH': 'gold',
            'ONLY_1019': 'lightcoral',
            'ONLY_65': 'lightblue'
        }
        sv_types = ['DEL', 'INS', 'DUP', 'INV', 'Alu', 'LINE', 'LTR', 'SVA']
        region_types = ['SD', 'SR', 'RM', 'Unique']
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(24, 14))
        x = np.arange(len(sv_types) * len(region_types))
        width = 0.8
        # 计算y轴最大值
        max_nw = nyuwa_df['TOTAL_COUNT'].max() if len(nyuwa_df) > 0 else 0
        max_kg = kg_df['TOTAL_COUNT'].max() if len(kg_df) > 0 else 0
        y_max = max(max_nw, max_kg) * 1.15
        y_max = max(y_max, 100)
        # --- 绘制NyuWa图 ---
        print("绘制NyuWa数据...")
        nw_data = {}
        for sv_type in sv_types:
            nw_data[sv_type] = {}
            for region in region_types:
                row = nyuwa_df[(nyuwa_df['SV_TYPE'] == sv_type) & (nyuwa_df['REGION'] == region)]
                if len(row) > 0:
                    nw_data[sv_type][region] = row.iloc[0]
                else:
                    nw_data[sv_type][region] = {'TOTAL_COUNT': 0, 'UNIQUE': 0, 'BOTH': 0, 'ONLY_1019': 0, 'ONLY_65': 0}
        # 堆叠柱状图
        for i, sv_type in enumerate(sv_types):
            for j, region in enumerate(region_types):
                idx = i * len(region_types) + j
                data = nw_data[sv_type][region]
                bottom = 0
                parts = ['BOTH', 'ONLY_1019', 'ONLY_65', 'UNIQUE']
                for part in parts:
                    value = data[part]
                    if value > 0:
                        ax1.bar(idx, value, width, bottom=bottom,
                                color=overlap_colors[part],
                                edgecolor='white', linewidth=0.5)
                        bottom += value
                total = data['TOTAL_COUNT']
                if total > 0:
                    ax1.text(idx, total + y_max * 0.01, f'{int(total):,}',
                            ha='center', va='bottom', fontsize=8, fontweight='bold')
        ax1.set_ylabel('NyuWa Variant Count', fontsize=12, fontweight='bold')
        ax1.set_ylim(0, y_max)
        ax1.set_title('NyuWa SV Type Distribution by Region and Dataset Overlap', fontsize=14, fontweight='bold')
        x_labels = []
        for sv_type in sv_types:
            for region in region_types:
                x_labels.append(f'{sv_type}\n{region}')
        ax1.set_xticks(x)
        ax1.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=9)
        ax1.grid(True, alpha=0.3, axis='y')
        # 添加Recall曲线（右轴）
        ax1_right = ax1.twinx()
        recall_values_nw = []
        for i, sv_type in enumerate(sv_types):
            for j, region in enumerate(region_types):
                ratio = nw_ratios.get((sv_type, region), 0.0)
                recall_values_nw.append(ratio)
        ax1_right.plot(x, recall_values_nw, 'r-o', linewidth=2, markersize=6, label='Recall (65 dataset)')
        ax1_right.set_ylabel('Recall', fontsize=12, fontweight='bold', color='red')
        ax1_right.set_ylim(0, 1.05)
        ax1_right.tick_params(axis='y', labelcolor='red')
        ax1_right.legend(loc='upper left', fontsize=10)
        # --- 绘制1019图 ---
        print("绘制1019数据...")
        kg_data = {}
        for sv_type in sv_types:
            kg_data[sv_type] = {}
            for region in region_types:
                row = kg_df[(kg_df['SV_TYPE'] == sv_type) & (kg_df['REGION'] == region)]
                if len(row) > 0:
                    kg_data[sv_type][region] = row.iloc[0]
                else:
                    kg_data[sv_type][region] = {'TOTAL_COUNT': 0, 'UNIQUE': 0, 'BOTH': 0, 'ONLY_1019': 0, 'ONLY_65': 0}
        for i, sv_type in enumerate(sv_types):
            for j, region in enumerate(region_types):
                idx = i * len(region_types) + j
                data = kg_data[sv_type][region]
                bottom = 0
                parts = ['BOTH', 'ONLY_1019', 'ONLY_65', 'UNIQUE']
                for part in parts:
                    value = data[part]
                    if value > 0:
                        ax2.bar(idx, value, width, bottom=bottom,
                                color=overlap_colors[part],
                                edgecolor='white', linewidth=0.5)
                        bottom += value
                total = data['TOTAL_COUNT']
                if total > 0:
                    ax2.text(idx, total + y_max * 0.01, f'{int(total):,}',
                            ha='center', va='bottom', fontsize=8, fontweight='bold')
        ax2.set_ylabel('1019 Variant Count', fontsize=12, fontweight='bold')
        ax2.set_ylim(0, y_max)
        ax2.set_title('1019 SV Type Distribution by Region and Dataset Overlap', fontsize=14, fontweight='bold')
        ax2.set_xticks(x)
        ax2.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=9)
        ax2.grid(True, alpha=0.3, axis='y')
        # 添加Recall曲线
        ax2_right = ax2.twinx()
        recall_values_kg = []
        for i, sv_type in enumerate(sv_types):
            for j, region in enumerate(region_types):
                ratio = kg_ratios.get((sv_type, region), 0.0)
                recall_values_kg.append(ratio)
        ax2_right.plot(x, recall_values_kg, 'r-o', linewidth=2, markersize=6, label='Recall (65 dataset)')
        ax2_right.set_ylabel('Recall', fontsize=12, fontweight='bold', color='red')
        ax2_right.set_ylim(0, 1.05)
        ax2_right.tick_params(axis='y', labelcolor='red')
        ax2_right.legend(loc='upper left', fontsize=10)
        # 统一图例（保留原有的堆叠图例）
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=overlap_colors['BOTH'], label='Share with both'),
            Patch(facecolor=overlap_colors['ONLY_1019'], label='Share with 1019'),
            Patch(facecolor=overlap_colors['ONLY_65'], label='Share with 65'),
            Patch(facecolor=overlap_colors['UNIQUE'], label='Unique')
        ]
        ax1.legend(handles=legend_elements, loc='upper right', fontsize=10)
        ax2.legend(handles=legend_elements, loc='upper right', fontsize=10)
        plt.suptitle('SV Type Distribution by Region and Dataset Overlap (NyuWa vs 1019 vs 65)',
                    fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(output_file, bbox_inches='tight', facecolor='white', format='pdf')
        print(f"综合区域分布对比图（含Recall）已保存到: {output_file}")
    except Exception as e:
        print(f"错误：创建综合区域分布对比图失败 - {e}")
        import traceback
        traceback.print_exc()

def save_detailed_region_statistics_with_recall(nyuwa_df: pd.DataFrame, kg_df: pd.DataFrame,
                                                nw_ratios: Dict, kg_ratios: Dict,
                                                output_file: str):
    """保存详细的区域统计结果，包含Recall信息"""
    print("正在保存详细区域统计结果（含Recall）...")
    try:
        with open(output_file, 'w') as f:
            f.write("SV_TYPE\tREGION\tDATASET\tTOTAL_COUNT\tUNIQUE\tSHARE_BOTH\tSHARE_1019\tSHARE_65\tRECALL\n")
            # NyuWa数据
            for _, row in nyuwa_df.iterrows():
                sv_type = row['SV_TYPE']
                region = row['REGION']
                recall = nw_ratios.get((sv_type, region), 0.0)
                f.write(f"{sv_type}\t{region}\t{row['DATASET']}\t{row['TOTAL_COUNT']}\t{row['UNIQUE']}\t{row['BOTH']}\t{row['ONLY_1019']}\t{row['ONLY_65']}\t{recall:.4f}\n")
            # 1019数据
            for _, row in kg_df.iterrows():
                sv_type = row['SV_TYPE']
                region = row['REGION']
                recall = kg_ratios.get((sv_type, region), 0.0)
                f.write(f"{sv_type}\t{region}\t{row['DATASET']}\t{row['TOTAL_COUNT']}\t{row['UNIQUE']}\t{row['BOTH']}\t{row['ONLY_1019']}\t{row['ONLY_65']}\t{recall:.4f}\n")
        print(f"详细区域统计结果（含Recall）已保存到: {output_file}")
    except Exception as e:
        print(f"错误：保存详细区域统计结果失败 - {e}")

# ==================== 主函数 ====================
def main():
    # 文件路径配置
    sr_t2t_file = "/Parastor300s_G30S/zhangp/pub/resources/SR_t2t.tsv.gz"
    sd_file = "/Parastor300s_G30S/haoxp/third/cmp/1019/extract/chm13v2.0_SD_merge.bed"
    rm_file = "/home3/phase2/02_tidyRes/REFERENCE/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed"
    kg_all_file = "/Parastor300s_G30S/haoxp/third/cmp/65/1019.svall.inv.addsvtype.rmsvtr.notduptoins_svlenge50.rm05mesv.rm05trsv_newsvtype.txt"
    nw_all_file = "/Parastor300s_G30S/haoxp/third/cmp/65/nyuwa.svall.addsvtype.svlenge50.rm05mesv.rm05trsv_newsvtype.txt"
    cmp_1019_nw_file = "/Parastor300s_G30S/haoxp/third/cmp/65/svcmp.1019.nw.sv.inv.over05.rmnan.addsvtype.rmdismatchINV.notduptoins_svlenge50.rm05mesv.rm05trsv.newsvtype.txt"
    cmp_65_nw_file = "/Parastor300s_G30S/haoxp/third/cmp/65/svcmp.65.nw.sv.inv.over05.rmnan.addsvtype.notduptoins_svlenge50.rm05mesv.rm05trsv.newsvtype.txt"
    cmp_65_1019_file = "/Parastor300s_G30S/haoxp/third/cmp/65/svcmp.65.1019.sv.inv.over05.rmnan.addsvtype.rmdismatchINV.notduptoins_svlenge50.rm05mesv.rm05trsv.newsvtype.txt"
    recall_file = "/Parastor300s_G30S/haoxp/third/cmp/65/16.sv.me.inv.info.txt"
    output_stat_file = "/Parastor300s_G30S/haoxp/third/cmp/65/65allsv_1019sv_1019inv_nwsv_nwme_allsvtype_region_detail_over05.notduptoins_svlenge50.rm05mesv.rm05trsv.recall.txt"
    output_plot_file = "/Parastor300s_G30S/haoxp/third/cmp/65/65allsv_1019sv_1019inv_nwsv_nwme_allsvtype_region_distribution_over05.notduptoins_svlenge50.rm05mesv.rm05trsv.recall.pdf"

    print("开始综合统计SV变异类型区域分布和数据集重叠关系（含Recall）...")
    total_start_time = time.time()

    # 1. 加载区域文件并构建区间树
    start_time = time.time()
    sr_t2t_regions_list, sr_t2t_tree = load_sr_t2t_file(sr_t2t_file)
    sd_regions_list, sd_tree = load_sd_regions(sd_file)
    sr_regions_list, sr_tree, rm_regions_list, rm_tree = load_rm_regions(rm_file)
    print(f"区域文件加载总耗时: {time.time() - start_time:.2f} 秒")

    # 2. 读取并分类所有变异
    start_time = time.time()
    kg_all_stats, kg_all_df = read_and_classify_variants(kg_all_file, "1019", sr_t2t_tree, sr_tree, sd_tree, rm_tree, sr_t2t_regions_list, sr_regions_list, sd_regions_list, rm_regions_list)
    nw_all_stats, nw_all_df = read_and_classify_variants(nw_all_file, "NyuWa", sr_t2t_tree, sr_tree, sd_tree, rm_tree, sr_t2t_regions_list, sr_regions_list, sd_regions_list, rm_regions_list)
    print(f"所有变异分类总耗时: {time.time() - start_time:.2f} 秒")

    # 3. 从比较文件中提取变异ID
    start_time = time.time()
    nw_1019_ids, kg_1019_nw_ids, nw_65_ids, kg_65_1019_ids = extract_variant_ids_from_cmp_files(
        cmp_1019_nw_file, cmp_65_nw_file, cmp_65_1019_file)
    print(f"变异ID提取总耗时: {time.time() - start_time:.2f} 秒")

    # 4. 分析重叠关系（使用已分类的区域数据）
    start_time = time.time()
    nyuwa_region_df, kg_region_df = analyze_overlap_by_region_with_classified_data(
        nw_all_df, kg_all_df, nw_1019_ids, kg_1019_nw_ids, nw_65_ids, kg_65_1019_ids)
    print(f"重叠关系分析总耗时: {time.time() - start_time:.2f} 秒")

    # 5. 计算Recall比例
    sv_types = ['DEL', 'INS', 'DUP', 'INV', 'Alu', 'LINE', 'LTR', 'SVA']
    region_types = ['SD', 'SR', 'RM', 'Unique']
    nw_ratios, kg_ratios = compute_recall_ratios(recall_file, cmp_65_nw_file, cmp_65_1019_file,
                                                 sv_types, region_types)

    # 6. 保存详细统计结果（含Recall）
    save_detailed_region_statistics_with_recall(nyuwa_region_df, kg_region_df,
                                                nw_ratios, kg_ratios, output_stat_file)

    # 7. 创建综合区域分布对比图（含Recall曲线）
    create_combined_region_plot_with_recall(nyuwa_region_df, kg_region_df,
                                           nw_ratios, kg_ratios, output_plot_file)

    total_time = time.time() - total_start_time
    print(f"处理完成！总耗时: {total_time:.2f} 秒 ({total_time/60:.2f} 分钟)")

if __name__ == "__main__":
    main()
