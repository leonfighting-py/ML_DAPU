import pandas as pd
import numpy as np
import re
import os
import warnings

# 忽略警告
warnings.filterwarnings('ignore')

# ==========================================
# 1. 基础工具函数 (保持不变)
# ==========================================
file_path = r"C:\Users\leon\OneDrive\Desktop\Project_DAPU\DATA_WASHING"
file_name = 'DA-PU数据整合.xlsx'

def clean_numeric(value, distinct_nan=False):
    if pd.isna(value) or str(value).lower() in ['nan', 'none', '', 'null']:
        return np.nan if distinct_nan else 0.0
    s_val = str(value).replace('%', '')
    match = re.search(r"(\-?\d+\.?\d*)", s_val)
    if match:
        return float(match.group(1))
    return np.nan if distinct_nan else 0.0

def parse_da_strategy(x):
    s = '' if pd.isna(x) else str(x).strip().lower()
    if '不' in s and '含' in s: return 0
    if 'non' in s and 'contain' in s: return 0
    if '含' in s or 'contain' in s: return 1
    m = re.search(r'(-?\d+\.?\d*)', s)
    if m: return 1 if float(m.group(1)) != 0 else 0
    return 0

def parse_side_chain_flag(x):
    s = '' if pd.isna(x) else str(x).lower()
    true_keywords = ['侧链', 'side', '1', 'pendant']
    for tk in true_keywords:
        if tk in s: return 1
    return 0

def parse_flag_like(x, true_keywords=('交联', 'crosslink', '1')):
    s = '' if pd.isna(x) else str(x).lower()
    if '不' in s or 'non' in s or '0' in s: return 0
    extended_keywords = list(true_keywords) + ['yes', 'y', '有', 'true']
    for tk in extended_keywords:
        if tk in s: return 1
    return 0

# ==========================================
# 2. 化学分类与评分引擎 (核心升级)
# ==========================================
def classify_role(row):
    """基础角色分类"""
    g_type = str(row.get('group_type', '')).lower()
    name = str(row.get('component_name', '')).lower()
    mw = clean_numeric(row.get('M_i', 0))
    
    # DA单体识别
    da_keywords = ['furan', 'maleimide', '马来', '呋喃', 'bmi', 'ka', 'kc']
    if any(k in g_type or k in name for k in da_keywords): return 'DA_monomer'
        
    # 异氰酸酯识别
    iso_keywords = ['isocyanate', '异氰酸酯', 'nco']
    if any(k in g_type for k in iso_keywords) or any(k in name for k in ['mdi','tdi','hdi','ipdi','ndi','cdi','hmdi']): return 'Isocyanate'
        
    # 软段识别
    if 'polyol' in g_type or '软段' in g_type: return 'Polyol_Soft'
    if mw > 600: return 'Polyol_Soft'
    
    # 剩下的都是扩链剂/交联剂
    return 'Chain_Extender'

def infer_functionality(row):
    """推断官能度"""
    name = str(row['component_name']).lower().strip()
    # 三官能度关键词
    trimer_kws = ['tri', 'tmp', 'glycerol', 'triol', '三', 'isocyanurate', 'crosslinker', '交联剂']
    if any(k in name for k in trimer_kws): return 3.0
    # 四官能度关键词
    if 'tetra' in name or 'pentaerythritol' in name or '四' in name: return 4.0
    return 2.0

def get_da_linker_type(row):
    """DA类型：1=胺(KA), 2=醇(KC), 0=无"""
    if row['role'] != 'DA_monomer': return 0
    name = str(row['component_name']).lower()
    # 胺类关键词
    if any(k in name for k in ['amine', 'ka', '胺']): return 1
    # 醇类关键词
    if any(k in name for k in ['alcohol', 'kc', '醇']): return 2
    # 默认归为胺(反应活性高)或根据实际情况调整
    return 1 

def get_iso_class(row):
    """硬段类型：0=脂肪族(IPDI), 1=芳香族(MDI)"""
    if row['role'] != 'Isocyanate': return np.nan
    name = str(row['component_name']).upper()
    aromatic = ['MDI', 'TDI', 'NDI', 'PPDI', 'XDI']
    if any(k in name for k in aromatic): return 1 # 芳香族
    return 0 # 脂肪族 (IPDI, HDI, HMDI)

def get_polyol_class(row):
    """软段类型：0=聚醚, 1=聚酯"""
    if row['role'] != 'Polyol_Soft': return np.nan
    name = str(row['component_name']).upper()
    ester_kws = ['PBA', 'PCL', 'PADG', 'CAPA', 'POLYESTER', 'PCDL', 'PLA', '酯']
    if any(k in name for k in ester_kws): return 1 # 聚酯
    return 0 # 聚醚 (PPG, PTMG, PEG...)

def get_hard_symmetry(row):
    """硬段对称性：0=不对称, 1=对称"""
    if row['role'] != 'Isocyanate': return np.nan
    name = str(row['component_name']).upper()
    # 对称分子
    symmetric = ['MDI', 'HDI', 'NDI', 'CDI', 'PPDI', 'CHDI']
    if any(k in name for k in symmetric): return 1.0
    return 0.0 # 不对称 (IPDI, TDI)

def get_extender_class(row):
    # 只处理扩链剂角色
    if row['role'] != 'Chain_Extender': 
        return 0.0 # 非扩链剂默认为 0，不影响特征计算
    
    name = str(row['component_name']).upper().strip()
    
    # === 1. 刚性扩链剂关键词 (Rigid) ===
    # BQDO: 苯醌结构 (环状)
    # HQEE/HER: 苯环结构
    # CHDM/ISOSORBIDE: 脂环/杂环结构
    rigid_kws = [
        'BQDO', 'HQEE', 'HER', 'MOCA', 'CHDM', 
        'ISOSORBIDE', '环', 'BENZ', 'CYCLO', 'PHE'
    ]
    
    if any(k in name for k in rigid_kws):
        return 1.0
        
    # === 2. 柔性扩链剂 (Flexible) ===
    # DAG (二甘醇) 在这里，默认返回 0
    # BDO, HDO, EG 也在这里
    return 0.0

def get_soft_cryst_flag(row):
    """软段结晶性：0=否, 1=是"""
    if row['role'] != 'Polyol_Soft': return 0.0
    name = str(row['component_name']).upper()
    # 可结晶软段列表
    cryst_types = ['PTMG', 'PCL', 'PBA', 'PEO', 'PEG', 'PLLA', 'PCDL']
    if any(k in name for k in cryst_types):
        # 简单分子量判断：太小的不结晶
        mw = clean_numeric(row.get('M_i', 0))
        if mw > 0 and mw < 600: return 0.0
        return 1.0
    return 0.0

# ==========================================
# 3. 聚合计算逻辑 (Feature Engineering Core)
# ==========================================
def calculate_features(group):
    # 1. 基础物理量计算
    group['mass_i'] = group['n_i'] * group['M_i']
    group['f_i'] = group['f_i'].fillna(2.0)
    
    total_mass = group['mass_i'].sum()
    if total_mass <= 1e-9: total_mass = 1.0

    # 2. 组分质量分类
    # X1: 软段
    m_polyol = group.loc[group['role'] == 'Polyol_Soft', 'mass_i'].sum()
    X1 = m_polyol / total_mass
    
    # X2: DA单体 (强制识别)
    m_da = group.loc[group['role'] == 'DA_monomer', 'mass_i'].sum()
    X2 = m_da / total_mass

    # X5: 额外非DA交联剂 (关键修正)
    # 逻辑：角色是 Chain_Extender 且 官能度 > 2 的组分
    cross_mask = (group['role'] == 'Chain_Extender') & (group['f_i'] > 2.0)
    m_add_cross = group.loc[cross_mask, 'mass_i'].sum()
    X5 = m_add_cross / total_mass
    
    # X3: 硬段 (异氰酸酯 + 线性扩链剂)
    # 逻辑：剩下的都是硬段 (或者用 1 - X1 - X2 - X5)
    m_iso = group.loc[group['role'] == 'Isocyanate', 'mass_i'].sum()
    # 线性扩链剂 = 扩链剂角色 且 官能度 <= 2
    linear_ext_mask = (group['role'] == 'Chain_Extender') & (group['f_i'] <= 2.0)
    m_linear_ext = group.loc[linear_ext_mask, 'mass_i'].sum()
    
    X3 = (m_iso + m_linear_ext) / total_mass
    
    # X4: R值 (NCO/OH)
    # 计算 NCO 摩尔数
    n_nco = (group.loc[group['role'] == 'Isocyanate', 'n_i'] * group.loc[group['role'] == 'Isocyanate', 'f_i']).sum()
    # 计算 活性氢 摩尔数 (除去异氰酸酯的所有组分)
    # 注意：这里假设 DA 单体也带活性氢参与反应
    n_active = (group.loc[group['role'] != 'Isocyanate', 'n_i'] * group.loc[group['role'] != 'Isocyanate', 'f_i']).sum()
    X4 = n_nco / n_active if n_active > 0 else 0.0

    # 3. 身份特征提取 (Identity Layer)
    # DA类型 (取最大值，优先识别出有的)
    da_type_series = group.apply(get_da_linker_type, axis=1)
    DA_Linker_Type = da_type_series.max() # 0, 1(KA), 2(KC)
    
    # 网络拓扑 (0=主链, 1=侧链)
    # 优先读取 DA_class 标签，如果没有则默认为主链(0)
    Network_Topology = 1 if group['DA_class'].max() == 1 else 0
    
    # 是否有额外交联剂
    Add_Crosslinker_Flag = 1 if X5 > 0 else 0
    
    # 软段/硬段大类 (加权平均后取整，或取众数)
    # 这里简化：取含量最大的那个的类型
    poly_idx = group.loc[group['role'] == 'Polyol_Soft', 'mass_i'].idxmax() if m_polyol > 0 else None
    Polyol_Class = get_polyol_class(group.loc[poly_idx]) if poly_idx is not None else 0
    
    iso_idx = group.loc[group['role'] == 'Isocyanate', 'mass_i'].idxmax() if m_iso > 0 else None
    Iso_Class = get_iso_class(group.loc[iso_idx]) if iso_idx is not None else 0

    ext_scores = group.apply(get_extender_class, axis=1)
    Extender_Class = 1 if ext_scores.max() > 0 else 0

    # 4. 机理特征提取 (Mechanism Layer)
    # 软段分子量 Soft_Mw
    if m_polyol > 0:
        Soft_Mw = np.average(group.loc[group['role']=='Polyol_Soft', 'M_i'], 
                             weights=group.loc[group['role']=='Polyol_Soft', 'mass_i'])
    else:
        Soft_Mw = 0.0 # 或者 np.nan, 但建议0以避免报错
        
    # 硬段对称性 Hard_Symmetry
    if m_iso > 0:
        sym_scores = group.apply(get_hard_symmetry, axis=1)
        # 加权平均
        w_sym = np.average(sym_scores[group['role']=='Isocyanate'], 
                           weights=group.loc[group['role']=='Isocyanate', 'mass_i'])
        Hard_Symmetry = 1 if w_sym > 0.5 else 0
    else:
        Hard_Symmetry = 0
        
    # 软段结晶性 Soft_Cryst
    if m_polyol > 0:
        cryst_scores = group.apply(get_soft_cryst_flag, axis=1)
        # 只要有能结晶的成分且含量不低，就视为能结晶
        Soft_Cryst = 1 if cryst_scores.max() > 0 else 0
    else:
        Soft_Cryst = 0
        
    # 协同因子 Synergy_Feature
    Synergy_Feature = X1 * Soft_Cryst * Hard_Symmetry
    
    # 冻结因子 Constraint_Factor
    # 防止分母为0
    Constraint_Factor = X3 / Soft_Mw if Soft_Mw > 100 else 0.0 # 假设Mw极小视为完全冻结或无软段

    # 5. 结果提取 (保持原样)
    def get_val(col):
        v = group[col].dropna()
        return v.iloc[0] if not v.empty else np.nan

    return pd.Series({
        # --- 用量层 ---
        'X1_Soft_Content': X1,
        'X2_DA_Content': X2,
        'X3_Hard_Content': X3,
        'X4_R_Ratio': X4,
        'X5_Add_Crosslink': X5,

        # --- 身份层 ---
        'DA_Linker_Type': DA_Linker_Type,
        'Network_Topology': Network_Topology,
        'Add_Crosslinker_Flag': Add_Crosslinker_Flag,
        'Polyol_Class': Polyol_Class,
        'Iso_Class': Iso_Class,
        'Extender_Class': Extender_Class,

        # --- 机理层 ---
        'Soft_Mw': Soft_Mw,
        'Constraint_Factor': Constraint_Factor,
        'Hard_Symmetry': Hard_Symmetry,
        'Soft_Cryst': Soft_Cryst,
        'Synergy_Feature': Synergy_Feature,
        
        # --- 实验条件与结果 (保留) ---
        'poly_tem': get_val('poly_tem'),
        'strain_rate': get_val('strain_rate'),
        'healing_temperature': get_val('healing_temperature'),
        'healing_time': get_val('healing_time'),
        'tensile_strength': get_val('tensile_strength'),
        'elongation': get_val('elongation'),
        'healing_eff': get_val('healing_eff')
    })

# ==========================================
# 4. 主程序执行
# ==========================================
def main():
    file_path = r"C:\Users\leon\OneDrive\Desktop\Project_DAPU\DATA_WASHING\DA-PU数据整合.xlsx"
    
    print(f"正在处理文件: {file_name}")
    try:
        if file_path.endswith('.csv'): df = pd.read_csv(file_path)
        else: df = pd.read_excel(file_path)
    except Exception as e:
        print(f"读取失败，请检查文件路径: {e}")
        return

    # 1. 预处理列名
    col_map = {
        '名称sample_id': 'sample_id',
        '含DA和不含(DA_strategy)': 'DA_strategy', # 仅用于辅助判断
        '合成温度(poly_tem)': 'poly_tem',
        '交联/线性(cross_class)': 'cross_class', # 仅用于辅助
        '组分名称component_name': 'component_name',
        '组分摩尔质量M_i（g/mol）': 'M_i',
        '组分角色group_type（Isocyanate Hydroxyl Amine ': 'group_type',
        'DA位置（DA_class）': 'DA_class', # 0=主链, 1=侧链
        '组分摩尔用量（n_i）': 'n_i',
        '拉伸测试条件(strain_rate mm/min)': 'strain_rate',
        '自愈合温度（healing_temperature ℃)': 'healing_temperature', 
        '自愈合时间（healing time h）': 'healing_time',
        '输出：原拉伸强度MPa（tensile_strength）': 'tensile_strength',
        '输出：初始拉伸率%(elongation)': 'elongation',
        '输出：自愈合率%（healing_eff）': 'healing_eff'
    }
    
    df.columns = [str(c).strip() for c in df.columns]
    df = df.rename(columns={k.strip(): v for k, v in col_map.items()})
    df['sample_id'] = df['sample_id'].ffill()
    
    # 填充实验条件 (不仅是bfill，还要ffill确保组内一致)
    cond_cols = ['DA_strategy', 'DA_class', 'poly_tem', 'strain_rate', 
                 'healing_temperature', 'healing_time', 
                 'tensile_strength', 'elongation', 'healing_eff']
    for c in cond_cols:
        if c in df.columns:
            df[c] = df.groupby('sample_id')[c].transform(lambda x: x.ffill().bfill())

    df = df.dropna(subset=['component_name'])

    # 2. 清洗数值列
    num_cols = ['M_i', 'n_i', 'poly_tem', 'strain_rate', 
                'healing_temperature', 'healing_time',
                'tensile_strength', 'elongation', 'healing_eff']
    for c in num_cols:
        if c in df.columns:
            # 结果列保留NaN，条件列转0 (视情况而定，这里保留NaN更安全，但在计算时fillna)
            distinct = True if c in ['tensile_strength', 'elongation', 'healing_eff'] else False
            df[c] = df[c].apply(lambda x: clean_numeric(x, distinct_nan=distinct))

    # 3. 辅助标签清洗
    if 'DA_class' in df.columns:
        df['DA_class'] = df['DA_class'].apply(parse_side_chain_flag)
    else:
        df['DA_class'] = 0

    # 4. 角色识别
    print("正在识别化学组分角色...")
    df['role'] = df.apply(classify_role, axis=1)
    df['f_i'] = df.apply(infer_functionality, axis=1)

    # 5. 聚合计算
    print("正在计算全能型特征仓库...")
    final_df = df.groupby('sample_id', sort=False).apply(calculate_features).reset_index()
    final_df = final_df.round(4)
    
    # 6. 保存
    out_file = 'Final_Features_Ultimate.csv'
    final_df.to_csv(out_file, index=False, encoding='utf-8-sig')
    
    print("-" * 30)
    print(f"特征工程完成！")
    print(f"生成的特征数: {final_df.shape[1]}")
    print(f"文件已保存至: {out_file}")
    print(f"包含 DA_Linker_Type 分布:\n{final_df['DA_Linker_Type'].value_counts()}")

if __name__ == "__main__":
    main()