import pandas as pd
import numpy as np
import re
import os
import warnings

# 忽略警告
warnings.filterwarnings('ignore')

# 1. 基础工具函数
# ==========================================
file_name = 'DA-PU数据整合.xlsx' 

def clean_numeric(value, distinct_nan=False):
    """
    基础数值清洗
    distinct_nan=True 时，保留 np.nan，不转为 0.0
    """
    if pd.isna(value) or str(value).lower() in ['nan', 'none', '', 'null']:
        return np.nan if distinct_nan else 0.0
    s_val = str(value).replace('%', '')
    match = re.search(r"(\-?\d+\.?\d*)", s_val)
    if match:
        return float(match.group(1))
    return np.nan if distinct_nan else 0.0

def clean_time_to_hours(value):
    if pd.isna(value): return 0.0
    s_val = str(value).lower()
    match = re.search(r"(\-?\d+\.?\d*)", s_val)
    if not match: return 0.0
    num = float(match.group(1))
    if 'min' in s_val or '分' in s_val: return num / 60.0
    if 'day' in s_val or '天' in s_val: return num * 24.0
    return num

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
    true_keywords = ['侧链', 'side', '1']
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

# 2. 化学规则引擎
# ==========================================
def classify_role(row):
    g_type = str(row.get('group_type', '')).lower()
    name = str(row.get('component_name', '')).lower()
    mw = clean_numeric(row.get('M_i', 0))
    
    # 保持原来的强力DA识别逻辑
    da_keywords = ['furan', 'maleimide', '马来', '呋喃', 'bmi']
    if any(k in g_type or k in name for k in da_keywords): return 'DA_monomer'
        
    iso_keywords = ['isocyanate', '异氰酸酯', 'nco']
    if any(k in g_type for k in iso_keywords) or any(k in name for k in ['mdi','tdi','hdi','ipdi','ndi','cdi','hmdi']): return 'Isocyanate'
        
    if 'filler' in g_type or '填料' in g_type: return 'Other'
    if 'polyol' in g_type or '软段' in g_type: return 'Polyol_Soft'
    if 'chain_extender' in g_type or '扩链剂' in g_type: return 'Chain_Extender'
    
    if mw > 600: return 'Polyol_Soft'
    return 'Chain_Extender'

def infer_functionality(row):
    name = str(row['component_name']).lower().strip()
    trimer_kws = ['tri-hdi', 'tmn', 'nf', 'tmp', 'glycerol', 'triol', '三', 'isocyanurate']
    if any(k in name for k in trimer_kws): return 3.0
    if 'tetra' in name or 'pentaerythritol' in name or '四' in name: return 4.0
    if 'crosslinker' in name or '交联剂' in name: return 3.0
    return 2.0

# 3. 特征描述符计算
# ==========================================
def get_iso_score(row):
    if row['role'] != 'Isocyanate': return np.nan
    name = str(row['component_name']).upper()
    if 'NDI' in name or 'PPDI' in name: return 5.0
    if 'MDI' in name: return 4.0
    if 'TDI' in name: return 3.0
    if 'IPDI' in name or 'HMDI' in name: return 2.0
    if 'HDI' in name or 'HT' in name or 'CDI' in name: return 1.0
    return np.nan 

def get_hard_symmetry(row):
    if row['role'] != 'Isocyanate': return np.nan
    name = str(row['component_name']).upper()
    symmetric_list = ['MDI', 'HDI', 'NDI','CDI', 'PPDI']
    asymmetric_list = ['IPDI', 'TDI']
    if any(k in name for k in symmetric_list): return 1.0
    if any(k in name for k in asymmetric_list): return 0.0
    return 0.0 

def get_polyol_type_score(row):
    if row['role'] != 'Polyol_Soft': return np.nan
    name = str(row['component_name']).upper()
    ether_kws = ['PPG', 'PTMG', 'TMN', 'DMN', 'DL', 'PEG', 'POLYETHER']
    if any(k in name for k in ether_kws): return 0.0
    ester_kws = ['PBA', 'PCL', 'PADG', 'CAPA', 'POLYESTER', 'PCDL']
    if any(k in name for k in ester_kws): return 1.0
    return np.nan

def get_mw_score(row):
    if row['role'] != 'Polyol_Soft': return np.nan
    mw = clean_numeric(row.get('M_i', 0))
    if mw <= 1000: return 0.0
    if mw >= 3000: return 1.0
    return (mw - 1000) / 2000.0

def get_crystallinity_score(row):
    if row['role'] != 'Polyol_Soft': return np.nan
    name = str(row['component_name']).upper()
    cryst_types = ['PTMG', 'PCL', 'PBA', 'PEO', 'PEG']
    if not any(k in name for k in cryst_types): return 0.0
    mn = clean_numeric(row.get('M_i', 0))
    if mn < 10: 
        mn_match = re.search(r'(\d{3,})', name)
        mn = int(mn_match.group(1)) if mn_match else 0
    if mn < 800: return 0.0
    elif 800 <= mn < 1500: return 0.5
    else: return 1.0

# 4. 聚合计算逻辑 (修正了聚合结果提取)
# ==========================================
def calculate_features(group):
    # --- 基础物理量 ---
    group['mass_i'] = group['n_i'] * group['M_i']
    group['f_i'] = group['f_i'].fillna(2.0)
    group['moles_group'] = group['n_i'] * group['f_i']
    total_mass = group['mass_i'].sum()
    if total_mass <= 1e-9: total_mass = 1.0
    
    m_polyol = group.loc[group['role'] == 'Polyol_Soft', 'mass_i'].sum()
    m_iso = group.loc[group['role'] == 'Isocyanate', 'mass_i'].sum()
    m_extender = group.loc[group['role'] == 'Chain_Extender', 'mass_i'].sum()
    m_da = group.loc[group['role'] == 'DA_monomer', 'mass_i'].sum()
    
    X1 = m_polyol / total_mass
    
    # 修复版 DA Strategy
    da_strat = group['DA_strategy'].max() 
    
    # 修复版 X2
    if da_strat == 1:
        X2 = m_da / total_mass
    else:
        X2 = 0.0
    
    X3 = (m_iso + m_extender + m_da) / total_mass
    
    n_nco = group.loc[group['role'] == 'Isocyanate', 'moles_group'].sum()
    n_active_H = group.loc[group['role'].isin(['Polyol_Soft', 'Chain_Extender', 'DA_monomer']), 'moles_group'].sum()
    X4 = n_nco / n_active_H if n_active_H > 1e-9 else 0.0
    
    # 修复版 X5 (侧链逻辑)
    is_side_chain = group['DA_class'].max() == 1
    term_chem_cross = (group['n_i'] * (group['f_i'] - 2)).sum() / total_mass
    if term_chem_cross < 0: term_chem_cross = 0 
    term_da_side = (m_da / total_mass) if (da_strat == 1 and is_side_chain) else 0.0
    X5 = term_chem_cross + term_da_side
    
    # --- 特征加权 ---
    def weighted_avg(sub_df, score_col):
        valid = sub_df.dropna(subset=[score_col])
        if valid.empty or valid['mass_i'].sum() == 0: return np.nan 
        return np.average(valid[score_col], weights=valid['mass_i'])

    soft_df = group[group['role'] == 'Polyol_Soft']
    w_poly_type = weighted_avg(soft_df, 'raw_poly_type')
    w_mw_score = weighted_avg(soft_df, 'raw_mw_score')
    poly_crystallinity = soft_df['raw_cryst_score'].max() if not soft_df.empty else 0.0

    iso_df = group[group['role'] == 'Isocyanate']
    w_iso_score = weighted_avg(iso_df, 'raw_iso_score')
    
    if not iso_df.empty and iso_df['mass_i'].sum() > 0:
        raw_sym = np.average(iso_df['raw_hard_sym'], weights=iso_df['mass_i'])
        X_Hard_Sym = 1.0 if raw_sym > 0.5 else 0.0
    else:
        X_Hard_Sym = 0.0

    Interact_Cryst = X1 * poly_crystallinity 
    X_Synergy = X1 * poly_crystallinity * X_Hard_Sym

    # 安全提取函数
    def get_max_valid(series):
        valid = series.dropna()
        return valid.max() if not valid.empty else np.nan

    final_heal = get_max_valid(group['healing_eff'])
    final_str = get_max_valid(group['tensile_strength'])
    final_elo = get_max_valid(group['elongation'])
    final_tem = get_max_valid(group['healing_temperature'])
    final_time = get_max_valid(group['healing_time'])
    
    return pd.Series({
        'X1_SoftSeg': X1,
        'X2_DA_Content': X2,
        'X3_HardSeg': X3,
        'X4_R_Ratio': X4,
        'X5_Crosslink': X5,
        'Iso_Type': w_iso_score,
        'Hard_Symmetry': X_Hard_Sym,
        'Polyol_Type': w_poly_type,
        'Polyol_Mw_Score': w_mw_score,
        'Soft_Cryst': poly_crystallinity,
        'Interact_Cryst_Content': Interact_Cryst,
        'Synergy_Feature': X_Synergy,
        'DA_strategy': da_strat,
        'cross_class': 1 if is_side_chain else 0, 
        'poly_tem': group['poly_tem'].max(),
        'healing_temperature': final_tem,
        'healing_time': final_time,
        'strain_rate': group['strain_rate'].max(),
        'tensile_strength': final_str,
        'elongation': final_elo,
        'healing_eff': final_heal
    })

# 5. 主程序执行
# ==========================================
def main():
    base_path = os.path.dirname(os.path.abspath(__file__)) if '__file__' in locals() else os.getcwd()
    file_path = os.path.join(base_path, file_name)
    
    print(f"正在读取文件: {file_path}")
    if not os.path.exists(file_path):
        csv_path = file_path.replace('.xlsx', '.csv')
        if os.path.exists(csv_path): file_path = csv_path
    
    try:
        if file_path.endswith('.csv'): df = pd.read_csv(file_path)
        else: df = pd.read_excel(file_path)
    except Exception as e:
        print(f"读取失败: {e}"); return

    # 列名映射
    col_map = {
        '名称sample_id': 'sample_id',
        '含DA和不含(DA_strategy)': 'DA_strategy',
        '合成温度(poly_tem)': 'poly_tem',
        '交联/线性(cross_class)': 'cross_class', 
        '组分名称component_name': 'component_name',
        '组分摩尔质量M_i（g/mol）': 'M_i',
        '组分角色group_type（Isocyanate Hydroxyl Amine ': 'group_type',
        'DA位置（DA_class）': 'DA_class', 
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
    
    # --- 1. 数据填充 (核心修正点) ---
    # 【修正】：移除了 result_cols (tensile_strength, elongation, healing_eff)
    # 只有“实验条件”才需要填充，“实验结果”如果缺失就保持缺失 (NaN)
    fill_cols = ['DA_strategy', 'poly_tem', 'cross_class', 'DA_class', 
                 'strain_rate']
    
    valid_fill_cols = [c for c in fill_cols if c in df.columns]
    for col in valid_fill_cols:
        df[col] = df.groupby('sample_id')[col].ffill().bfill()

    df = df.dropna(subset=['component_name'])

    # --- 2. 标签清洗 ---
    print("1. 正在清洗标签 (DA策略与位置)...")
    if 'DA_strategy' in df.columns: 
        df['DA_strategy'] = df['DA_strategy'].apply(parse_da_strategy)
    if 'DA_class' in df.columns: 
        df['DA_class'] = df['DA_class'].apply(parse_side_chain_flag)
    else:
        df['DA_class'] = 0 
    if 'cross_class' in df.columns:
         df['cross_class'] = df['cross_class'].apply(lambda x: parse_flag_like(x, ['交联', 'cross', '1']))

    # 3. 数值清洗
    clean_cols = ['M_i', 'n_i', 'poly_tem', 'strain_rate']
    for c in clean_cols:
        if c in df.columns: df[c] = df[c].apply(clean_numeric)
    
    # 【修正】：使用 distinct_nan=True，确保 NaN 被保留，而不是转为 0
    result_cols = ['tensile_strength', 'elongation', 'healing_eff','healing_time','healing_temperature']
    for c in result_cols:
        if c in df.columns: df[c] = df[c].apply(lambda x: clean_numeric(x, distinct_nan=True))
    
    print("2. 正在识别组分角色与计算微观属性...")
    df['role'] = df.apply(classify_role, axis=1)
    df['f_i'] = df.apply(infer_functionality, axis=1)
    
    df['raw_iso_score'] = df.apply(get_iso_score, axis=1)
    df['raw_hard_sym'] = df.apply(get_hard_symmetry, axis=1)
    df['raw_poly_type'] = df.apply(get_polyol_type_score, axis=1)
    df['raw_mw_score'] = df.apply(get_mw_score, axis=1)
    df['raw_cryst_score'] = df.apply(get_crystallinity_score, axis=1)

    print("3. 正在聚合样本特征...")
    final_features = df.groupby('sample_id', sort=False).apply(calculate_features).reset_index()
    final_features = final_features.round(4)
    
    # 检查输出
    nan_heal_count = final_features['healing_eff'].isna().sum()
    print(f"\n[检查] healing_eff 为 NaN 的样本数 (预期应保留空值): {nan_heal_count}")
    print(f"[检查] DA_strategy=1 (含DA) 的样本数: {(final_features['DA_strategy'] == 1).sum()}")

    output_file = 'Final_Features_12.7.csv'
    final_features.to_csv(output_file, index=False, encoding='utf-8-sig')
    print("-" * 30)
    print(f"处理完成！文件已保存至: {output_file}")

if __name__ == "__main__":
    main()