import pandas as pd
import numpy as np
import re
import os

# 1. 读取数据
# ==========================================
file_name = 'DA-PU数据整合 - Sheet1.csv'
base_path = os.path.dirname(os.path.abspath(__file__)) if '__file__' in locals() else os.getcwd()
file_path = os.path.join(base_path, file_name)

print(f"正在读取文件: {file_path}")

try:
    df = pd.read_excel(file_path)
except:
    print("尝试读取 CSV...")
    df = pd.read_csv(file_path.replace('.xlsx', '.csv'))

# 2. 列名映射
# ==========================================
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
    '自愈合温度（healing_temperature K)': 'healing_temperature', 
    '自愈合时间（healing time h）': 'healing_time',
    '输出：原拉伸强度MPa（tensile_strength）': 'tensile_strength',
    '输出：初始拉伸率%(elongation)': 'elongation',
    '输出：自愈合率%（healing_eff）': 'healing_eff'
}

# 清理列名空格
df.columns = [str(c).strip() for c in df.columns]
col_map_stripped = {k.strip(): v for k, v in col_map.items()}
df = df.rename(columns=col_map_stripped)

# --- 基础填充 ---
df['sample_id'] = df['sample_id'].ffill()
identity_cols = [
    'DA_strategy', 'poly_tem', 'cross_class', 
    'DA_class', 'strain_rate',
    'healing_temperature', 'healing_time' 
]
df[identity_cols] = df[identity_cols].ffill()
df = df.dropna(subset=['component_name'])

# 3. 数值清洗
# ==========================================
def clean_numeric_fill_zero(value):
    if pd.isna(value): return 0.0
    match = re.search(r"(\-?\d+\.?\d*)", str(value))
    if match: return float(match.group(1))
    return 0.0

def clean_numeric_keep_nan(value):
    if pd.isna(value) or str(value).lower() == 'nan': return np.nan
    match = re.search(r"(\-?\d+\.?\d*)", str(value))
    if match: return float(match.group(1))
    return np.nan

def clean_time_to_hours(value):
    if pd.isna(value): return 0.0
    s_val = str(value).lower()
    match = re.search(r"(\-?\d+\.?\d*)", s_val)
    if not match: return 0.0
    num = float(match.group(1))
    if 'min' in s_val or '分' in s_val: return num / 60.0
    return num

# 应用清洗
cols_fill_zero = ['M_i', 'n_i', 'poly_tem', 'strain_rate', 'healing_temperature']
for col in cols_fill_zero:
    if col in df.columns: df[col] = df[col].apply(clean_numeric_fill_zero)

if 'healing_time' in df.columns:
    df['healing_time'] = df['healing_time'].apply(clean_time_to_hours)

if 'tensile_strength' in df.columns:
    df['tensile_strength'] = df['tensile_strength'].apply(clean_numeric_keep_nan)
if 'elongation' in df.columns:
    df['elongation'] = df['elongation'].apply(clean_numeric_keep_nan)
if 'healing_eff' in df.columns:
    df['healing_eff'] = df['healing_eff'].apply(clean_numeric_keep_nan)

# 状态列清洗
def parse_da_strategy(x):
    s = '' if pd.isna(x) else str(x).strip().lower()
    if '不' in s and '含' in s: return 0
    if '含' in s: return 1
    m = re.search(r'(-?\d+\.?\d*)', s)
    if m: return 1 if float(m.group(1)) != 0 else 0
    return 0
df['DA_strategy'] = df['DA_strategy'].apply(parse_da_strategy)

def parse_flag_like(x, true_keywords=('交联','crosslink','1')):
    s = '' if pd.isna(x) else str(x).lower()
    if '不' in s or '0' in s: return 0 # 简单处理，若包含否定词则为0
    for tk in true_keywords:
        if tk in s: return 1
    return 0
df['cross_class'] = df['cross_class'].apply(lambda x: parse_flag_like(x, true_keywords=('交联','crosslink','1')))
df['DA_class'] = df['DA_class'].apply(lambda x: parse_flag_like(x, true_keywords=('侧链','side','1')))

# 4. 角色与官能度
# ==========================================
def classify_role(row):
    g_type = str(row.get('group_type', '')).lower()
    name = str(row['component_name']).lower()
    mw = row['M_i']
    if 'furan' in g_type or 'maleimide' in g_type or '马来' in g_type or '呋喃' in g_type: return 'DA_monomer'
    if 'da' in name or 'furan' in name or 'bmi' in name: return 'DA_monomer'
    if 'isocyanate' in g_type or '异氰酸酯' in g_type or 'nco' in g_type: return 'Isocyanate'
    if 'filler' in g_type or '填料' in g_type: return 'Other'
    if 'polyol' in g_type or '软段' in g_type: return 'Polyol_Soft'
    if 'chain_extender' in g_type or '扩链剂' in g_type: return 'Chain_Extender'
    if mw > 600: return 'Polyol_Soft'
    return 'Chain_Extender'
df['role'] = df.apply(classify_role, axis=1)

def infer_functionality(row):
    name = str(row['component_name']).lower().strip()
    g_type = str(row.get('group_type', '')).lower()
    trimer_keywords = ['tri-hdi', 'tmn', 'nf', 'tmp', 'glycerol', 'triol', '三', 'isocyanurate']
    for kw in trimer_keywords:
        if kw in name: return 3.0
    if 'tetra' in name or 'pentaerythritol' in name or '四' in name: return 4.0
    if 'crosslinker' in name or '交联剂' in name or 'crosslinker' in g_type: return 3.0
    return 2.0
df['f_i'] = df.apply(infer_functionality, axis=1)

df['mass_i'] = df['n_i'] * df['M_i']
df['moles_group'] = df['n_i'] * df['f_i']

# 5. 特征识别 (Iso_rank, Polyol_type, Mw_score)
# ==========================================

## 1. 异氰酸酯种类+硬度
def get_iso_score(row):
    if row['role'] != 'Isocyanate': return np.nan # 非硬段不参与打分
    name = str(row['component_name']).upper()
    
    # 优先级：NDI/PPDI > MDI > TDI > IPDI/HMDI > HDI
    if 'NDI' in name or 'PPDI' in name: return 5
    if 'MDI' in name: return 4
    if 'TDI' in name: return 3
    if 'IPDI' in name or 'HMDI' in name: return 2
    if 'HDI' in name or 'HT-100' in name: return 1
    return np.nan # 未知异氰酸酯，不参与加权

# 2. 软段类型 (0=聚醚/软, 1=聚酯/硬)
def get_polyol_type_score(row):
    if row['role'] != 'Polyol_Soft': return np.nan
    name = str(row['component_name']).upper()
    
    # 聚醚类 (Soft -> 0)
    ether_kws = ['PPG', 'PTMG', 'TMN', 'DMN', 'DL']
    if any(kw in name for kw in ether_kws): return 0.0
    
    # 聚酯类 (Hard -> 1)
    ester_kws = ['PBA', 'PCL', 'PADG']
    if any(kw in name for kw in ester_kws): return 1.0
    
    return np.nan

# 3. 软段分子量 (<=1000 -> 0, >=2000 -> 1, 中间线性插值)
def get_mw_score(row):
    if row['role'] != 'Polyol_Soft': return np.nan
    mw = row['M_i']
    if mw <= 1000: return 0.0
    if mw >= 2000: return 1.0
    # 1000-2000 之间记为 0.5 (代表中等相分离)
    return 0.5

# 应用到每一行，暂存为辅助列
df['raw_iso_score'] = df.apply(get_iso_score, axis=1)
df['raw_poly_type'] = df.apply(get_polyol_type_score, axis=1)
df['raw_mw_score'] = df.apply(get_mw_score, axis=1)

# 6. 特征聚合
# ==========================================
def calculate_features(group):
    # --- 基础质量计算 ---
    m_polyol = group.loc[group['role'] == 'Polyol_Soft', 'mass_i'].sum()
    m_iso = group.loc[group['role'] == 'Isocyanate', 'mass_i'].sum()
    m_extender = group.loc[group['role'] == 'Chain_Extender', 'mass_i'].sum()
    m_da = group.loc[group['role'] == 'DA_monomer', 'mass_i'].sum()
    
    calc_total_mass = m_polyol + m_iso + m_extender + m_da
    if calc_total_mass == 0: calc_total_mass = 1e-9
    
    X1 = m_polyol / calc_total_mass
    
    da_strat = group['DA_strategy'].max()
    if da_strat == 0:
        X2 = 0.0  
    else:
        X2 = m_da / calc_total_mass 
        
    X3 = (m_iso + m_da + m_extender) / calc_total_mass
    
    n_nco = group.loc[group['role'] == 'Isocyanate', 'moles_group'].sum()
    n_active_H = group.loc[group['role'].isin(['Polyol_Soft', 'Chain_Extender', 'DA_monomer']), 'moles_group'].sum()
    X4 = n_nco / n_active_H if n_active_H > 0 else 0.0
        
    term1 = (group['n_i'] * (group['f_i'] - 2)).sum() / calc_total_mass
    is_side_chain = group['DA_class'].max() == 1 
    term2 = (group.loc[group['role'] == 'DA_monomer', 'n_i'].sum() / calc_total_mass) if is_side_chain else 0.0
    X5 = term1 + term2

    # --- Iso Type 加权计算 ---
    iso_rows = group.loc[(group['role'] == 'Isocyanate') & (group['raw_iso_score'].notna())]
    if not iso_rows.empty and iso_rows['mass_i'].sum() > 0:
        # 加权平均：(m1*s1 + m2*s2) / (m1+m2)
        weighted_iso = np.average(iso_rows['raw_iso_score'], weights=iso_rows['mass_i'])
    else:
        weighted_iso = 0 # 没有识别出异氰酸酯类型
        
    # --- Polyol Type 加权计算 ---
    poly_rows = group.loc[(group['role'] == 'Polyol_Soft') & (group['raw_poly_type'].notna())]
    if not poly_rows.empty and poly_rows['mass_i'].sum() > 0:
        # 结果接近0表示全是聚醚，接近1表示全是聚酯，0.5表示混用
        weighted_poly_type = np.average(poly_rows['raw_poly_type'], weights=poly_rows['mass_i'])
    else:
        weighted_poly_type = 0 # 默认为0或标记缺失
        
    # --- Mw Score 加权计算 ---
    mw_rows = group.loc[(group['role'] == 'Polyol_Soft') & (group['raw_mw_score'].notna())]
    if not mw_rows.empty and mw_rows['mass_i'].sum() > 0:
        weighted_mw_score = np.average(mw_rows['raw_mw_score'], weights=mw_rows['mass_i'])
    else:
        weighted_mw_score = 0

    # --- 结果提取 ---
    valid_heals = group['healing_eff'].dropna()
    if not valid_heals.empty:
        raw_heal = valid_heals.iloc[0]
    else:
        raw_heal = np.nan
    
    if da_strat == 0:
        final_heal = 0.0
    else:
        final_heal = raw_heal if not pd.isna(raw_heal) else np.nan

    return pd.Series({
        'X1_SoftSeg': X1,
        'X2_DA_Content': X2,
        'X3_HardSeg': X3,
        'X4_R_Ratio': X4,
        'X5_Crosslink': X5,
        'Iso_Type': weighted_iso,  
        'Polyol_Type': weighted_poly_type, #
        'Polyol_Mw_Score': weighted_mw_score,   
        'DA_strategy': da_strat,
        'poly_tem': group['poly_tem'].max(),
        'cross_class': group['cross_class'].max(),
        'healing_temperature': group['healing_temperature'].max(),
        'healing_time': group['healing_time'].max(),
        'strain_rate': group['strain_rate'].max(),
        'tensile_strength': group['tensile_strength'].max(),
        'elongation': group['elongation'].max(),
        'healing_eff': final_heal
    })

# 执行聚合
final_features = df.groupby('sample_id', sort=False).apply(calculate_features).reset_index()

# 四舍五入 (新特征保留2位小数，因为加权平均后可能是小数，这比整数包含更多信息)
round_config = {
    'X1_SoftSeg': 2, 
    'X2_DA_Content': 2,
    'X3_HardSeg': 2,
    'X4_R_Ratio': 2,
    'X5_Crosslink': 6, 
    'Iso_Type': 2, 
    'Polyol_Type': 2,
    'Polyol_Mw_Score': 2,
    'tensile_strength': 2,
    'elongation': 2,
    'healing_eff': 2,
    'healing_time': 2,
    'healing_temperature': 2
}
final_features = final_features.round(round_config)

# 7. 输出
# ==========================================
print("-" * 30)
print(f"处理完成！共提取 {len(final_features)} 条样品。")
output_file = 'Final_Features_NewRules_1121.csv'
final_features.to_csv(output_file, index=False, encoding='utf-8-sig')
print(f"\n文件已保存为: {output_file}，请核对")