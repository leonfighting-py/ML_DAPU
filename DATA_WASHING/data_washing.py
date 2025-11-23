import pandas as pd
import numpy as np
import re
import os

# ==========================================
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

# ==========================================
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

# --- 【关键修正 1：防止错误填充】 ---
# 仅填充 Sample ID，确保我们可以分组
df['sample_id'] = df['sample_id'].ffill()

# 定义“固有配方属性” (可以向下填充的列)
# 注意：【绝对不要】在这里填充 healing_eff，否则上一个样品的数值会覆盖下一个空样品
identity_cols = [
    'DA_strategy', 'poly_tem', 'cross_class', 
    'DA_class', 'strain_rate',
    'healing_temperature', 'healing_time' 
]
df[identity_cols] = df[identity_cols].ffill()

# 删除无效行
df = df.dropna(subset=['component_name'])

# ==========================================
# 3. 数值清洗 (差异化处理)
# ==========================================

# A. 通用清洗 (空值转0) - 用于配方量、分子量等
def clean_numeric_fill_zero(value):
    if pd.isna(value): return 0.0
    match = re.search(r"(\-?\d+\.?\d*)", str(value))
    if match: return float(match.group(1))
    return 0.0

# B. 特殊清洗 (空值保留 NaN) - 用于 healing_eff
def clean_numeric_keep_nan(value):
    # 如果本身就是空，或者字符串是 'nan'/'NaN'
    if pd.isna(value) or str(value).lower() == 'nan': 
        return np.nan
    
    # 尝试提取数字
    match = re.search(r"(\-?\d+\.?\d*)", str(value))
    if match: 
        return float(match.group(1))
    
    # 如果是无法解析的字符（且不是空），视情况而定
    # 这里假设无法解析等同于缺失
    return np.nan

def clean_time_to_hours(value):
    if pd.isna(value): return 0.0
    s_val = str(value).lower()
    match = re.search(r"(\-?\d+\.?\d*)", s_val)
    if not match: return 0.0
    num = float(match.group(1))
    if 'min' in s_val or '分' in s_val: return num / 60.0
    return num

# --- 应用清洗 ---
# 1. 配方和条件列 (填0)
cols_fill_zero = ['M_i', 'n_i', 'poly_tem', 'strain_rate', 'healing_temperature']
for col in cols_fill_zero:
    if col in df.columns: df[col] = df[col].apply(clean_numeric_fill_zero)

if 'healing_time' in df.columns:
    df['healing_time'] = df['healing_time'].apply(clean_time_to_hours)

# 2. 结果列 (拉伸强度等填0或保留NaN看需求，这里healing_eff必须保留NaN)
# 为了安全，拉伸数据如果有值就提取，没值先暂存NaN，后面聚合时max会自动忽略NaN
if 'tensile_strength' in df.columns:
    df['tensile_strength'] = df['tensile_strength'].apply(clean_numeric_keep_nan)
if 'elongation' in df.columns:
    df['elongation'] = df['elongation'].apply(clean_numeric_keep_nan)

# 3. 【重点】清洗 healing_eff (必须保留 NaN)
if 'healing_eff' in df.columns:
    df['healing_eff'] = df['healing_eff'].apply(clean_numeric_keep_nan)

# 4. 状态列清洗 (0/1)
## DA状态清洗
def parse_da_strategy(x):
    s = '' if pd.isna(x) else str(x).strip()
    s_low = s.lower()
    # 明确表示不含的情况
    if '不' in s_low and '含' in s_low:
        return 0
    # 含关键字
    if '含' in s_low:
        return 1
    # 数字形式：1/0/1.0 等
    m = re.search(r'(-?\d+\.?\d*)', s_low)
    if m:
        try:
            num = float(m.group(1))
            return 1 if num != 0 else 0
        except:
            pass
    # 默认认为不含
    return 0
df['DA_strategy'] = df['DA_strategy'].apply(parse_da_strategy)

## 交联与否清洗及DA位置清洗
def parse_flag_like(x, true_keywords=('交联','crosslink','crosslinker','1'), false_keywords=('不', '无', '0')):
    s = '' if pd.isna(x) else str(x).lower()
    for fk in false_keywords:
        if fk in s and ('交联' in s or 'crosslink' in s or '1' in s):
            return 0
    for tk in true_keywords:
        if tk in s:
            return 1
    # 尝试数值
    m = re.search(r'(-?\d+\.?\d*)', s)
    if m:
        try:
            return 1 if float(m.group(1)) != 0 else 0
        except:
            pass
    return 0
df['cross_class'] = df['cross_class'].apply(lambda x: parse_flag_like(x, true_keywords=('交联','crosslink','1')))
df['DA_class'] = df['DA_class'].apply(lambda x: parse_flag_like(x, true_keywords=('侧链','side','1')))

# ==========================================
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
    trimer_keywords = ['HT-100','tri-hdi', 'TMN', 'NF', 'TMP', 'glycerol', 'triol', '三', 'isocyanurate']
    for kw in trimer_keywords:
        if kw in name: return 3.0
    if 'tetra' in name or 'pentaerythritol' in name or '四' in name: return 4.0
    if 'crosslinker' in name or '交联剂' in name or 'crosslinker' in g_type: return 3.0
    return 2.0
df['f_i'] = df.apply(infer_functionality, axis=1)

df['mass_i'] = df['n_i'] * df['M_i']
df['moles_group'] = df['n_i'] * df['f_i']

# ==========================================
# 5. 特征聚合 (逻辑修正版)
# ==========================================
def calculate_features(group):
    # --- X1-X5 计算部分 ---
    m_polyol = group.loc[group['role'] == 'Polyol_Soft', 'mass_i'].sum()
    m_iso = group.loc[group['role'] == 'Isocyanate', 'mass_i'].sum()
    m_extender = group.loc[group['role'] == 'Chain_Extender', 'mass_i'].sum()
    m_da = group.loc[group['role'] == 'DA_monomer', 'mass_i'].sum()
    calc_total_mass = m_polyol + m_iso + m_extender + m_da
    if calc_total_mass == 0: calc_total_mass = 1e-9
    
    X1 = m_polyol / calc_total_mass
    da_strat = group['DA_strategy'].max() # 获取该配方的策略
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
    
    # --- 【核心修复逻辑】 ---
    da_strat = group['DA_strategy'].max() # 取该组的策略 (0或1)
    
    # 获取原始愈合率
    # 逻辑：取该 sample_id 组内所有非空的愈合率。
    # 如果整个组都是 NaN，则 valid_heals 为空，raw_heal 为 NaN
    valid_heals = group['healing_eff'].dropna()
    if not valid_heals.empty:
        raw_heal = valid_heals.iloc[0] # 取第一个有效值
    else:
        raw_heal = np.nan # 全是空，说明未测试
    
    # 应用规则
    if da_strat == 0:
        # 规则 1: 无 DA 策略 -> 强制为 0
        final_heal = 0.0
    else:
        # 规则 2: 含 DA 策略
        if pd.isna(raw_heal):
            # 原始为空 -> 保持 NaN (表示未测试)
            final_heal = np.nan
        else:
            # 原始有值 -> 继承原始值
            final_heal = raw_heal
    
    return pd.Series({
        'X1_SoftSeg': X1,
        'X2_DA_Content': X2,
        'X3_HardSeg': X3,
        'X4_R_Ratio': X4,
        'X5_Crosslink': X5,
        'DA_strategy': da_strat,
        'poly_tem': group['poly_tem'].max(),
        'cross_class': group['cross_class'].max(),
        'healing_temperature': group['healing_temperature'].max(),
        'healing_time': group['healing_time'].max(),
        'strain_rate': group['strain_rate'].max(),
        'tensile_strength': group['tensile_strength'].max(), # max 会自动忽略 NaN
        'elongation': group['elongation'].max(),
        'healing_eff': final_heal # <--- 修正后的值
    })

# 执行聚合
final_features = df.groupby('sample_id', sort=False).apply(calculate_features).reset_index()

# 四舍五入
round_config = {
    'X1_SoftSeg': 2, 
    'X2_DA_Content': 2,
    'X3_HardSeg': 2,
    'X4_R_Ratio': 2,
    'X5_Crosslink': 6, 
    'tensile_strength': 2,
    'elongation': 2,
    'healing_eff': 2,
    'healing_time': 2,
    'healing_temperature': 2
}
final_features = final_features.round(round_config)

# ==========================================
# 6. 输出与检查
# ==========================================
print("-" * 30)
print(f"处理完成！共提取 {len(final_features)} 条样品。")
output_file = 'Final_Features_Fixed_1121.csv'
final_features.to_csv(output_file, index=False, encoding='utf-8-sig')
print(f"\n文件已保存为: {output_file}")