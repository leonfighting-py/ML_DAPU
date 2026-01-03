import pandas as pd
import numpy as np
import re
import os
import warnings

warnings.filterwarnings('ignore')

# ======================================================
# 0. 文件名
# ======================================================
file_name = 'DA-PU数据整合.xlsx'

# ======================================================
# 1. 基础清洗函数
# ======================================================
def clean_numeric(x, distinct_nan=False):
    if pd.isna(x) or str(x).lower() in ['nan', 'none', '', 'null']:
        return np.nan if distinct_nan else 0.0
    s = str(x).replace('%', '')
    m = re.search(r"(\-?\d+\.?\d*)", s)
    return float(m.group(1)) if m else (np.nan if distinct_nan else 0.0)

def parse_da_strategy(x):
    s = '' if pd.isna(x) else str(x).lower()
    if '不' in s or 'non' in s or '0' in s:
        return 0
    return 1

def parse_side_chain_flag(x):
    s = '' if pd.isna(x) else str(x).lower()
    return 1 if ('侧链' in s or 'side' in s or s == '1') else 0

def parse_cross_class(x):
    s = '' if pd.isna(x) else str(x).lower()
    return 1 if ('交联' in s or 'cross' in s or s == '1') else 0

# ======================================================
# 2. 角色与官能度
# ======================================================
def classify_role(row):
    name = str(row['component_name']).lower()
    mw = clean_numeric(row['M_i'], True)

    if any(k in name for k in ['furan', 'maleimide', 'bmi', '呋喃', '马来']):
        return 'DA_monomer'
    if any(k in name for k in ['mdi', 'tdi', 'hdi', 'ipdi', 'ndi', 'cdi', 'hmdi']):
        return 'Isocyanate'
    if mw > 600 or 'polyol' in name:
        return 'Polyol_Soft'
    return 'Chain_Extender'

def infer_functionality(row):
    name = str(row['component_name']).lower()
    if any(k in name for k in ['tri', 'tmp', 'glycerol', '三']):
        return 3.0
    if any(k in name for k in ['tetra', 'penta', '四']):
        return 4.0
    return 2.0

# ======================================================
# 3. 微观描述符
# ======================================================
def get_iso_score(row):
    if row['role'] != 'Isocyanate':
        return np.nan
    name = row['component_name'].upper()
    if 'NDI' in name: return 5
    if 'MDI' in name: return 4
    if 'TDI' in name: return 3
    if 'IPDI' in name or 'HMDI' in name: return 2
    if 'HDI' in name or 'CDI' in name: return 1
    return np.nan

def get_hard_sym(row):
    if row['role'] != 'Isocyanate':
        return np.nan
    name = row['component_name'].upper()
    if any(k in name for k in ['MDI','HDI','NDI','CDI']):
        return 1
    return 0

def get_poly_type(row):
    if row['role'] != 'Polyol_Soft':
        return np.nan
    name = row['component_name'].upper()
    if any(k in name for k in ['PPG', 'PTMG', 'TMN', 'DMN', 'DL', 'PEG', 'POLYETHER']):
        return 0
    if any(k in name for k in ['PBA', 'PCL', 'PADG', 'CAPA', 'POLYESTER', 'PCDL']):
        return 1
    return np.nan

def get_poly_mw_score(row):
    if row['role'] != 'Polyol_Soft':
        return np.nan
    mw = clean_numeric(row['M_i'], True)
    if mw <= 1000: return 0
    if mw >= 3000: return 1
    return (mw - 1000) / 2000

def get_soft_cryst(row):
    if row['role'] != 'Polyol_Soft':
        return 0
    name = row['component_name'].upper()
    return 1 if any(k in name for k in ['PTMG','PCL','PBA','PEG']) else 0

# ======================================================
# 4. one-hot（新增）
# ======================================================
def add_onehot(df):
    name = df['component_name'].str.lower().fillna('')
    df['DA_KA'] = name.apply(lambda x: 1 if ('ka' in x or '糠胺' in x) else 0)
    df['DA_KC'] = name.apply(lambda x: 1 if ('kc' in x or '糠醇' in x) else 0)
    df['Has_Non_DA_Extender'] = name.apply(
        lambda x: 1 if any(k in x for k in ['bdo','dag','edo','hdo']) else 0
    )
    return df

# ======================================================
# 5. 样本级聚合
# ======================================================
def calculate_features(g):
    g['mass_i'] = g['n_i'] * g['M_i']
    g['f_i'] = g['f_i'].fillna(2)
    g['mol'] = g['n_i'] * g['f_i']

    total_mass = g['mass_i'].sum()
    if total_mass <= 0: total_mass = 1

    m_soft = g.loc[g['role']=='Polyol_Soft','mass_i'].sum()
    m_hard = g.loc[g['role']!='Polyol_Soft','mass_i'].sum()
    m_da = g.loc[g['role']=='DA_monomer','mass_i'].sum()

    X1 = m_soft / total_mass
    X3 = m_hard / total_mass
    X2 = m_da / total_mass if g['DA_strategy'].max()==1 else 0

    n_nco = g.loc[g['role']=='Isocyanate','mol'].sum()
    n_act = g.loc[g['role']!='Isocyanate','mol'].sum()
    X4 = n_nco / n_act if n_act>0 else 0

    X5 = max((g['n_i']*(g['f_i']-2)).sum()/total_mass,0)

    Soft_Mw = np.average(
        g.loc[g['role']=='Polyol_Soft','M_i'],
        weights=g.loc[g['role']=='Polyol_Soft','mass_i']
    ) if m_soft>0 else np.nan

    Constraint_Factor = X3 / Soft_Mw if Soft_Mw>0 else np.nan

    Soft_Cryst = g['raw_cryst'].max()
    Interact = X1 * Soft_Cryst
    Synergy = Interact * g['raw_sym'].max()

    return pd.Series({
        'X1_SoftSeg':X1,'X2_DA_Content':X2,'X3_HardSeg':X3,
        'X4_R_Ratio':X4,'X5_Crosslink':X5,
        'Iso_Type':np.nanmean(g['raw_iso']),
        'Hard_Symmetry':g['raw_sym'].max(),
        'Polyol_Type':np.nanmean(g['raw_poly']),
        'Polyol_Mw_Score':np.nanmean(g['raw_mw']),
        'Soft_Cryst':Soft_Cryst,
        'Interact_Cryst_Content':Interact,
        'Synergy_Feature':Synergy,
        'Soft_Mw':Soft_Mw,
        'Constraint_Factor':Constraint_Factor,

        'DA_KA':g['DA_KA'].max(),
        'DA_KC':g['DA_KC'].max(),
        'Has_Non_DA_Extender':g['Has_Non_DA_Extender'].max(),
        'cross_class':g['cross_class'].max(),

        'tensile_strength':g['tensile_strength'].max(),
        'elongation':g['elongation'].max(),
        'healing_eff':g['healing_eff'].max()
    })

# ======================================================
# 6. 主程序
# ======================================================
def main():
    base = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(base, file_name)
    df = pd.read_excel(path)

    col_map = {
        '名称sample_id':'sample_id',
        '含DA和不含(DA_strategy)':'DA_strategy',
        '交联/线性(cross_class)':'cross_class',
        '组分名称component_name':'component_name',
        '组分摩尔质量M_i（g/mol）':'M_i',
        '组分摩尔用量（n_i）':'n_i',
        '输出：原拉伸强度MPa（tensile_strength）':'tensile_strength',
        '输出：初始拉伸率%(elongation)':'elongation',
        '输出：自愈合率%（healing_eff）':'healing_eff'
    }

    df = df.rename(columns=col_map)
    df['sample_id'] = df['sample_id'].ffill()

    df['DA_strategy'] = df['DA_strategy'].apply(parse_da_strategy)
    df['cross_class'] = df['cross_class'].apply(parse_cross_class)

    for c in ['M_i','n_i']:
        df[c] = df[c].apply(clean_numeric)

    for c in ['tensile_strength','elongation','healing_eff']:
        df[c] = df[c].apply(lambda x: clean_numeric(x,True))

    df['role'] = df.apply(classify_role,axis=1)
    df['f_i'] = df.apply(infer_functionality,axis=1)

    df['raw_iso'] = df.apply(get_iso_score,axis=1)
    df['raw_sym'] = df.apply(get_hard_sym,axis=1)
    df['raw_poly'] = df.apply(get_poly_type,axis=1)
    df['raw_mw'] = df.apply(get_poly_mw_score,axis=1)
    df['raw_cryst'] = df.apply(get_soft_cryst,axis=1)

    df = add_onehot(df)

    final = df.groupby('sample_id',sort=False).apply(calculate_features).reset_index()
    final = final.round(4)
    final.to_csv('Final_Features_FULL_OneHot.csv',index=False,encoding='utf-8-sig')
    print('完成：Final_Features_FULL_OneHot.csv')

if __name__ == '__main__':
    main()
