import pandas as pd
import numpy as np
import re
import warnings
from pathlib import Path


warnings.filterwarnings('ignore')


file_name = 'DA-PU数据整合-new_furan_cleaned.csv'
BASE_DIR = Path(__file__).resolve().parent
WORKSPACE_DIR = BASE_DIR.parent
INPUT_FILE = BASE_DIR / file_name
OUTPUT_DIR = WORKSPACE_DIR / 'forecast_ML_3.8'
OUTPUT_FILE = OUTPUT_DIR / 'Final_Features_3.8.csv'


FEATURE_NAME_MAP = {
    'X1_Soft_Content': 'Soft_Segment_Content',
    'X2_DA_Content': ('Furan_Content', 'Maleimide_Content'),
    'X3_Hard_Content': 'Hard_Segment_Content',
    'X4_R_Ratio': 'NCO_OH_Ratio',
    'X5_Add_Crosslink': 'Extra_Crosslinker_Content',
}


ORIGINAL_14_FEATURES = [
    'X1_Soft_Content',
    'X2_DA_Content',
    'X3_Hard_Content',
    'X4_R_Ratio',
    'X5_Add_Crosslink',
    'Network_Topology',
    'Polyol_Class',
    'Iso_Class',
    'Soft_Mw',
    'Constraint_Factor',
    'Hard_Symmetry',
    'Soft_Cryst',
    'Synergy_Feature',
]


RECOMMENDED_FEATURES_V38 = [
    'Soft_Segment_Content',
    'Furan_Content',
    'Maleimide_Content',
    'Hard_Segment_Content',
    'NCO_OH_Ratio',
    'Extra_Crosslinker_Content',
    'Furan_Type',
    'Network_Topology',
    'Polyol_Class',
    'Iso_Class',
    'Soft_Mw',
    'Constraint_Factor',
    'Hard_Symmetry',
    'Soft_Cryst',
    'Synergy_Feature',
]


def clean_numeric(value, distinct_nan=False):
    if pd.isna(value) or str(value).lower() in ['nan', 'none', '', 'null']:
        return np.nan if distinct_nan else 0.0
    s_val = str(value).replace('%', '')
    match = re.search(r'(\-?\d+\.?\d*)', s_val)
    if match:
        return float(match.group(1))
    return np.nan if distinct_nan else 0.0


def parse_da_strategy(x):
    s = '' if pd.isna(x) else str(x).strip().lower()
    if '不' in s and '含' in s:
        return 0
    if 'non' in s and 'contain' in s:
        return 0
    if '含' in s or 'contain' in s:
        return 1
    match = re.search(r'(-?\d+\.?\d*)', s)
    if match:
        return 1 if float(match.group(1)) != 0 else 0
    return 0


def parse_side_chain_flag(x):
    s = '' if pd.isna(x) else str(x).lower()
    true_keywords = ['侧链', 'side', '1', 'pendant']
    for keyword in true_keywords:
        if keyword in s:
            return 1
    return 0


def classify_role(row):
    g_type = str(row.get('group_type', '')).lower()
    name = str(row.get('component_name', '')).lower()
    mw = clean_numeric(row.get('M_i', 0))

    da_keywords = ['furan', 'maleimide', '马来', '呋喃', 'bmi', 'ka', 'kc']
    if any(keyword in g_type or keyword in name for keyword in da_keywords):
        return 'DA_monomer'

    iso_keywords = ['isocyanate', '异氰酸酯', 'nco']
    if any(keyword in g_type for keyword in iso_keywords) or any(keyword in name for keyword in ['mdi', 'tdi', 'hdi', 'ipdi', 'ndi', 'cdi', 'hmdi']):
        return 'Isocyanate'

    if 'polyol' in g_type or '软段' in g_type:
        return 'Polyol_Soft'
    if mw > 600:
        return 'Polyol_Soft'

    return 'Chain_Extender'


def infer_functionality(row):
    name = str(row['component_name']).lower().strip()
    trimer_keywords = ['tri', 'tmp', 'glycerol', 'triol', '三', 'isocyanurate', 'crosslinker', '交联剂']
    if any(keyword in name for keyword in trimer_keywords):
        return 3.0
    if 'tetra' in name or 'pentaerythritol' in name or '四' in name:
        return 4.0
    return 2.0


def normalize_furan_type_label(value):
    label = '' if pd.isna(value) else str(value).strip()
    label_lower = label.lower()

    if label in ['呋喃二甲醇', '2,5呋喃二甲醇'] or 'fdm' in label_lower:
        return '呋喃二甲醇'
    if label == '糠醇' or 'kc' in label_lower:
        return '糠醇'
    if label == '糠胺' or any(keyword in label_lower for keyword in ['ka', 'fge', 'fdo', 'fpdp']):
        return '糠胺'
    return ''


def infer_furan_type_from_row(row):
    cleaned_value = normalize_furan_type_label(row.get('furan_diol_clean', ''))
    if cleaned_value:
        return cleaned_value

    g_type = str(row.get('group_type', '')).lower()
    if 'furan_diol' not in g_type:
        return ''

    name = str(row.get('component_name', '')).strip()
    return normalize_furan_type_label(name)


def is_maleimide_component(row):
    name = str(row.get('component_name', '')).lower().strip()
    g_type = str(row.get('group_type', '')).lower().strip()
    maleimide_keywords = ['maleimide', 'bmi', 'tmi', 'dhpm', '1,8-bmi', '马来']
    return any(keyword in name or keyword in g_type for keyword in maleimide_keywords)


def is_furan_component(row):
    if is_maleimide_component(row):
        return False

    name = str(row.get('component_name', '')).lower().strip()
    g_type = str(row.get('group_type', '')).lower().strip()
    furan_keywords = ['furan', '呋喃', 'ka', 'kc', 'fge', 'fdo', 'fpdp', 'fhf', '糠胺', '糠醇']
    return any(keyword in name or keyword in g_type for keyword in furan_keywords)


def get_furan_type_code(row):
    if row['role'] != 'DA_monomer' or not is_furan_component(row):
        return 0
    furan_type = infer_furan_type_from_row(row)
    code_map = {
        '呋喃二甲醇': 1,
        '糠醇': 2,
        '糠胺': 3,
    }
    return code_map.get(furan_type, 0)


def get_iso_class(row):
    if row['role'] != 'Isocyanate':
        return np.nan
    name = str(row['component_name']).upper()
    aromatic = ['MDI', 'TDI', 'NDI', 'PPDI', 'XDI']
    if any(keyword in name for keyword in aromatic):
        return 1
    return 0


def get_polyol_class(row):
    if row['role'] != 'Polyol_Soft':
        return np.nan
    name = str(row['component_name']).upper()
    ester_keywords = ['PBA', 'PCL', 'PADG', 'CAPA', 'POLYESTER', 'PCDL', 'PLA', '酯']
    if any(keyword in name for keyword in ester_keywords):
        return 1
    return 0


def get_hard_symmetry(row):
    if row['role'] != 'Isocyanate':
        return np.nan
    name = str(row['component_name']).upper()
    symmetric = ['MDI', 'HDI', 'NDI', 'CDI', 'PPDI', 'CHDI']
    if any(keyword in name for keyword in symmetric):
        return 1.0
    return 0.0


def get_extender_class(row):
    if row['role'] != 'Chain_Extender':
        return 0.0

    name = str(row['component_name']).upper().strip()
    rigid_keywords = ['BQDO', 'HQEE', 'HER', 'MOCA', 'CHDM', 'ISOSORBIDE', '环', 'BENZ', 'CYCLO', 'PHE']
    if any(keyword in name for keyword in rigid_keywords):
        return 1.0
    return 0.0


def get_soft_cryst_flag(row):
    if row['role'] != 'Polyol_Soft':
        return 0.0
    name = str(row['component_name']).upper()
    cryst_types = ['PTMG', 'PCL', 'PBA', 'PEO', 'PEG', 'PLLA', 'PCDL']
    if any(keyword in name for keyword in cryst_types):
        mw = clean_numeric(row.get('M_i', 0))
        if 0 < mw < 600:
            return 0.0
        return 1.0
    return 0.0


def calculate_features(group):
    group = group.copy()
    group['mass_i'] = group['n_i'] * group['M_i']
    group['f_i'] = group['f_i'].fillna(2.0)

    total_mass = group['mass_i'].sum()
    if total_mass <= 1e-9:
        total_mass = 1.0

    m_polyol = group.loc[group['role'] == 'Polyol_Soft', 'mass_i'].sum()
    soft_segment_content = m_polyol / total_mass

    furan_mask = group.apply(is_furan_component, axis=1)
    maleimide_mask = group.apply(is_maleimide_component, axis=1)

    m_furan = group.loc[furan_mask, 'mass_i'].sum()
    furan_content = m_furan / total_mass

    m_maleimide = group.loc[maleimide_mask, 'mass_i'].sum()
    maleimide_content = m_maleimide / total_mass

    cross_mask = (group['role'] == 'Chain_Extender') & (group['f_i'] > 2.0)
    m_add_cross = group.loc[cross_mask, 'mass_i'].sum()
    extra_crosslinker_content = m_add_cross / total_mass

    m_iso = group.loc[group['role'] == 'Isocyanate', 'mass_i'].sum()
    linear_ext_mask = (group['role'] == 'Chain_Extender') & (group['f_i'] <= 2.0)
    m_linear_ext = group.loc[linear_ext_mask, 'mass_i'].sum()
    hard_segment_content = (m_iso + m_linear_ext) / total_mass

    n_nco = (group.loc[group['role'] == 'Isocyanate', 'n_i'] * group.loc[group['role'] == 'Isocyanate', 'f_i']).sum()
    n_active = (group.loc[group['role'] != 'Isocyanate', 'n_i'] * group.loc[group['role'] != 'Isocyanate', 'f_i']).sum()
    nco_oh_ratio = n_nco / n_active if n_active > 0 else 0.0

    furan_type_series = group.apply(get_furan_type_code, axis=1)
    furan_type = furan_type_series.max()

    network_topology = 1 if group['DA_class'].max() == 1 else 0
    add_crosslinker_flag = 1 if extra_crosslinker_content > 0 else 0

    poly_idx = group.loc[group['role'] == 'Polyol_Soft', 'mass_i'].idxmax() if m_polyol > 0 else None
    polyol_class = get_polyol_class(group.loc[poly_idx]) if poly_idx is not None else 0

    iso_idx = group.loc[group['role'] == 'Isocyanate', 'mass_i'].idxmax() if m_iso > 0 else None
    iso_class = get_iso_class(group.loc[iso_idx]) if iso_idx is not None else 0

    ext_scores = group.apply(get_extender_class, axis=1)
    extender_class = 1 if ext_scores.max() > 0 else 0

    if m_polyol > 0:
        soft_mw = np.average(
            group.loc[group['role'] == 'Polyol_Soft', 'M_i'],
            weights=group.loc[group['role'] == 'Polyol_Soft', 'mass_i']
        )
    else:
        soft_mw = 0.0

    if m_iso > 0:
        sym_scores = group.apply(get_hard_symmetry, axis=1)
        w_sym = np.average(
            sym_scores[group['role'] == 'Isocyanate'],
            weights=group.loc[group['role'] == 'Isocyanate', 'mass_i']
        )
        hard_symmetry = 1 if w_sym > 0.5 else 0
    else:
        hard_symmetry = 0

    if m_polyol > 0:
        cryst_scores = group.apply(get_soft_cryst_flag, axis=1)
        soft_cryst = 1 if cryst_scores.max() > 0 else 0
    else:
        soft_cryst = 0

    synergy_feature = soft_segment_content * soft_cryst * hard_symmetry
    constraint_factor = hard_segment_content / soft_mw if soft_mw > 100 else 0.0

    def get_val(col):
        values = group[col].dropna()
        return values.iloc[0] if not values.empty else np.nan

    return pd.Series({
        'Soft_Segment_Content': soft_segment_content,
        'Furan_Content': furan_content,
        'Maleimide_Content': maleimide_content,
        'Hard_Segment_Content': hard_segment_content,
        'NCO_OH_Ratio': nco_oh_ratio,
        'Extra_Crosslinker_Content': extra_crosslinker_content,
        'Furan_Type': furan_type,
        'Network_Topology': network_topology,
        'Add_Crosslinker_Flag': add_crosslinker_flag,
        'Polyol_Class': polyol_class,
        'Iso_Class': iso_class,
        'Extender_Class': extender_class,
        'Soft_Mw': soft_mw,
        'Constraint_Factor': constraint_factor,
        'Hard_Symmetry': hard_symmetry,
        'Soft_Cryst': soft_cryst,
        'Synergy_Feature': synergy_feature,
        'poly_tem': get_val('poly_tem'),
        'strain_rate': get_val('strain_rate'),
        'healing_temperature': get_val('healing_temperature'),
        'healing_time': get_val('healing_time'),
        'tensile_strength': get_val('tensile_strength'),
        'elongation': get_val('elongation'),
        'healing_eff': get_val('healing_eff')
    })


def load_raw_table(file_path):
    df = pd.read_csv(file_path)

    col_map = {
        '名称sample_id': 'sample_id',
        '含DA和不含(DA_strategy)': 'DA_strategy',
        '合成温度(poly_tem)': 'poly_tem',
        '交联/线性(cross_class)': 'cross_class',
        '组分名称component_name': 'component_name',
        '组分摩尔质量M_i（g/mol）': 'M_i',
        '组分角色group_type（Isocyanate Hydroxyl Amine': 'group_type',
        '组分角色group_type（Isocyanate Hydroxyl Amine ': 'group_type',
        'DA位置（DA_class）': 'DA_class',
        '组分摩尔用量（n_i）': 'n_i',
        '拉伸测试条件(strain_rate mm/min)': 'strain_rate',
        '自愈合温度（healing_temperature ℃)': 'healing_temperature',
        '自愈合时间（healing time h）': 'healing_time',
        '输出：原拉伸强度MPa（tensile_strength）': 'tensile_strength',
        '输出：初始拉伸率%(elongation)': 'elongation',
        '输出：自愈合率%（healing_eff）': 'healing_eff',
    }

    df.columns = [str(col).strip() for col in df.columns]
    df = df.rename(columns={key.strip(): value for key, value in col_map.items()})
    df['sample_id'] = df['sample_id'].ffill()

    cond_cols = [
        'DA_strategy', 'DA_class', 'poly_tem', 'strain_rate',
        'healing_temperature', 'healing_time',
        'tensile_strength', 'elongation', 'healing_eff'
    ]
    for col in cond_cols:
        if col in df.columns:
            df[col] = df.groupby('sample_id')[col].transform(lambda x: x.ffill().bfill())

    df = df.dropna(subset=['component_name']).copy()

    num_cols = [
        'M_i', 'n_i', 'poly_tem', 'strain_rate',
        'healing_temperature', 'healing_time',
        'tensile_strength', 'elongation', 'healing_eff'
    ]
    for col in num_cols:
        if col in df.columns:
            distinct = col in ['tensile_strength', 'elongation', 'healing_eff']
            df[col] = df[col].apply(lambda x: clean_numeric(x, distinct_nan=distinct))

    if 'DA_class' in df.columns:
        df['DA_class'] = df['DA_class'].apply(parse_side_chain_flag)
    else:
        df['DA_class'] = 0

    if 'DA_strategy' in df.columns:
        df['DA_strategy'] = df['DA_strategy'].apply(parse_da_strategy)

    if 'furan_diol_clean' not in df.columns:
        df['furan_diol_clean'] = ''

    df['role'] = df.apply(classify_role, axis=1)
    df['f_i'] = df.apply(infer_functionality, axis=1)
    return df


def build_feature_dataframe(file_path):
    df = load_raw_table(file_path)
    feature_df = df.groupby('sample_id', sort=False).apply(calculate_features).reset_index()
    return feature_df.round(4)


if __name__ == '__main__':
    print(f'正在处理文件: {INPUT_FILE.name}')
    feature_df = build_feature_dataframe(INPUT_FILE)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    feature_df.to_csv(OUTPUT_FILE, index=False, encoding='utf-8-sig')
    print('-' * 30)
    print(f'特征工程完成！样本数: {len(feature_df)}')
    print(f'生成的特征数: {feature_df.shape[1]}')
    print(f'文件已保存至: {OUTPUT_FILE}')